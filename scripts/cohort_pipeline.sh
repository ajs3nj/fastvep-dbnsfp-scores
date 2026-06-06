#!/usr/bin/env bash
# cohort_pipeline.sh -- end-to-end cohort variant annotation + tiering.
#
# Reads a manifest TSV of (sample_id, vcf_path) pairs and runs:
#   Stage 1: per-sample fastvep annotate + bcftools genotype extraction (parallel)
#   Stage 2: cohort merge -- dedupe by (chrom,pos,ref,alt,transcript) AND
#                            count distinct samples each variant appears in
#                            (added as `n_samples_annotated` column)
#   Stage 3: pre-filter -- drop common / no-signal variants
#   Stage 4: single tier_variants.R run on the cohort
#
# Outputs a single cohort.variants.tsv (one row per unique variant with tier,
# modifier flag, n_samples_annotated, cohort AC/AN/n_carriers, etc.) and
# cohort.genes.tsv (per-gene rollup with sample-burden columns).
#
# Manifest format (TSV, header required):
#   sample_id<TAB>vcf_path
#   S001<TAB>/data/vcfs/S001.vcf.gz
#   S002<TAB>/data/vcfs/S002.vcf.gz
#   ...

set -euo pipefail

# ============================== usage ========================================
usage() {
  cat >&2 <<EOF
Usage: $(basename "$0") --manifest FILE --out-dir DIR --data-dir DIR --sa-dir DIR [options]

Required:
  --manifest FILE        TSV with header 'sample_id<TAB>vcf_path'
  --out-dir DIR          Output directory (created if missing)
  --data-dir DIR         Holds GFF3, FASTA, transcript cache
  --sa-dir DIR           Holds built .osa/.osi/.oga files

Optional:
  --threads N            Parallel sample-annotation jobs (default: nproc/2)
  --gff3 FILE            Override default GFF3 path
                         (default: \$DATA_DIR/Homo_sapiens.GRCh38.115.gff3)
  --fasta FILE           Override default FASTA path
                         (default: \$DATA_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa)
  --cache FILE           Override default transcript cache path
                         (default: \$DATA_DIR/grch38_115.fvcache)
  --modifier-genes FILE  Gene list for modifier-flag scope
                         (default: <repo>/R/nf_modifier_genes.txt)
  --tier-script PATH     Path to tier_variants.R
                         (default: <repo>/R/tier_variants.R)
  --af-max FLOAT         Tiering rarity gate (default: 1e-4)
  --filter-af FLOAT      Pre-filter: drop variants with gnomAD popmax > this
                         (default: 0.01)
  --skip-annotate        Skip Stage 1; reuse existing \$OUT/per_sample/*
  --skip-merge           Skip Stage 2; reuse existing \$OUT/cohort.{annotated,genotypes}.*
  --include-sex-chroms   Keep chrX/chrY variants (default: drop them so cohort
                         tiering doesn't have to reason about hemizygous dosage
                         or split by sex). Mitochondrial (chrM/MT) is always kept.
  -h, --help             This message

Manifest example: see scripts/cohort.manifest.example.tsv
EOF
  exit 1
}

# ============================== args =========================================
MANIFEST=""; OUT_DIR=""; DATA_DIR=""; SA_DIR=""
THREADS=$(( $(nproc 2>/dev/null || echo 4) / 2 ))
GFF3=""; FASTA=""; CACHE=""
MODIFIER_GENES=""; TIER_SCRIPT=""; CSQ_CONVERTER=""
PYTHON_BIN="${PYTHON_BIN:-python3}"
AF_MAX="1e-4"; FILTER_AF="0.01"
SKIP_ANNOTATE=0; SKIP_MERGE=0
INCLUDE_SEX_CHROMS=0   # default: drop X/Y to avoid hemizygous / dosage complexity

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)        MANIFEST="$2";       shift 2 ;;
    --out-dir)         OUT_DIR="$2";        shift 2 ;;
    --data-dir)        DATA_DIR="$2";       shift 2 ;;
    --sa-dir)          SA_DIR="$2";         shift 2 ;;
    --threads)         THREADS="$2";        shift 2 ;;
    --gff3)            GFF3="$2";           shift 2 ;;
    --fasta)           FASTA="$2";          shift 2 ;;
    --cache)           CACHE="$2";          shift 2 ;;
    --modifier-genes)  MODIFIER_GENES="$2"; shift 2 ;;
    --tier-script)     TIER_SCRIPT="$2";    shift 2 ;;
    --af-max)          AF_MAX="$2";         shift 2 ;;
    --filter-af)       FILTER_AF="$2";      shift 2 ;;
    --skip-annotate)   SKIP_ANNOTATE=1;     shift ;;
    --skip-merge)      SKIP_MERGE=1;        shift ;;
    --include-sex-chroms) INCLUDE_SEX_CHROMS=1; shift ;;
    -h|--help)         usage ;;
    *) echo "Unknown arg: $1" >&2; usage ;;
  esac
done

[[ -z "$MANIFEST" || -z "$OUT_DIR" || -z "$DATA_DIR" || -z "$SA_DIR" ]] && usage

# Resolve defaults relative to the repo this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
[[ -z "$GFF3"           ]] && GFF3="$DATA_DIR/Homo_sapiens.GRCh38.115.gff3"
[[ -z "$FASTA"          ]] && FASTA="$DATA_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
[[ -z "$CACHE"          ]] && CACHE="$DATA_DIR/grch38_115.fvcache"
[[ -z "$TIER_SCRIPT"    ]] && TIER_SCRIPT="$REPO_DIR/R/tier_variants.R"
[[ -z "$MODIFIER_GENES" ]] && MODIFIER_GENES="$REPO_DIR/R/nf_modifier_genes.txt"
[[ -z "$CSQ_CONVERTER"  ]] && CSQ_CONVERTER="$REPO_DIR/scripts/csq_to_wide_tab.py"

# ============================== helpers ======================================
log()     { printf '[%(%H:%M:%S)T] %s\n' -1 "$*"; }
section() { echo; log "=== $* ==="; }
die()     { echo "ERROR: $*" >&2; exit 1; }

# ============================== validation ===================================
section "Validation"
for f in "$MANIFEST" "$FASTA" "$TIER_SCRIPT" "$MODIFIER_GENES" "$CSQ_CONVERTER"; do
  [[ -e "$f" ]] || die "file not found: $f"
done
# Need either the binary transcript cache OR the GFF3 (cache wins if both).
if [[ -e "$CACHE" ]]; then
  log "Using transcript cache: $CACHE"
  GFF3=""   # signal: don't pass --gff3 to fastvep
elif [[ -e "$GFF3" ]]; then
  log "Using GFF3 (no cache found): $GFF3"
  CACHE=""
else
  die "neither transcript cache nor GFF3 found (looked for $CACHE and $GFF3)"
fi
[[ -d "$SA_DIR" ]] || die "SA dir not found: $SA_DIR"
command -v fastvep      >/dev/null || die "fastvep not on PATH"
command -v bcftools     >/dev/null || die "bcftools not on PATH"
command -v Rscript      >/dev/null || die "Rscript not on PATH"
command -v "$PYTHON_BIN" >/dev/null || die "$PYTHON_BIN not on PATH"

mkdir -p "$OUT_DIR/per_sample" "$OUT_DIR/logs"

# Validate manifest header. First two columns must be sample_id and vcf_path.
# Any additional columns are passed through as sample metadata.
header=$(head -1 "$MANIFEST")
first_two=$(echo "$header" | awk -F'\t' '{print $1"\t"$2}')
expected=$'sample_id\tvcf_path'
[[ "$first_two" == "$expected" ]] \
  || die "manifest first two columns must be 'sample_id<TAB>vcf_path'; got: '$first_two'"

# Pass-through sample metadata (cols 3+ if present) for downstream joins.
cp "$MANIFEST" "$OUT_DIR/sample_metadata.tsv"

n_samples=$(tail -n+2 "$MANIFEST" | grep -c .)
[[ "$n_samples" -gt 0 ]] || die "manifest has no sample rows"
log "Manifest OK: $n_samples samples"

# Validate VCF paths
missing=$(tail -n+2 "$MANIFEST" | awk -F'\t' '{ if (system("test -e " $2 ) != 0) print $1": "$2 }')
[[ -n "$missing" ]] && { echo "ERROR: missing VCFs:" >&2; echo "$missing" >&2; exit 1; }

log "Threads: $THREADS"
log "Pre-filter AF threshold: $FILTER_AF (gnomAD popmax)"
log "Tiering rarity gate:     $AF_MAX"

# ============================== Stage 1: per-sample =========================
if [[ "$SKIP_ANNOTATE" == "0" ]]; then
  section "Stage 1: annotating $n_samples samples (parallel x $THREADS)"

  process_one() {
    local sample_id="$1" vcf_path="$2"
    local out_vcf="$OUT_DIR/per_sample/${sample_id}.annotated.vcf"
    local out_tab="$OUT_DIR/per_sample/${sample_id}.annotated.tab"
    local out_gt="$OUT_DIR/per_sample/${sample_id}.genotypes.tsv"
    local log_file="$OUT_DIR/logs/${sample_id}.log"

    if [[ -s "$out_tab" && -s "$out_gt" ]]; then
      echo "[$sample_id] cached"
      return 0
    fi

    {
      echo "=== $(date) $sample_id annotating (VCF intermediate) ==="
      # Use VCF output so the full CSQ + FV_* INFO fields are preserved.
      # Pass either --transcript-cache or --gff3, whichever was validated upstream.
      backbone_args=()
      if [[ -n "$CACHE" ]]; then
        backbone_args=(--transcript-cache "$CACHE")
      else
        backbone_args=(--gff3 "$GFF3")
      fi
      fastvep annotate -i "$vcf_path" -o "$out_vcf" \
        "${backbone_args[@]}" \
        --fasta "$FASTA" \
        --sa-dir "$SA_DIR" \
        --output-format vcf \
        --hgvs --canonical

      echo "=== $(date) $sample_id converting CSQ -> wide tab ==="
      "$PYTHON_BIN" "$CSQ_CONVERTER" \
        --input "$out_vcf" \
        --output "$out_tab" \
        --sample "$sample_id"

      # Discard the VCF intermediate -- the wide tab carries everything we need
      # and the VCF is ~5-10x larger.
      rm -f "$out_vcf"

      echo "=== $(date) $sample_id extracting genotypes ==="
      (
        printf 'chrom\tpos\tref\talt\tsample\tgt\n'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$vcf_path" | \
          awk -v s="$sample_id" 'BEGIN{OFS="\t"} { print $1, $2, $3, $4, s, $5 }'
      ) > "$out_gt"

      echo "=== $(date) $sample_id done ==="
    } > "$log_file" 2>&1

    echo "[$sample_id] done"
  }
  export -f process_one
  export OUT_DIR CACHE GFF3 FASTA SA_DIR PYTHON_BIN CSQ_CONVERTER

  # Process in parallel. Use GNU parallel if available (cleaner tab handling),
  # else fall back to xargs -n 2 with a space-delimited list (we know vcf paths
  # in our manifest don't contain spaces; if they ever do, switch to parallel).
  if command -v parallel >/dev/null 2>&1; then
    tail -n+2 "$MANIFEST" | awk -F'\t' '{ print $1"\t"$2 }' | \
      parallel -j "$THREADS" --colsep '\t' process_one {1} {2}
  else
    tail -n+2 "$MANIFEST" | awk -F'\t' '{ print $1, $2 }' | \
      xargs -P "$THREADS" -n 2 bash -c 'process_one "$1" "$2"' _
  fi

  log "Stage 1 complete"
else
  log "Stage 1 skipped (--skip-annotate)"
fi

# ============================== Stage 2: merge ==============================
if [[ "$SKIP_MERGE" == "0" ]]; then
  section "Stage 2: cohort merge"

  # 2ab. Per-sample dedupe + cohort-wide count in ONE pass, no intermediate file.
  # The `seen` hash tracks per-sample uniqueness and is cleared between files via
  # `delete`. Only the global `count` hash persists (~10-30M variants ~ 1-1.5 GB).
  # The old write-to-disk-then-aggregate version produced a ~30+ GB intermediate
  # at cohort scale (one row per (variant, sample)) and was disk-bound.
  log "  2ab) counting samples per variant (single-pass)..."
  awk -F'\t' '
    FNR == 1 {
      for (k in seen) delete seen[k]   # reset per-sample dedupe set
      next                              # skip header
    }
    {
      key = $1"\t"$2"\t"$3"\t"$4
      if (!(key in seen)) {
        seen[key] = 1
        count[key]++
      }
    }
    END { for (k in count) print k"\t"count[k] }
  ' "$OUT_DIR/per_sample/"*.annotated.tab > "$OUT_DIR/_variant_counts.tsv"

  # 2c. Dedupe annotated rows by (chrom,pos,ref,alt,transcript).
  log "  2c) deduplicating annotated rows..."
  first_tab=$(ls "$OUT_DIR/per_sample/"*.annotated.tab | head -1)
  head -1 "$first_tab" > "$OUT_DIR/_cohort.uniq.tab"
  # Determine transcript column index from header (default 7 in fastVEP tab)
  TCOL=$(head -1 "$first_tab" | awk -F'\t' '{ for (i=1;i<=NF;i++) if ($i=="transcript") print i }')
  TCOL=${TCOL:-7}
  for f in "$OUT_DIR/per_sample/"*.annotated.tab; do
    tail -n+2 "$f"
  done | awk -F'\t' -v tc="$TCOL" '!seen[$1"\t"$2"\t"$3"\t"$4"\t"$tc]++' \
       >> "$OUT_DIR/_cohort.uniq.tab"

  # 2d. Join sample counts onto annotated rows as n_samples_annotated.
  log "  2d) attaching n_samples_annotated column..."
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR == FNR { count[$1"\t"$2"\t"$3"\t"$4] = $5; next }
    FNR == 1  { print $0, "n_samples_annotated"; next }
              { key = $1"\t"$2"\t"$3"\t"$4
                print $0, (key in count ? count[key] : 0) }
  ' "$OUT_DIR/_variant_counts.tsv" "$OUT_DIR/_cohort.uniq.tab" \
    > "$OUT_DIR/cohort.annotated.tab"

  # 2e. Concatenate per-sample genotype TSVs (no dedupe needed).
  log "  2e) concatenating genotype TSVs..."
  first_gt=$(ls "$OUT_DIR/per_sample/"*.genotypes.tsv | head -1)
  head -1 "$first_gt" > "$OUT_DIR/cohort.genotypes.tsv"
  for f in "$OUT_DIR/per_sample/"*.genotypes.tsv; do
    tail -n+2 "$f"
  done >> "$OUT_DIR/cohort.genotypes.tsv"

  rm -f "$OUT_DIR/_sample_variants.tsv" \
        "$OUT_DIR/_variant_counts.tsv" \
        "$OUT_DIR/_cohort.uniq.tab"

  n_unique=$(($(wc -l < "$OUT_DIR/cohort.annotated.tab") - 1))
  n_gt=$(($(wc -l < "$OUT_DIR/cohort.genotypes.tsv") - 1))
  log "Merged: $n_unique unique variant-transcript rows; $n_gt genotype rows"
else
  log "Stage 2 skipped (--skip-merge)"
fi

# ============================== Stage 3: pre-filter =========================
section "Stage 3: pre-filter"
[[ "$INCLUDE_SEX_CHROMS" == "0" ]] && log "  excluding chrX/chrY (override: --include-sex-chroms)"

awk -F'\t' -v thresh="$FILTER_AF" -v drop_xy="$((1 - INCLUDE_SEX_CHROMS))" '
  NR == 1 {
    for (i=1; i<=NF; i++) col[$i] = i
    print
    next
  }
  {
    # Drop sex chromosomes unless --include-sex-chroms is set.
    if (drop_xy) {
      chrom = $col["chrom"]
      if (chrom == "X" || chrom == "Y" || chrom == "chrX" || chrom == "chrY") next
    }

    keep = 0
    if (col["impact"] && ($col["impact"] == "HIGH" || $col["impact"] == "MODERATE")) keep = 1
    else if (col["spliceai_ds_max"] && $col["spliceai_ds_max"] != "" && $col["spliceai_ds_max"] != "." && $col["spliceai_ds_max"]+0 >= 0.2) keep = 1
    else if (col["clin_sig"] && $col["clin_sig"] != "" && $col["clin_sig"] != ".") keep = 1
    else if (col["gnomad_popmax_af"]) {
      af = $col["gnomad_popmax_af"]
      if (af == "" || af == "." || af+0 <= thresh) keep = 1
    }
    if (keep) print
  }
' "$OUT_DIR/cohort.annotated.tab" > "$OUT_DIR/cohort.filtered.tab"

n_pre=$(($(wc -l < "$OUT_DIR/cohort.annotated.tab") - 1))
n_post=$(($(wc -l < "$OUT_DIR/cohort.filtered.tab") - 1))
[[ "$n_pre" -gt 0 ]] && pct=$(awk "BEGIN{printf \"%.1f\", $n_post*100/$n_pre}") || pct="0.0"
log "Pre-filter retained $n_post / $n_pre rows (${pct}%)"

# ============================== Stage 4: tier ===============================
section "Stage 4: tiering"
Rscript "$TIER_SCRIPT" \
  --input         "$OUT_DIR/cohort.filtered.tab" \
  --genotypes     "$OUT_DIR/cohort.genotypes.tsv" \
  --modifier-genes "$MODIFIER_GENES" \
  --af-max        "$AF_MAX" \
  --out-prefix    "$OUT_DIR/cohort"

# ============================== Summary =====================================
section "Summary"
log "Cohort outputs:"
log "  variants table:  $OUT_DIR/cohort.variants.tsv ($(( $(wc -l < "$OUT_DIR/cohort.variants.tsv") - 1 )) rows)"
log "  gene table:      $OUT_DIR/cohort.genes.tsv ($(( $(wc -l < "$OUT_DIR/cohort.genes.tsv") - 1 )) rows)"
log ""
log "Tier breakdown (per-variant):"
awk -F'\t' 'NR == 1 {
              for (i=1; i<=NF; i++) if ($i == "tier") tc=i
              for (i=1; i<=NF; i++) if ($i == "modifier_candidate") mc=i
              next
            }
            { t[$tc]++; if (mc && $mc == "TRUE") m++ }
            END {
              for (k=1; k<=5; k++) if (k in t) printf "  Tier %d: %d\n", k, t[k]
              if (m) printf "  modifier_candidate: %d\n", m
            }' "$OUT_DIR/cohort.variants.tsv"

log ""
log "Done. Top-tier candidates: head -20 $OUT_DIR/cohort.variants.tsv"
