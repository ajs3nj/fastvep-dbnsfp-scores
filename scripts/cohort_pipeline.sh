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
  --filter-allow-missing-popmax
                         (Stage 3) Treat variants with missing gnomAD popmax as
                         rare-enough-to-keep. Default OFF: missing popmax does
                         not bypass the rarity gate, since at cohort scale that
                         rule keeps ~90% of non-coding variants and OOMs R.
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
ALLOW_MISSING_POPMAX=0 # default: variants with no popmax need other signal to pass filter

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
    --filter-allow-missing-popmax) ALLOW_MISSING_POPMAX=1; shift ;;
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

# Point GNU parallel (and any other temp-file-creators) at a directory on the
# output volume rather than /tmp, which is often tiny on containerized hosts
# (we've seen 64M /dev/shm and similarly constrained /tmp). With 8 parallel
# fastvep workers, even modest per-job output buffering fills /tmp quickly.
export TMPDIR="$OUT_DIR/_parallel_tmp"
mkdir -p "$TMPDIR"
log "TMPDIR=$TMPDIR  ($(df -h "$TMPDIR" | awk 'NR==2{print $4" free"}'))"

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

# --- FASTA .fai index check -------------------------------------------------
# Without .fai, fastvep loads the entire ~3 GB FASTA into RAM per worker
# instead of memory-mapping, causing 30-90s startup hits and per-variant
# slowdown. Auto-generate if samtools is available; otherwise fail fast with
# the regen command so the user isn't surprised by silently-slow runs.
if [[ ! -e "${FASTA}.fai" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    log "FASTA .fai index missing; generating with samtools faidx (one-time, ~30s)..."
    samtools faidx "$FASTA" \
      || die "samtools faidx failed on $FASTA -- check file integrity"
  else
    die "FASTA .fai index missing at ${FASTA}.fai and samtools not on PATH.
       Generate the index once with:    samtools faidx $FASTA
       Without it, per-worker FASTA load adds 30-90s of overhead per sample
       and total cohort runtime can multiply 5-10x. Install samtools via
       'conda install -c bioconda samtools' and re-run."
  fi
fi

# --- SA database contents check ---------------------------------------------
# List which .osa/.osa.idx/.oga files are in --sa-dir, warn loudly if any of
# the tiering-critical ones are missing. Tiering can still produce output
# without these, but the result will be largely uninformative.
log "Supplementary annotation databases in $SA_DIR:"
sa_found=0
for expected_base in dbnsfp clinvar gnomad_genes spliceai; do
  matches=()
  for ext in osa osa.idx oga; do
    for f in "$SA_DIR/${expected_base}.${ext}" "$SA_DIR/${expected_base}.${ext}".idx; do
      [[ -e "$f" ]] && matches+=("$f")
    done
  done
  if [[ ${#matches[@]} -gt 0 ]]; then
    primary="${matches[0]}"
    log "  ✓ $expected_base → $(basename "$primary") ($(du -h "$primary" | cut -f1))"
    sa_found=$((sa_found + 1))
  else
    case "$expected_base" in
      dbnsfp|clinvar|gnomad_genes)
        log "  ⚠ $expected_base NOT FOUND  -- tiering will lack signal for this source"
        ;;
      spliceai)
        log "  ⓘ $expected_base not loaded   (acceptable; see docs/noncoding_v2_plan.md)"
        ;;
    esac
  fi
done
[[ $sa_found -eq 0 ]] && die "no recognized SA databases found in $SA_DIR -- tiering would be uninformative; aborting"

# --- Chromosome-naming check ------------------------------------------------
# Compare the chromosome naming convention of the first VCF against the FASTA.
# If they don't match, fastvep produces all-intergenic CSQ output silently and
# the entire cohort run is wasted. Fail fast with a clear remediation message.
#
# Skip entirely when --skip-annotate is set: the check exists to catch a
# problem that only manifests during Stage 1 annotation, and when we're
# re-running just Stages 2-5 (e.g. via run_cohort_stages.sh with a synthesized
# v1 manifest) the "first VCF" in the manifest may be a dummy path that
# doesn't reflect the real source data.
if [[ "$SKIP_ANNOTATE" == "1" ]]; then
  log "  ⓘ skipping chromosome-naming check (--skip-annotate is set)"
else
first_vcf=$(tail -n+2 "$MANIFEST" | awk -F'\t' 'NR==1 {print $2; exit}')
vcf_chrom=$(zcat -f "$first_vcf" 2>/dev/null | awk '!/^#/ {print $1; exit}')
fasta_chrom=$(grep '^>' "$FASTA" 2>/dev/null | head -1 | sed 's/^>//' | awk '{print $1}')

vcf_chr=$(echo "$vcf_chrom" | grep -c '^chr' || true)
fasta_chr=$(echo "$fasta_chrom" | grep -c '^chr' || true)

if [[ -z "$vcf_chrom" || -z "$fasta_chrom" ]]; then
  log "  ⚠ couldn't determine chromosome naming (vcf='$vcf_chrom' fasta='$fasta_chrom'); proceeding"
elif [[ "$vcf_chr" != "$fasta_chr" ]]; then
  die "CHROMOSOME NAMING MISMATCH detected!
       First VCF chrom:  '$vcf_chrom'  ($first_vcf)
       FASTA chrom:      '$fasta_chrom'  ($FASTA)
       fastvep will produce ALL-INTERGENIC annotation (no gene/transcript
       match for any variant), wasting the entire cohort run silently.

       Fix one or the other before re-running:
       a) Rename VCFs to match FASTA convention (Ensembl-style, no chr prefix):
            for vcf in /path/to/*.vcf.gz; do
              bcftools annotate --rename-chrs /tmp/chr_to_ensembl.txt \\
                -Oz -o renamed/\$(basename \$vcf) \$vcf
            done
       b) Rebuild the cache + use a chr-prefixed FASTA throughout.

       See docs/pipeline_methods.md \xc2\xa72.1 for chromosome-naming requirements."
else
  log "  ✓ chromosome naming consistent (VCF='$vcf_chrom' FASTA='$fasta_chrom')"
fi
fi   # end: skip chrom-naming check when --skip-annotate

log "Threads: $THREADS"
log "Pre-filter AF threshold: $FILTER_AF (gnomAD popmax)"
log "Tiering rarity gate:     $AF_MAX"

# --- Provenance log ---------------------------------------------------------
# Records git commit, tool versions, reference data paths + mtimes, all flags
# used. Written before Stage 1 starts so it survives any subsequent failure.
PROVENANCE="$OUT_DIR/run.provenance.tsv"
{
  printf "key\tvalue\n"
  printf "timestamp_utc\t%s\n"        "$(date -u +%FT%TZ)"
  printf "hostname\t%s\n"             "$(hostname)"
  printf "git_commit\t%s\n"           "$(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
  printf "git_branch\t%s\n"           "$(git -C "$REPO_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
  printf "git_dirty\t%s\n"            "$(git -C "$REPO_DIR" diff --quiet 2>/dev/null && echo no || echo yes)"
  printf "fastvep_version\t%s\n"      "$(fastvep --version 2>&1 | head -1)"
  printf "bcftools_version\t%s\n"     "$(bcftools --version 2>&1 | head -1)"
  printf "rscript_version\t%s\n"      "$(Rscript --version 2>&1 | head -1)"
  printf "python_version\t%s\n"       "$($PYTHON_BIN --version 2>&1 | head -1)"
  printf "samtools_version\t%s\n"     "$(samtools --version 2>&1 | head -1 || echo missing)"
  printf "manifest\t%s\n"             "$MANIFEST"
  printf "n_samples\t%s\n"            "$n_samples"
  printf "out_dir\t%s\n"              "$OUT_DIR"
  printf "data_dir\t%s\n"             "$DATA_DIR"
  printf "fasta_path\t%s\n"           "$FASTA"
  printf "fasta_mtime\t%s\n"          "$(stat -c %y "$FASTA" 2>/dev/null || stat -f '%Sm' "$FASTA" 2>/dev/null)"
  printf "fasta_size\t%s\n"           "$(du -h "$FASTA" | cut -f1)"
  printf "cache_path\t%s\n"           "${CACHE:-(unused)}"
  printf "gff3_path\t%s\n"            "${GFF3:-(unused)}"
  printf "sa_dir\t%s\n"               "$SA_DIR"
  for osa in "$SA_DIR"/*.osa "$SA_DIR"/*.oga; do
    [[ -e "$osa" ]] || continue
    printf "sa_file:%s\t%s | mtime=%s\n" \
      "$(basename "$osa")" \
      "$(du -h "$osa" | cut -f1)" \
      "$(stat -c %y "$osa" 2>/dev/null || stat -f '%Sm' "$osa" 2>/dev/null)"
  done
  printf "modifier_genes_file\t%s\n"  "$MODIFIER_GENES"
  printf "modifier_genes_count\t%s\n" "$(grep -vc '^#\|^$' "$MODIFIER_GENES" 2>/dev/null || echo unknown)"
  printf "filter_af\t%s\n"            "$FILTER_AF"
  printf "af_max_tier_gate\t%s\n"     "$AF_MAX"
  printf "include_sex_chroms\t%s\n"   "$INCLUDE_SEX_CHROMS"
  printf "allow_missing_popmax\t%s\n" "$ALLOW_MISSING_POPMAX"
  printf "threads\t%s\n"              "$THREADS"
  printf "skip_annotate\t%s\n"        "$SKIP_ANNOTATE"
  printf "skip_merge\t%s\n"           "$SKIP_MERGE"
} > "$PROVENANCE"
log "Wrote provenance log: $PROVENANCE"

# ============================== Stage 1: per-sample =========================
if [[ "$SKIP_ANNOTATE" == "0" ]]; then
  section "Stage 1: annotating $n_samples samples (parallel x $THREADS)"

  process_one() {
    local sample_id="$1" vcf_path="$2"
    local out_vcf="$OUT_DIR/per_sample/${sample_id}.annotated.vcf"
    # Write gzipped per-sample tabs (~80% smaller; downstream awks read via zcat
    # globbing on both .tab and .tab.gz). Resume check below accepts either.
    local out_tab="$OUT_DIR/per_sample/${sample_id}.annotated.tab.gz"
    local out_tab_uncompressed="$OUT_DIR/per_sample/${sample_id}.annotated.tab"
    local out_gt="$OUT_DIR/per_sample/${sample_id}.genotypes.tsv"
    local log_file="$OUT_DIR/logs/${sample_id}.log"

    # Accept either compressed (new) or uncompressed (already-done from a prior run).
    if [[ ( -s "$out_tab" || -s "$out_tab_uncompressed" ) && -s "$out_gt" ]]; then
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

  # Gather the per-sample annotated.tab files. Glob both .tab and .tab.gz so a
  # cohort with mixed compression (older samples uncompressed, newer ones gzipped)
  # works through Stage 2 transparently.
  shopt -s nullglob
  ann_files=("$OUT_DIR/per_sample/"*.annotated.tab "$OUT_DIR/per_sample/"*.annotated.tab.gz)
  shopt -u nullglob
  [[ ${#ann_files[@]} -gt 0 ]] || die "no annotated.tab(.gz) files found in $OUT_DIR/per_sample/"
  first_tab="${ann_files[0]}"

  # Helper to stream a file's contents (handles both plain and gzipped).
  stream_tab() {
    case "$1" in
      *.gz) zcat "$1" ;;
      *)    cat  "$1" ;;
    esac
  }
  export -f stream_tab 2>/dev/null || true

  # 2ab. Per-sample dedupe + cohort-wide count in ONE pass, no intermediate file.
  # The `seen` hash tracks per-sample uniqueness and is cleared between files via
  # `delete`. Only the global `count` hash persists (~10-30M variants ~ 1-1.5 GB).
  # Resumable: skip if a non-empty _variant_counts.tsv already exists.
  if [[ -s "$OUT_DIR/_variant_counts.tsv" ]]; then
    log "  2ab) skipping (existing _variant_counts.tsv: $(du -h "$OUT_DIR/_variant_counts.tsv" | cut -f1))"
  else
    log "  2ab) counting samples per variant (single-pass, gzip-aware)..."
    # Concat all per-sample files into one stream, marking file boundaries with
    # a sentinel line so awk can reset the per-sample dedupe hash between files.
    # Works for both plain .tab and gzipped .tab.gz inputs.
    (
      for f in "${ann_files[@]}"; do
        printf '@@FILE@@\n'
        case "$f" in
          *.gz) zcat "$f" ;;
          *)    cat  "$f" ;;
        esac
      done
    ) | awk -F'\t' '
        /^@@FILE@@$/ {
          for (k in seen) delete seen[k]
          skip_header = 1
          next
        }
        skip_header { skip_header = 0; next }
        {
          key = $1"\t"$2"\t"$3"\t"$4
          if (!(key in seen)) {
            seen[key] = 1
            count[key]++
          }
        }
        END { for (k in count) print k"\t"count[k] }
      ' > "$OUT_DIR/_variant_counts.tsv"
  fi

  # 2c. Dedupe annotated rows by (chrom,pos,ref,alt,transcript).
  # Resumable: skip if a non-empty _cohort.uniq.tab already exists.
  # Resume check: skip ONLY if the file is non-trivially large (a header-only
  # 410-byte stub from an interrupted Stage 2c shouldn't qualify as "done").
  if [[ -s "$OUT_DIR/_cohort.uniq.tab" ]] && [[ $(wc -l < "$OUT_DIR/_cohort.uniq.tab") -gt 1 ]]; then
    log "  2c) skipping (existing _cohort.uniq.tab: $(du -h "$OUT_DIR/_cohort.uniq.tab" | cut -f1))"
  else
    log "  2c) deduplicating annotated rows (gzip-aware)..."
    # Read the header without `head -1` -- that would SIGPIPE zcat and
    # pipefail+errexit would kill the script silently. The subshell with
    # pipefail disabled lets us read the first line and discard the rest.
    header_line=$( set +o pipefail; stream_tab "$first_tab" | head -1 )
    printf '%s\n' "$header_line" > "$OUT_DIR/_cohort.uniq.tab"
    TCOL=$(printf '%s\n' "$header_line" | awk -F'\t' \
      '{ for (i=1;i<=NF;i++) if ($i=="transcript") print i }')
    TCOL=${TCOL:-7}

    # Stream all per-sample files (decompressed if needed) through a cohort-wide
    # dedupe on (chrom,pos,ref,alt,transcript). No per-file reset needed -- this
    # is a cohort-level dedupe, not per-sample. Skip headers via FILE markers.
    (
      for f in "${ann_files[@]}"; do
        printf '@@FILE@@\n'
        case "$f" in
          *.gz) zcat "$f" ;;
          *)    cat  "$f" ;;
        esac
      done
    ) | awk -F'\t' -v tc="$TCOL" '
        /^@@FILE@@$/ { skip_header = 1; next }
        skip_header  { skip_header = 0; next }
        { if (!seen[$1"\t"$2"\t"$3"\t"$4"\t"$tc]++) print }
      ' >> "$OUT_DIR/_cohort.uniq.tab"
  fi

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
  # Resumable: skip if a non-empty cohort.genotypes.tsv already exists.
  if [[ -s "$OUT_DIR/cohort.genotypes.tsv" ]]; then
    log "  2e) skipping (existing cohort.genotypes.tsv: $(du -h "$OUT_DIR/cohort.genotypes.tsv" | cut -f1))"
  else
    log "  2e) concatenating genotype TSVs..."
    shopt -s nullglob
    gt_files=("$OUT_DIR/per_sample/"*.genotypes.tsv)
    shopt -u nullglob
    [[ ${#gt_files[@]} -gt 0 ]] || die "no genotype TSVs found in $OUT_DIR/per_sample/"
    head -1 "${gt_files[0]}" > "$OUT_DIR/cohort.genotypes.tsv"
    awk 'FNR == 1 { next } { print }' "${gt_files[@]}" >> "$OUT_DIR/cohort.genotypes.tsv"
  fi

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

# Keep rule (any one is sufficient):
#   - HIGH or MODERATE consequence (coding hits — always)
#   - SpliceAI ds_max >= 0.2 (cryptic splice signal — always)
#   - clin_sig non-empty (ClinVar opinion — always)
#   - gnomad popmax EXPLICITLY <= threshold (genuinely rare)
#
# Note: missing popmax does NOT bypass the rarity gate. Earlier versions of
# this script treated "missing popmax = keep", which kept >90% of non-coding
# variants and produced a 100+ GB filtered cohort that OOM-killed R. Variants
# without functional signal AND without explicit rarity evidence are dropped.
# Use --filter-allow-missing-popmax to restore the loose behaviour if needed.
awk -F'\t' -v thresh="$FILTER_AF" -v drop_xy="$((1 - INCLUDE_SEX_CHROMS))" \
    -v allow_missing_pmax="$ALLOW_MISSING_POPMAX" '
  NR == 1 {
    for (i=1; i<=NF; i++) col[$i] = i
    print
    next
  }
  {
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
      if (af != "" && af != "." && af+0 <= thresh) keep = 1
      else if (allow_missing_pmax && (af == "" || af == ".")) keep = 1
    }
    if (keep) print
  }
' "$OUT_DIR/cohort.annotated.tab" > "$OUT_DIR/cohort.filtered.tab"

n_pre=$(($(wc -l < "$OUT_DIR/cohort.annotated.tab") - 1))
n_post=$(($(wc -l < "$OUT_DIR/cohort.filtered.tab") - 1))
[[ "$n_pre" -gt 0 ]] && pct=$(awk "BEGIN{printf \"%.1f\", $n_post*100/$n_pre}") || pct="0.0"
log "Pre-filter retained $n_post / $n_pre rows (${pct}%)"

# ============================== Stage 3.5: cohort summary + join =============
# Pre-aggregate the genotype long table into a per-variant summary, then attach
# those columns to the filtered table. This keeps tier_variants.R's input small
# (no need to fread the 30+ GB long genotype table inside R, which OOM-killed
# the script at cohort scale). The R script natively handles inputs that already
# have cohort_ac/cohort_an/cohort_af/n_het/n_hom/n_carriers columns -- it just
# skips its own per_variant_cohort() step when --genotypes is omitted.
section "Stage 3.5: per-variant cohort summary"

# 3.5a: aggregate genotypes
if [[ -s "$OUT_DIR/cohort.variant_summary.tsv" ]]; then
  log "  3.5a) skipping (existing cohort.variant_summary.tsv: $(du -h "$OUT_DIR/cohort.variant_summary.tsv" | cut -f1))"
else
  log "  3.5a) aggregating cohort.genotypes.tsv into per-variant summary..."
  awk -F'\t' '
    NR == 1 { next }
    {
      key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      if      ($6 == "0/0" || $6 == "0|0")                                { called[key]++ }
      else if ($6 == "0/1" || $6 == "1/0" || $6 == "0|1" || $6 == "1|0")  { called[key]++; het[key]++; car[key]++ }
      else if ($6 == "1/1" || $6 == "1|1")                                { called[key]++; hom[key]++; car[key]++ }
    }
    END {
      print "chrom\tpos\tref\talt\tcohort_ac\tcohort_an\tcohort_af\tn_het\tn_hom\tn_carriers"
      for (k in called) {
        split(k, a, SUBSEP)
        h = het[k]+0; H = hom[k]+0
        ac = h + 2*H; an = 2*called[k]; af = (an>0) ? ac/an : 0
        printf "%s\t%s\t%s\t%s\t%d\t%d\t%.6g\t%d\t%d\t%d\n", a[1], a[2], a[3], a[4], ac, an, af, h, H, car[k]+0
      }
    }
  ' "$OUT_DIR/cohort.genotypes.tsv" > "$OUT_DIR/cohort.variant_summary.tsv"
  log "  3.5a) summary rows: $(($(wc -l < "$OUT_DIR/cohort.variant_summary.tsv") - 1))"
fi

# 3.5b: join summary into filtered table
if [[ -s "$OUT_DIR/cohort.filtered.with_cohort.tab" ]]; then
  log "  3.5b) skipping (existing cohort.filtered.with_cohort.tab: $(du -h "$OUT_DIR/cohort.filtered.with_cohort.tab" | cut -f1))"
else
  log "  3.5b) joining cohort summary into filtered table..."
  awk -F'\t' 'BEGIN { OFS="\t" }
    NR == FNR {
      if (FNR == 1) {
        hdr = ""; for (i=5; i<=NF; i++) hdr = hdr "\t" $i
        n_extra = NF - 4
        next
      }
      key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      rec = ""; for (i=5; i<=NF; i++) rec = rec "\t" $i
      summary[key] = rec
      next
    }
    FNR == 1 { print $0 hdr; next }
    {
      key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      if (key in summary) print $0 summary[key]
      else { pad = ""; for (i=1;i<=n_extra;i++) pad = pad "\t"; print $0 pad }
    }
  ' "$OUT_DIR/cohort.variant_summary.tsv" "$OUT_DIR/cohort.filtered.tab" \
    > "$OUT_DIR/cohort.filtered.with_cohort.tab"
fi

# ============================== Stage 4: tier ===============================
section "Stage 4: tiering"
# Note: --genotypes is intentionally omitted here. The cohort columns are already
# attached to cohort.filtered.with_cohort.tab by Stage 3.5; tier_variants.R skips
# its own per_variant_cohort() step when --genotypes is unset, avoiding the OOM
# that happens when R tries to fread the 30+ GB long genotype table.
Rscript "$TIER_SCRIPT" \
  --input         "$OUT_DIR/cohort.filtered.with_cohort.tab" \
  --modifier-genes "$MODIFIER_GENES" \
  --af-max        "$AF_MAX" \
  --out-prefix    "$OUT_DIR/cohort"

# ============================== Stage 5: per-gene burden ====================
# tier_variants.R skips gene_sample_burden when --genotypes is omitted, so we
# compute it externally here via awk: stream the long genotype table, look up
# each carrier's variant in the tiered table, count distinct (gene, sample)
# pairs per qualifying criterion. Join the result back into cohort.genes.tsv.
section "Stage 5: per-gene sample burden"

if [[ -s "$OUT_DIR/cohort.gene_burden.tsv" ]]; then
  log "  5a) skipping (existing cohort.gene_burden.tsv)"
else
  log "  5a) computing per-gene sample burden..."
  awk -F'\t' '
    NR == FNR {
      if (FNR == 1) { for (i=1;i<=NF;i++) c[$i] = i; next }
      key = $c["chrom"] SUBSEP $c["pos"] SUBSEP $c["ref"] SUBSEP $c["alt"]
      g[key]    = $c["gene"]
      tier[key] = $c["tier"]+0
      mod[key]  = $c["modifier_candidate"]
      imp[key]  = $c["impact"]
      vc[key]   = $c["variant_class"]
      next
    }
    FNR == 1 { next }
    $6 ~ /^(0\/1|1\/0|0\|1|1\|0|1\/1|1\|1)$/ {
      key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      if (!(key in g)) next
      gene = g[key]; sample = $5; t = tier[key]
      seen_any[gene, sample] = 1
      if (imp[key] == "HIGH")                seen_hi [gene, sample] = 1
      if (t == 1)                            seen_t1 [gene, sample] = 1
      if (t == 1 || t == 2)                  seen_t12[gene, sample] = 1
      if ((t == 1 || t == 2) || (vc[key] == "pLoF" && t > 0 && t <= 3))
                                             seen_q  [gene, sample] = 1
      if (mod[key] == "TRUE")                seen_m  [gene, sample] = 1
    }
    END {
      print "gene\tn_samples_high_impact\tn_samples_tier1\tn_samples_tier12\tn_samples_qualifying\tn_samples_with_modifier_candidate"
      for (k in seen_any) { split(k, a, SUBSEP); genes[a[1]] = 1 }
      for (gn in genes) {
        hi = t1 = t12 = q = m = 0
        for (k in seen_any) {
          split(k, a, SUBSEP)
          if (a[1] != gn) continue
          if ((gn, a[2]) in seen_hi)  hi++
          if ((gn, a[2]) in seen_t1)  t1++
          if ((gn, a[2]) in seen_t12) t12++
          if ((gn, a[2]) in seen_q)   q++
          if ((gn, a[2]) in seen_m)   m++
        }
        printf "%s\t%d\t%d\t%d\t%d\t%d\n", gn, hi, t1, t12, q, m
      }
    }
  ' "$OUT_DIR/cohort.variants.tsv" "$OUT_DIR/cohort.genotypes.tsv" \
    > "$OUT_DIR/cohort.gene_burden.tsv"
fi

# 5b: join burden columns into cohort.genes.tsv
log "  5b) joining burden into cohort.genes.tsv..."
awk -F'\t' 'BEGIN { OFS="\t" }
  NR == FNR {
    if (FNR == 1) {
      hdr = ""; for (i=2; i<=NF; i++) hdr = hdr "\t" $i
      n_extra = NF - 1
      next
    }
    rec = ""; for (i=2; i<=NF; i++) rec = rec "\t" $i
    burden[$1] = rec
    next
  }
  FNR == 1 { print $0 hdr; next }
  {
    if ($1 in burden) print $0 burden[$1]
    else { pad = ""; for (i=1;i<=n_extra;i++) pad = pad "\t"; print $0 pad }
  }
' "$OUT_DIR/cohort.gene_burden.tsv" "$OUT_DIR/cohort.genes.tsv" \
  > "$OUT_DIR/cohort.genes.with_burden.tsv"

mv "$OUT_DIR/cohort.genes.with_burden.tsv" "$OUT_DIR/cohort.genes.tsv"
log "  5b) merged burden columns into cohort.genes.tsv"

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
