#!/usr/bin/env bash
# run_batch.sh -- process one batch of samples through Stages 0 and 1.
#
# Called from cohort_pipeline.sh in a loop over batches; can also be run
# standalone to (re-)process a single batch.
#
# Per-sample flow inside a batch:
#   1. Download VCF (if vcf_url is set in manifest; else use vcf_path).
#   2. Stage 0 normalize: bcftools norm -m- + optional chrom rename + bgzip + tabix.
#   3. Stage 1 annotate: fastvep on the normalized VCF -> annotated VCF.
#   4. Stage 1b csq_to_wide_tab: extract the wide tab + gzip it.
#   5. (Optional) LOFTEE pass via ensembl-vep (post-fastvep VCF).
#   6. After all samples in the batch are annotated: delete the raw and
#      normalized VCFs to reclaim disk before the next batch starts.
#
# State sentinels (under $OUT_DIR/state/):
#   batch_<N>.downloaded
#   batch_<N>.normalized
#   batch_<N>.annotated
#   batch_<N>.cleaned
#
# Each stage is skipped if its sentinel already exists; --retry-batch <N> on
# the parent script (or state_clear_batch from lib/batch_state.sh) forces a
# full re-run.
#
# Usage:
#   run_batch.sh \
#       --manifest      manifest.tsv \
#       --batch-id      3 \
#       --out-dir       /data/v2 \
#       --data-dir      /data/v2/ref \
#       --sa-dir        /data/v2/sa \
#       --fasta         /data/v2/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#       --chr-remap     /data/v2/chr_remap.tsv \
#       [--keep-vcfs]                  # skip the cleanup step (debugging)
#       [--threads N]                  # parallel sample-annotation jobs
#       [--loftee]                     # run the post-VEP LOFTEE pass
#       [--loftee-dir DIR]             # LOFTEE plugin directory
#       [--vep-cache DIR]              # Ensembl VEP cache (for LOFTEE)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# ============================== args ==========================================
MANIFEST=""; BATCH_ID=""
OUT_DIR=""; DATA_DIR=""; SA_DIR=""
FASTA=""; GFF3=""; CACHE=""
CHR_REMAP=""
KEEP_VCFS=0
THREADS=$(( $(nproc 2>/dev/null || echo 4) / 2 ))
LOFTEE=0; LOFTEE_DIR=""; VEP_CACHE=""
CSQ_CONVERTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)      MANIFEST="$2";      shift 2 ;;
    --batch-id)      BATCH_ID="$2";      shift 2 ;;
    --out-dir)       OUT_DIR="$2";       shift 2 ;;
    --data-dir)      DATA_DIR="$2";      shift 2 ;;
    --sa-dir)        SA_DIR="$2";        shift 2 ;;
    --fasta)         FASTA="$2";         shift 2 ;;
    --gff3)          GFF3="$2";          shift 2 ;;
    --cache)         CACHE="$2";         shift 2 ;;
    --chr-remap)     CHR_REMAP="$2";     shift 2 ;;
    --keep-vcfs)     KEEP_VCFS=1;        shift ;;
    --threads)       THREADS="$2";       shift 2 ;;
    --loftee)        LOFTEE=1;           shift ;;
    --loftee-dir)    LOFTEE_DIR="$2";    shift 2 ;;
    --vep-cache)     VEP_CACHE="$2";     shift 2 ;;
    --csq-converter) CSQ_CONVERTER="$2"; shift 2 ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
  esac
done

for required in MANIFEST BATCH_ID OUT_DIR FASTA; do
  [[ -z "${!required}" ]] && { echo "ERROR: --${required,,} required" >&2; exit 1; }
done
[[ -z "$CSQ_CONVERTER" ]] && CSQ_CONVERTER="$REPO_DIR/scripts/csq_to_wide_tab.py"

# State tracking
STATE_DIR="$OUT_DIR/state"
source "$SCRIPT_DIR/lib/batch_state.sh"

# Workspace per batch.
# - inputs/, normalized/, annotated/ are per-batch and get deleted on cleanup.
# - per_sample/ is cumulative across batches (same convention as v1
#   cohort_pipeline.sh). The downstream stages glob per_sample/*.annotated.tab(.gz)
#   and per_sample/*.genotypes.tsv -- by writing into that directory we can hand
#   off to run_cohort_stages.sh (or the legacy cohort_pipeline.sh --skip-annotate)
#   without further plumbing.
BATCH_DIR_VCF="$OUT_DIR/inputs/batch_$BATCH_ID"
BATCH_DIR_NORM="$OUT_DIR/normalized/batch_$BATCH_ID"
BATCH_DIR_ANNOT="$OUT_DIR/annotated/batch_$BATCH_ID"
PER_SAMPLE_DIR="$OUT_DIR/per_sample"   # cumulative across batches; not deleted
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$BATCH_DIR_VCF" "$BATCH_DIR_NORM" "$BATCH_DIR_ANNOT" "$PER_SAMPLE_DIR" "$LOG_DIR"

# ============================== helpers ======================================
log()     { printf '[%(%H:%M:%S)T] [batch %s] %s\n' -1 "$BATCH_ID" "$*"; }
section() { echo; log "=== $* ==="; }
die()     { echo "ERROR: $*" >&2; exit 1; }

# Read manifest header to locate columns by name. The manifest is TSV with a
# required header line; columns can appear in any order.
read_header() {
  local header
  header=$(head -1 "$MANIFEST")
  local IFS=$'\t' i col
  i=1
  for col in $header; do
    case "$col" in
      sample_id)     COL_SAMPLE=$i ;;
      batch_id)      COL_BATCH=$i ;;
      vcf_url)       COL_URL=$i ;;
      vcf_path)      COL_PATH=$i ;;
    esac
    i=$((i+1))
  done
  [[ -z "${COL_SAMPLE:-}" ]] && die "manifest missing 'sample_id' column"
}
read_header

# Filter manifest to just this batch's rows. If no batch_id column exists,
# treat all rows as batch 1.
batch_rows() {
  local IFS=
  if [[ -n "${COL_BATCH:-}" ]]; then
    awk -F'\t' -v b="$BATCH_ID" -v ci="$COL_BATCH" \
        'NR>1 && $ci==b { print }' "$MANIFEST"
  elif [[ "$BATCH_ID" == "1" ]]; then
    awk -F'\t' 'NR>1 { print }' "$MANIFEST"
  else
    return 0
  fi
}

n_in_batch=$(batch_rows | wc -l)
[[ "$n_in_batch" -eq 0 ]] && die "no samples found for batch $BATCH_ID in manifest"
log "samples in batch: $n_in_batch"

# ============================== Stage 0: download =============================
if state_done "$BATCH_ID" downloaded; then
  log "skip: download (sentinel exists)"
else
  section "Stage 0a: download"
  if [[ -n "${COL_URL:-}" ]]; then
    while IFS=$'\t' read -r -a fields; do
      sid="${fields[$((COL_SAMPLE-1))]}"
      url="${fields[$((COL_URL-1))]}"
      out="$BATCH_DIR_VCF/$sid.vcf.gz"
      if [[ -s "$out" ]]; then
        log "  $sid: already present, skipping download"
      else
        log "  $sid: downloading from $url"
        # Use the appropriate tool depending on URL scheme. Add more as needed.
        # Synapse (Sage Bionetworks) URLs are supported natively -- entity IDs
        # look like "syn12345678" and the synapse CLI handles authenticated
        # download. Requires the SYNAPSE_AUTH_TOKEN env var (or ~/.synapseConfig).
        case "$url" in
          s3://*)      aws s3 cp "$url" "$out" --no-progress ;;
          gs://*)      gcloud storage cp "$url" "$out" --quiet ;;
          http*|ftp*)  curl -fsSL "$url" -o "$out" ;;
          syn[0-9]*|syn://*)
            # synapse get downloads with the entity's native filename into a
            # target directory; we then rename to the sample_id-based path so
            # downstream stages find it at a predictable location.
            local syn_id="${url#syn://}"
            local tmp_dl
            tmp_dl=$(mktemp -d "$BATCH_DIR_VCF/.dl.$sid.XXXXXX")
            if ! synapse get "$syn_id" --downloadLocation "$tmp_dl" >/dev/null; then
              rm -rf "$tmp_dl"
              die "synapse get failed for $sid ($url) -- check SYNAPSE_AUTH_TOKEN and entity permissions"
            fi
            local downloaded
            downloaded=$(find "$tmp_dl" -maxdepth 1 -type f \( -name '*.vcf.gz' -o -name '*.vcf.bgz' -o -name '*.vcf' \) | head -1)
            [[ -n "$downloaded" ]] || die "synapse get for $sid returned no VCF-shaped file: $(ls "$tmp_dl")"
            mv "$downloaded" "$out"
            # Bring the .tbi along if Synapse stored it too (many Synapse
            # projects store index files as separate entities, in which case
            # tabix is re-created after Stage 0 normalization anyway).
            local idx
            idx=$(find "$tmp_dl" -maxdepth 1 -type f -name '*.tbi' | head -1)
            [[ -n "$idx" ]] && mv "$idx" "$out.tbi"
            rm -rf "$tmp_dl"
            ;;
          *) die "unknown URL scheme for $sid: $url (supported: s3://, gs://, http(s)://, ftp://, syn<id>, syn://<id>)" ;;
        esac
      fi
    done < <(batch_rows)
  elif [[ -n "${COL_PATH:-}" ]]; then
    # Local path mode: symlink instead of copying
    while IFS=$'\t' read -r -a fields; do
      sid="${fields[$((COL_SAMPLE-1))]}"
      src="${fields[$((COL_PATH-1))]}"
      out="$BATCH_DIR_VCF/$sid.vcf.gz"
      [[ -e "$src" ]] || die "vcf not found for $sid: $src"
      [[ -L "$out" || -e "$out" ]] || ln -s "$(readlink -f "$src")" "$out"
    done < <(batch_rows)
  else
    die "manifest has neither 'vcf_url' nor 'vcf_path' column"
  fi
  state_mark "$BATCH_ID" downloaded "samples: $n_in_batch"
fi

# ============================== Stage 0: normalize ============================
if state_done "$BATCH_ID" normalized; then
  log "skip: normalize (sentinel exists)"
else
  section "Stage 0b: normalize (bcftools norm -m-)"
  norm_script="$SCRIPT_DIR/stage0_normalize_vcf.sh"
  [[ -x "$norm_script" ]] || chmod +x "$norm_script"
  export FASTA CHR_REMAP norm_script BATCH_DIR_VCF BATCH_DIR_NORM
  norm_one() {
    local sid="$1"
    local in="$BATCH_DIR_VCF/$sid.vcf.gz"
    local out="$BATCH_DIR_NORM/$sid.norm.vcf.gz"
    [[ -s "$out" ]] && return 0
    local args=(--input "$in" --output "$out" --fasta "$FASTA")
    [[ -n "$CHR_REMAP" ]] && args+=(--chr-remap "$CHR_REMAP")
    "$norm_script" "${args[@]}" >/dev/null
  }
  export -f norm_one
  batch_rows | cut -f"$COL_SAMPLE" \
    | xargs -n 1 -P "$THREADS" -I {} bash -c 'norm_one "$@"' _ {}
  state_mark "$BATCH_ID" normalized "samples: $n_in_batch"
fi

# ============================== Stage 1: annotate =============================
# Per-sample: fastvep annotate -> csq_to_wide_tab.py -> bcftools query genotypes.
# Outputs land in the cumulative $PER_SAMPLE_DIR/, matching the v1 file
# convention so downstream stages can be invoked unchanged.
#   $PER_SAMPLE_DIR/$sid.annotated.tab.gz
#   $PER_SAMPLE_DIR/$sid.genotypes.tsv
if state_done "$BATCH_ID" annotated; then
  log "skip: annotate (sentinel exists)"
else
  section "Stage 1: fastVEP annotate + csq_to_wide_tab + bcftools genotypes"
  export SA_DIR CACHE GFF3 FASTA CSQ_CONVERTER BATCH_DIR_NORM BATCH_DIR_ANNOT \
         PER_SAMPLE_DIR LOG_DIR PYTHON_BIN
  annotate_one() {
    local sid="$1"
    local norm="$BATCH_DIR_NORM/$sid.norm.vcf.gz"
    local ann_vcf="$BATCH_DIR_ANNOT/$sid.vep.vcf"
    local tab="$PER_SAMPLE_DIR/$sid.annotated.tab.gz"
    local gt="$PER_SAMPLE_DIR/$sid.genotypes.tsv"
    local logf="$LOG_DIR/$sid.log"
    # Cached output -- skip if both files are non-empty
    if [[ -s "$tab" && -s "$gt" ]]; then
      echo "[$sid] cached"
      return 0
    fi
    {
      echo "=== $(date) $sid fastvep annotate ==="
      local fastvep_args=(annotate -i "$norm" -o "$ann_vcf" --fasta "$FASTA" --sa-dir "$SA_DIR"
                          --output-format vcf --hgvs --canonical)
      if [[ -n "$CACHE" ]]; then fastvep_args+=(--transcript-cache "$CACHE")
      else                       fastvep_args+=(--gff3 "$GFF3"); fi
      fastvep "${fastvep_args[@]}"

      echo "=== $(date) $sid csq_to_wide_tab ==="
      "${PYTHON_BIN:-python3}" "$CSQ_CONVERTER" \
        --input  "$ann_vcf" \
        --output "$tab" \
        --sample "$sid"

      # The annotated VCF is ~5-10x larger than the wide tab and isn't needed
      # by downstream stages once the tab is written. Drop it now to keep the
      # batch's intermediate disk usage in check.
      rm -f "$ann_vcf"

      echo "=== $(date) $sid bcftools query genotypes ==="
      (
        printf 'chrom\tpos\tref\talt\tsample\tgt\n'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$norm" \
          | awk -v s="$sid" 'BEGIN{OFS="\t"} { print $1, $2, $3, $4, s, $5 }'
      ) > "$gt"

      echo "=== $(date) $sid done ==="
    } > "$logf" 2>&1
    echo "[$sid] done"
  }
  export -f annotate_one
  batch_rows | cut -f"$COL_SAMPLE" \
    | xargs -n 1 -P "$THREADS" -I {} bash -c 'annotate_one "$@"' _ {}
  state_mark "$BATCH_ID" annotated "samples: $n_in_batch"
fi

# ============================== Stage 1c: LOFTEE (optional) ===================
if [[ "$LOFTEE" -eq 1 ]] && ! state_done "$BATCH_ID" loftee; then
  section "Stage 1c: LOFTEE post-VEP pass"
  [[ -z "$LOFTEE_DIR" || -z "$VEP_CACHE" ]] && die "--loftee requires --loftee-dir and --vep-cache"
  export VEP_CACHE LOFTEE_DIR BATCH_DIR_ANNOT BATCH_DIR_NORM PER_SAMPLE_DIR \
         CSQ_CONVERTER LOG_DIR PYTHON_BIN
  loftee_one() {
    local sid="$1"
    local norm="$BATCH_DIR_NORM/$sid.norm.vcf.gz"
    local lof_vcf="$BATCH_DIR_ANNOT/$sid.vep.loftee.vcf"
    local tab="$PER_SAMPLE_DIR/$sid.annotated.tab.gz"
    local logf="$LOG_DIR/$sid.loftee.log"
    {
      echo "=== $(date) $sid LOFTEE pass ==="
      # Run VEP+LOFTEE on the normalized VCF directly (not on the fastvep
      # output, which we already deleted). LOFTEE only needs the LOF column,
      # which we'll join into the wide tab via csq_to_wide_tab.py's existing
      # loftee_lof extraction.
      vep --input_file "$norm" --output_file "$lof_vcf" --vcf \
          --plugin "LoF,loftee_path:$LOFTEE_DIR" \
          --offline --cache --dir_cache "$VEP_CACHE" --no_stats --fork 1

      # Re-extract wide tab. csq_to_wide_tab.py picks up loftee_lof from the
      # CSQ string automatically -- no script changes needed for this step.
      "${PYTHON_BIN:-python3}" "$CSQ_CONVERTER" \
        --input "$lof_vcf" \
        --output "$tab" \
        --sample "$sid"
      rm -f "$lof_vcf"
    } > "$logf" 2>&1
  }
  export -f loftee_one
  batch_rows | cut -f"$COL_SAMPLE" \
    | xargs -n 1 -P "$THREADS" -I {} bash -c 'loftee_one "$@"' _ {}
  state_mark "$BATCH_ID" loftee "samples: $n_in_batch"
fi

# ============================== Stage 0: cleanup ==============================
# Delete the raw and normalized VCFs to reclaim disk. The wide tabs in
# $BATCH_DIR_TAB are cumulative and stay.
if [[ "$KEEP_VCFS" -eq 1 ]]; then
  log "skip: cleanup (--keep-vcfs)"
else
  if state_done "$BATCH_ID" cleaned; then
    log "skip: cleanup (sentinel exists)"
  else
    section "Stage 0z: cleanup (reclaim disk)"
    # SAFETY: only proceed if the annotated sentinel exists -- never delete
    # source VCFs unless we've confirmed every sample in the batch has a wide tab.
    state_done "$BATCH_ID" annotated || die "refusing cleanup: annotate sentinel missing"
    # Verify every sample has a non-empty wide tab AND genotype TSV before
    # deleting source VCFs. If either is missing, the run is incomplete and
    # we must not lose the source data.
    while IFS=$'\t' read -r -a fields; do
      sid="${fields[$((COL_SAMPLE-1))]}"
      [[ -s "$PER_SAMPLE_DIR/$sid.annotated.tab.gz" ]] || die "refusing cleanup: $sid.annotated.tab.gz missing or empty"
      [[ -s "$PER_SAMPLE_DIR/$sid.genotypes.tsv"    ]] || die "refusing cleanup: $sid.genotypes.tsv missing or empty"
    done < <(batch_rows)
    # Disk before/after
    before=$(du -sh "$BATCH_DIR_VCF" "$BATCH_DIR_NORM" "$BATCH_DIR_ANNOT" 2>/dev/null | awk '{ sum+=$1 } END { print sum }')
    rm -rf "$BATCH_DIR_VCF" "$BATCH_DIR_NORM" "$BATCH_DIR_ANNOT"
    log "  reclaimed disk: $before"
    state_mark "$BATCH_ID" cleaned "reclaimed: $before"
  fi
fi

log "batch $BATCH_ID complete"
