#!/usr/bin/env bash
# cohort_pipeline_batched.sh -- v2 top-level batched orchestrator.
#
# Reads a v2 manifest (sample_id, batch_id, vcf_url, [ancestry, severity, ...])
# and runs the pipeline in batches:
#
#   for each batch:
#       download VCFs -> normalize -> annotate (+ optional LOFTEE) -> cleanup
#   once all batches done:
#       Stage 2 cohort merge
#       Stage 3 pre-filter
#       Stage 3.5 cohort summary + genotypes
#       Stage 4 R tier
#       Stage 5 per-gene burden
#
# Disk discipline: each batch reclaims its raw + normalized VCFs after
# annotation, so peak disk stays bounded by the largest single batch (~1-2 TB
# for a 50-sample batch) rather than the whole cohort.
#
# State sentinels under $OUT_DIR/state/ make the run resumable; killing it
# mid-batch and restarting picks up where it left off.
#
# Usage:
#   cohort_pipeline_batched.sh \
#       --manifest    manifest.v2.tsv \
#       --out-dir     /data/v2 \
#       --data-dir    /data/v2/ref \
#       --sa-dir      /data/v2/sa \
#       --fasta       /data/v2/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#       [--chr-remap  /data/v2/chr_remap.tsv]
#       [--threads    16]
#       [--loftee] [--loftee-dir DIR] [--vep-cache DIR]
#       [--keep-vcfs]                # disable disk-reclamation step (debug)
#       [--retry-batch N]            # force re-run of batch N (clears sentinels)
#       [--only-batch N]             # process JUST batch N then stop
#       [--skip-cohort-stages]       # batches only; skip Stages 2-5 (debug)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# ============================== args =========================================
MANIFEST=""
OUT_DIR=""; DATA_DIR=""; SA_DIR=""
FASTA=""; GFF3=""; CACHE=""; CHR_REMAP=""
THREADS=$(( $(nproc 2>/dev/null || echo 4) / 2 ))
LOFTEE=0; LOFTEE_DIR=""; VEP_CACHE=""
KEEP_VCFS=0
RETRY_BATCH=""
ONLY_BATCH=""
SKIP_COHORT_STAGES=0
PASSTHROUGH=()     # unrecognized args forwarded to run_cohort_stages.sh

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)            MANIFEST="$2";        shift 2 ;;
    --out-dir)             OUT_DIR="$2";         shift 2 ;;
    --data-dir)            DATA_DIR="$2";        shift 2 ;;
    --sa-dir)              SA_DIR="$2";          shift 2 ;;
    --fasta)               FASTA="$2";           shift 2 ;;
    --gff3)                GFF3="$2";            shift 2 ;;
    --cache)               CACHE="$2";           shift 2 ;;
    --chr-remap)           CHR_REMAP="$2";       shift 2 ;;
    --threads)             THREADS="$2";         shift 2 ;;
    --loftee)              LOFTEE=1;             shift ;;
    --loftee-dir)          LOFTEE_DIR="$2";      shift 2 ;;
    --vep-cache)           VEP_CACHE="$2";       shift 2 ;;
    --keep-vcfs)           KEEP_VCFS=1;          shift ;;
    --retry-batch)         RETRY_BATCH="$2";     shift 2 ;;
    --only-batch)          ONLY_BATCH="$2";      shift 2 ;;
    --skip-cohort-stages)  SKIP_COHORT_STAGES=1; shift ;;
    # Forward everything else (e.g. --af-max, --filter-af, --include-sex-chroms,
    # --filter-allow-missing-popmax, --modifier-genes) to run_cohort_stages.sh
    *)                     PASSTHROUGH+=("$1");  shift ;;
  esac
done

for required in MANIFEST OUT_DIR FASTA; do
  [[ -z "${!required}" ]] && { echo "ERROR: --${required,,} required" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
STATE_DIR="$OUT_DIR/state"
source "$SCRIPT_DIR/lib/batch_state.sh"

# ============================== helpers ======================================
log()     { printf '[%(%H:%M:%S)T] %s\n' -1 "$*"; }
section() { echo; log "=== $* ==="; }

# Extract distinct batch IDs from the manifest. If no batch_id column exists,
# every sample is treated as batch 1.
list_batches() {
  awk -F'\t' '
    NR == 1 {
      for (i=1;i<=NF;i++) col[$i] = i
      if (!("batch_id" in col)) { print "1"; exit }
      next
    }
    { print $col["batch_id"] }
  ' "$MANIFEST" | sort -u -V
}

# ============================== validation ==================================
section "Validation"
[[ -e "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }
[[ -e "$FASTA"    ]] || { echo "ERROR: fasta not found: $FASTA" >&2; exit 1; }

# Sanity check that the manifest has the required column
sample_col=$(awk -F'\t' 'NR==1 { for(i=1;i<=NF;i++) if ($i=="sample_id") print i }' "$MANIFEST")
[[ -n "$sample_col" ]] || { echo "ERROR: manifest missing 'sample_id' column" >&2; exit 1; }

batches=($(list_batches))
n_batches=${#batches[@]}
log "batches detected: $n_batches  (${batches[*]})"

# ============================== retry/only handling ==========================
if [[ -n "$RETRY_BATCH" ]]; then
  log "retry-batch: clearing all sentinels for batch $RETRY_BATCH"
  state_clear_batch "$RETRY_BATCH"
fi
if [[ -n "$ONLY_BATCH" ]]; then
  log "only-batch: limiting to batch $ONLY_BATCH"
  batches=("$ONLY_BATCH")
fi

# ============================== per-batch loop ===============================
section "Per-batch processing"
for batch in "${batches[@]}"; do
  log ">>> batch $batch starting"
  batch_args=(--manifest "$MANIFEST" --batch-id "$batch"
              --out-dir "$OUT_DIR" --data-dir "$DATA_DIR" --sa-dir "$SA_DIR"
              --fasta "$FASTA" --threads "$THREADS")
  [[ -n "$GFF3"       ]] && batch_args+=(--gff3 "$GFF3")
  [[ -n "$CACHE"      ]] && batch_args+=(--cache "$CACHE")
  [[ -n "$CHR_REMAP"  ]] && batch_args+=(--chr-remap "$CHR_REMAP")
  [[ "$KEEP_VCFS" -eq 1 ]] && batch_args+=(--keep-vcfs)
  if [[ "$LOFTEE" -eq 1 ]]; then
    batch_args+=(--loftee --loftee-dir "$LOFTEE_DIR" --vep-cache "$VEP_CACHE")
  fi
  "$SCRIPT_DIR/run_batch.sh" "${batch_args[@]}"
  log "<<< batch $batch done"
done

# ============================== cohort stages 2-5 ============================
if [[ "$SKIP_COHORT_STAGES" -eq 1 ]]; then
  log "skipping cohort stages 2-5 (--skip-cohort-stages)"
  exit 0
fi
if [[ -n "$ONLY_BATCH" ]]; then
  log "only-batch: skipping cohort stages 2-5"
  exit 0
fi

section "Cohort stages 2-5"
cohort_args=(--manifest "$MANIFEST" --out-dir "$OUT_DIR"
             --data-dir "$DATA_DIR" --sa-dir "$SA_DIR"
             --fasta "$FASTA" --threads "$THREADS")
[[ -n "$GFF3"  ]] && cohort_args+=(--gff3 "$GFF3")
[[ -n "$CACHE" ]] && cohort_args+=(--cache "$CACHE")
cohort_args+=("${PASSTHROUGH[@]}")
"$SCRIPT_DIR/run_cohort_stages.sh" "${cohort_args[@]}"

section "Pipeline complete"
log "state summary:"
state_summary
log "outputs: $OUT_DIR/cohort.variants.tsv, $OUT_DIR/cohort.genes.tsv, $OUT_DIR/figures/"
