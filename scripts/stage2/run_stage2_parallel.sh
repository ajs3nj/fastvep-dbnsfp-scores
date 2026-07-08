#!/usr/bin/env bash
# run_stage2_parallel.sh -- v3 chrom-sharded parallel Stage 2 orchestrator.
#
# Replaces v1/v2's single-threaded cohort_pipeline.sh Stage 2 (which spent
# ~30 hours streaming through 647 samples serially). Three-phase:
#
#   A. Split each per_sample/<sid>.annotated.tab.gz by chromosome, parallel
#      across samples (16-way default). Writes per_sample_sharded/<sid>/chr*.tab.gz.
#   B. Per-chrom dedup + count-attach, parallel across chromosomes (25-way).
#      Writes stage2_shards/cohort.annotated.chr*.tab. Each chrom's job reads
#      only that chrom's shards from all 647 samples (~1/25th of the I/O of
#      the whole-cohort awk).
#   C. Concatenate the 25 per-chrom outputs into cohort.annotated.tab.
#
# Compatible with v1's cohort_pipeline.sh Stage 2 output schema, so
# downstream stages (3, 3.5, 4, 5) work unchanged. Genotype concat (Stage 2e
# in v1) is done unchanged -- it's already fast enough at cohort scale.
#
# Usage:
#   run_stage2_parallel.sh \
#     --out-dir     /data/nf1/outputs \
#     [--threads    16]                 # samples in Phase A + chroms in Phase B
#     [--phase-a-only]                  # stop after splitting per_sample
#     [--phase-b-only]                  # skip Phase A (assume shards exist)
#     [--phase-c-only]                  # only concat existing chrom outputs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUT_DIR=""
THREADS=$(nproc 2>/dev/null || echo 8)
DO_A=1
DO_B=1
DO_C=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir)      OUT_DIR="$2"; shift 2 ;;
    --threads)      THREADS="$2"; shift 2 ;;
    --phase-a-only) DO_B=0; DO_C=0; shift ;;
    --phase-b-only) DO_A=0;         shift ;;
    --phase-c-only) DO_A=0; DO_B=0; shift ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$OUT_DIR" ]] && { echo "ERROR: --out-dir required" >&2; exit 1; }
[[ -d "$OUT_DIR/per_sample" ]] || { echo "ERROR: $OUT_DIR/per_sample missing" >&2; exit 1; }

PER_SAMPLE_DIR="$OUT_DIR/per_sample"
SHARDED_DIR="$OUT_DIR/per_sample_sharded"
CHROM_STAGE_DIR="$OUT_DIR/stage2_shards"
mkdir -p "$SHARDED_DIR" "$CHROM_STAGE_DIR"

log()     { printf '[%(%H:%M:%S)T] %s\n' -1 "$*"; }
section() { echo; log "=== $* ==="; }

# =========================================================================
# Phase A -- split per-sample wide tabs by chromosome
# =========================================================================
if [[ "$DO_A" -eq 1 ]]; then
  section "Phase A: splitting per-sample wide tabs by chromosome (${THREADS}-way parallel)"

  # Build the sample list from what's actually on disk
  mapfile -t samples < <(
    find "$PER_SAMPLE_DIR" -maxdepth 1 -name '*.annotated.tab.gz' -printf '%f\n' \
      | sed 's/\.annotated\.tab\.gz$//'
  )
  n_samples=${#samples[@]}
  log "found $n_samples per-sample wide tabs to shard"

  export SHARDED_DIR PER_SAMPLE_DIR
  export SPLITTER="$SCRIPT_DIR/split_wide_tabs_by_chrom.sh"
  chmod +x "$SPLITTER" 2>/dev/null || true

  split_one() {
    local sid="$1"
    local input="$PER_SAMPLE_DIR/${sid}.annotated.tab.gz"
    local outdir="$SHARDED_DIR/$sid"
    mkdir -p "$outdir"
    "$SPLITTER" "$input" "$outdir" 1 > /dev/null 2>&1
    echo "[splitA] done $sid"
  }
  export -f split_one

  printf '%s\n' "${samples[@]}" \
    | xargs -n 1 -P "$THREADS" -I {} bash -c 'split_one "$@"' _ {}

  log "Phase A complete: $(ls "$SHARDED_DIR" | wc -l) sample subdirs written"
fi

# =========================================================================
# Phase B -- per-chrom dedup
# =========================================================================
if [[ "$DO_B" -eq 1 ]]; then
  section "Phase B: per-chrom cohort dedup (25-way parallel across chromosomes)"

  chroms=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
          chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
          chr21 chr22 chrX chrY chrMT)

  export SHARDED_DIR CHROM_STAGE_DIR
  export DEDUPER="$SCRIPT_DIR/dedup_one_chrom.sh"
  chmod +x "$DEDUPER" 2>/dev/null || true

  dedup_one() {
    local chrom="$1"
    "$DEDUPER" "$chrom" "$SHARDED_DIR" "$CHROM_STAGE_DIR" > /dev/null 2>&1
    echo "[dedupB] done $chrom"
  }
  export -f dedup_one

  # Oversubscribe: 25 workers on ~16 cores, since dedup is I/O-bound.
  # Adjust downward if the machine has less headroom.
  printf '%s\n' "${chroms[@]}" \
    | xargs -n 1 -P 25 -I {} bash -c 'dedup_one "$@"' _ {}

  log "Phase B complete: $(ls "$CHROM_STAGE_DIR" | wc -l) per-chrom outputs written"
fi

# =========================================================================
# Phase C -- concat into cohort.annotated.tab
# =========================================================================
if [[ "$DO_C" -eq 1 ]]; then
  section "Phase C: concat per-chrom outputs into cohort.annotated.tab"

  # Deterministic chromosome order for the concat: 1..22, X, Y, MT
  chroms_ordered=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
                  chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
                  chr21 chr22 chrX chrY chrMT)

  # Header from first non-empty chrom output
  first_chrom=""
  for c in "${chroms_ordered[@]}"; do
    if [[ -s "$CHROM_STAGE_DIR/cohort.annotated.${c}.tab" ]]; then
      first_chrom="$c"; break
    fi
  done
  [[ -z "$first_chrom" ]] && { echo "ERROR: no chrom outputs found in $CHROM_STAGE_DIR" >&2; exit 1; }

  head -1 "$CHROM_STAGE_DIR/cohort.annotated.${first_chrom}.tab" > "$OUT_DIR/cohort.annotated.tab"

  for c in "${chroms_ordered[@]}"; do
    if [[ -s "$CHROM_STAGE_DIR/cohort.annotated.${c}.tab" ]]; then
      tail -n+2 "$CHROM_STAGE_DIR/cohort.annotated.${c}.tab"
    fi
  done >> "$OUT_DIR/cohort.annotated.tab"

  n_rows=$(( $(wc -l < "$OUT_DIR/cohort.annotated.tab") - 1 ))
  log "Phase C complete: cohort.annotated.tab has $n_rows unique variant-transcript rows"
fi

section "Stage 2 parallel done"
log "Output: $OUT_DIR/cohort.annotated.tab"
