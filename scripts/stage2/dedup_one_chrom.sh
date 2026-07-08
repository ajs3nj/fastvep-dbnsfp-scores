#!/usr/bin/env bash
# dedup_one_chrom.sh -- per-chrom Phase B of v3 chrom-sharded Stage 2.
#
# Runs the equivalent of v1/v2 Stages 2ab + 2c + 2d against a single
# chromosome's shards from all samples, producing cohort.annotated.chr<N>.tab.
#
# Called in parallel across chromosomes (25 workers oversubscribing 16 cores
# is fine since each is largely I/O-bound).
#
# Usage:
#   dedup_one_chrom.sh <chrom_shard_name> <sharded_input_root> <output_dir>
#
# chrom_shard_name    : "chr1", "chr2", ..., "chrMT"
# sharded_input_root  : per_sample_sharded/  (contains one subdir per sample)
# output_dir          : where cohort.annotated.<chrom>.tab lands
#
# Behavior:
#   - Enumerates $sharded_input_root/<sid>/<chrom_shard_name>.tab.gz across
#     all samples that have that shard
#   - Uses a single awk pass to build BOTH the count hash (2ab) and the
#     transcript-level dedup (2c), producing (chrom_pos_ref_alt, count)
#     and dedupe rows in one go. Two-hash single-pass -- cuts I/O in half
#     compared to v1's split 2ab+2c approach.
#   - Emits cohort.annotated.<chrom>.tab with the n_samples_annotated column
#     appended, matching v1's Stage 2d output schema exactly.

set -euo pipefail

CHROM="${1:?usage: dedup_one_chrom.sh <chrom> <sharded_root> <output_dir>}"
SHARD_ROOT="${2:?usage: dedup_one_chrom.sh <chrom> <sharded_root> <output_dir>}"
OUTDIR="${3:?usage: dedup_one_chrom.sh <chrom> <sharded_root> <output_dir>}"

mkdir -p "$OUTDIR"
OUTFILE="$OUTDIR/cohort.annotated.${CHROM}.tab"

# Idempotency check
if [[ -s "$OUTFILE" && $(wc -l < "$OUTFILE") -gt 1 ]]; then
  echo "[dedup] $CHROM: output already present (${OUTFILE}), skipping" >&2
  exit 0
fi

# Gather all shard files for this chrom, across samples
shopt -s nullglob
shards=("$SHARD_ROOT"/*/"${CHROM}.tab.gz")
shopt -u nullglob

if [[ ${#shards[@]} -eq 0 ]]; then
  echo "[dedup] $CHROM: no shards found under $SHARD_ROOT/*/${CHROM}.tab.gz" >&2
  # Emit an empty file (with just a header) so downstream concat has something
  # to work with. This case happens for e.g. chrY if the cohort is all-female.
  zcat "$SHARD_ROOT"/*/chr1.tab.gz 2>/dev/null | head -1 > "$OUTFILE" || true
  exit 0
fi

n_shards=${#shards[@]}
echo "[dedup] $CHROM: found $n_shards sample shards" >&2

# Determine the transcript column index from the first shard's header.
first_shard="${shards[0]}"
TCOL=$(zcat "$first_shard" | head -1 | awk -F'\t' '{
  for (i=1;i<=NF;i++) if ($i=="transcript") { print i; exit }
}')
: "${TCOL:=7}"

# Single-pass awk that maintains TWO hashes in one stream:
#   count[k4]     : how many samples emit variant k4 = chrom|pos|ref|alt
#   seen[k5]      : uniqueness on (chrom,pos,ref,alt,transcript) --
#                   emit the row the first time we see it, count only
# Then attach n_samples_annotated to the first-seen row.
#
# We also need to remember the FIRST row for each k5 so we can emit it once
# with the correct n_samples_annotated at END. But that's just a per-key
# in-memory copy, which for a single chrom is tractable (millions, not
# hundreds of millions).
#
# Header handling: prepend one header (from the first shard) with
# n_samples_annotated appended. Skip headers in subsequent shards via
# the @@FILE@@ sentinel pattern (same trick v1 uses).

(
  # Header first, with n_samples_annotated appended
  zcat "$first_shard" | head -1 | awk 'BEGIN{OFS="\t"} { print $0, "n_samples_annotated" }'

  # Now stream all shards with sentinel markers between them so per-sample
  # dedup can be applied for the count-per-variant computation.
  for f in "${shards[@]}"; do
    printf '@@FILE@@\n'
    zcat "$f"
  done | awk -F'\t' -v OFS='\t' -v tc="$TCOL" '
    /^@@FILE@@$/ {
      # New sample -- reset per-sample seen hash so we count each (chrom,pos,
      # ref,alt) at most once per sample
      for (k in per_sample_seen) delete per_sample_seen[k]
      skip_header = 1
      next
    }
    skip_header { skip_header = 0; next }
    {
      k4 = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      if (!(k4 in per_sample_seen)) {
        per_sample_seen[k4] = 1
        count[k4]++
      }
      k5 = k4 SUBSEP $tc
      if (!(k5 in seen_row)) {
        seen_row[k5] = $0
      }
    }
    END {
      for (k5 in seen_row) {
        # Recover k4 from k5 (they share the leading 4 SUBSEP-separated fields)
        n = split(k5, parts, SUBSEP)
        k4 = parts[1] SUBSEP parts[2] SUBSEP parts[3] SUBSEP parts[4]
        print seen_row[k5], count[k4]
      }
    }
  '
) > "$OUTFILE"

n_rows=$(( $(wc -l < "$OUTFILE") - 1 ))
echo "[dedup] $CHROM: wrote $n_rows unique variant-transcript rows to $OUTFILE" >&2
