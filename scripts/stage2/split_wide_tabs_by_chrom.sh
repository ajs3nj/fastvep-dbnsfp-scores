#!/usr/bin/env bash
# split_wide_tabs_by_chrom.sh -- per-sample Phase A of v3 chrom-sharded Stage 2.
#
# Reads one per_sample/<sid>.annotated.tab.gz and writes 25 chrom-specific
# gzipped shards under per_sample_sharded/<sid>/. Each shard has the same
# header as the input.
#
# Usage:
#   split_wide_tabs_by_chrom.sh <input.tab.gz> <output_dir> [chrom_col=1]
#
# input.tab.gz : full-cohort-per-sample wide tab (produced by csq_to_wide_tab.py)
# output_dir   : where to write chr*.tab.gz shards; created if missing
# chrom_col    : which TSV column holds the chromosome (default 1)
#
# The output_dir will end up with 25 files:
#   chr1.tab.gz  chr2.tab.gz ... chr22.tab.gz  chrX.tab.gz  chrY.tab.gz  chrMT.tab.gz
#
# Only primary contigs are kept (matches Stage 0 filtering). Rows with any
# other chrom value are silently dropped -- consistent with what fastvep saw
# after Stage 0 normalized the input VCFs.
#
# This script is idempotent: skips if all 25 output files already exist.

set -euo pipefail

INPUT="${1:?usage: split_wide_tabs_by_chrom.sh <input.tab.gz> <output_dir> [chrom_col=1]}"
OUTDIR="${2:?usage: split_wide_tabs_by_chrom.sh <input.tab.gz> <output_dir> [chrom_col=1]}"
CHROM_COL="${3:-1}"

[[ -s "$INPUT" ]] || { echo "ERROR: input not found or empty: $INPUT" >&2; exit 1; }

mkdir -p "$OUTDIR"

# Idempotency check: if all 25 expected outputs exist and are non-empty, skip.
expected=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13
          chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrMT)
all_present=1
for c in "${expected[@]}"; do
  [[ -s "$OUTDIR/${c}.tab.gz" ]] || { all_present=0; break; }
done
if [[ "$all_present" -eq 1 ]]; then
  echo "[split] $INPUT: all 25 shards already present in $OUTDIR, skipping" >&2
  exit 0
fi

# Streaming split via awk. For each row:
#   - line 1: capture header, write to every shard
#   - subsequent: emit to `<outdir>/chr<chrom>.tab.gz` if chrom is primary
#
# We open one pipe per output file. The output files are gzipped inline
# via `| gzip`. That's 25 concurrent gzip processes per sample; on 16-way
# parallel that's ~400 gzip processes, comfortably within default ulimit.
#
# Chromosomes are matched literally: post-Stage 0 the input uses Ensembl-style
# (1, 2, ..., X, Y, MT) with no chr prefix. Output files use `chr<N>` naming
# for clarity of which shard is which.

zcat "$INPUT" | awk -F'\t' -v OUTDIR="$OUTDIR" -v COL="$CHROM_COL" '
  BEGIN {
    OFS = "\t"
    # Map input chrom name -> output shard name
    primary["1"]="chr1"; primary["2"]="chr2"; primary["3"]="chr3"; primary["4"]="chr4"
    primary["5"]="chr5"; primary["6"]="chr6"; primary["7"]="chr7"; primary["8"]="chr8"
    primary["9"]="chr9"; primary["10"]="chr10"; primary["11"]="chr11"; primary["12"]="chr12"
    primary["13"]="chr13"; primary["14"]="chr14"; primary["15"]="chr15"; primary["16"]="chr16"
    primary["17"]="chr17"; primary["18"]="chr18"; primary["19"]="chr19"; primary["20"]="chr20"
    primary["21"]="chr21"; primary["22"]="chr22"; primary["X"]="chrX"
    primary["Y"]="chrY"; primary["MT"]="chrMT"
  }
  NR == 1 {
    # Capture the header. Write it to every primary shard exactly once.
    for (c in primary) {
      cmd = "gzip -c > " OUTDIR "/" primary[c] ".tab.gz"
      print $0 | cmd
    }
    next
  }
  {
    ch = $COL
    if (ch in primary) {
      cmd = "gzip -c > " OUTDIR "/" primary[ch] ".tab.gz"
      print $0 | cmd
    }
    # Silently drop non-primary contigs -- Stage 0 already filtered them
    # from the source VCFs, but be defensive in case any slipped through.
  }
  END {
    # Explicitly close each pipe so gzip processes flush and exit cleanly
    for (c in primary) {
      cmd = "gzip -c > " OUTDIR "/" primary[c] ".tab.gz"
      close(cmd)
    }
  }
'

echo "[split] $INPUT: wrote 25 shards to $OUTDIR" >&2
