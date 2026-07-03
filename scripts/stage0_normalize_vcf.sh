#!/usr/bin/env bash
# stage0_normalize_vcf.sh -- per-sample VCF normalization.
#
# Runs as the first step in a batch loop, before fastVEP annotation.
# Produces a bi-allelic-split, left-anchored, Ensembl-chrom-named, indexed VCF
# suitable for downstream pipeline stages.
#
# This script addresses three real bugs we hit in v1:
#
#   1. Multi-allelic VCF rows ("alt=C,G") don't match the bi-allelic keys our
#      tiering output uses. bcftools norm -m- splits them into one row per alt.
#
#   2. Indel representations differ between encoding conventions. We left-anchor
#      via bcftools norm so downstream variant_keys match cohort.genotypes.tsv
#      and cohort.variants.tsv consistently.
#
#   3. UCSC-style chrom names ("chr1") vs Ensembl-style ("1") -- if the input
#      uses chr-prefixed names but downstream tools expect Ensembl, every
#      annotation looks intergenic. Apply a rename map up front.
#
# Usage:
#   stage0_normalize_vcf.sh \
#       --input  S001.vcf.gz \
#       --output S001.normalized.vcf.gz \
#       --fasta  GRCh38.fa \
#       [--chr-remap chr_remap.tsv]
#
# The chr-remap file is optional; if omitted and the input has chr-prefixed
# chroms, the script logs a warning. Defaults are Ensembl-style; if you're
# using UCSC-style downstream, your tools should handle it consistently.

set -euo pipefail

INPUT=""; OUTPUT=""; FASTA=""; CHR_REMAP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)     INPUT="$2";     shift 2 ;;
    --output)    OUTPUT="$2";    shift 2 ;;
    --fasta)     FASTA="$2";     shift 2 ;;
    --chr-remap) CHR_REMAP="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | sed 's/^# //; s/^#//'
      exit 0
      ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$INPUT"  ]] && { echo "ERROR: --input required"  >&2; exit 1; }
[[ -z "$OUTPUT" ]] && { echo "ERROR: --output required" >&2; exit 1; }
[[ -z "$FASTA"  ]] && { echo "ERROR: --fasta required (for left-anchoring)" >&2; exit 1; }
[[ -e "$INPUT"  ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[[ -e "$FASTA"  ]] || { echo "ERROR: fasta not found: $FASTA" >&2; exit 1; }
[[ -e "$FASTA.fai" ]] || { echo "ERROR: fasta index missing: $FASTA.fai (run samtools faidx)" >&2; exit 1; }

# Detect chrom-naming convention in input. Used both for the warning if no
# chr-remap was supplied AND for a sanity-check against the fasta convention.
first_chrom=$(bcftools view -h "$INPUT" 2>/dev/null \
  | awk '/^##contig=<ID=/ { sub(/^##contig=<ID=/, ""); sub(/,.*/, ""); print; exit }')
fasta_first=$(awk '/^>/ { sub(/^>/, ""); sub(/[[:space:]].*/, ""); print; exit }' "$FASTA")
case "$first_chrom" in
  chr*) input_style="ucsc" ;;
  *)    input_style="ensembl" ;;
esac
case "$fasta_first" in
  chr*) fasta_style="ucsc" ;;
  *)    fasta_style="ensembl" ;;
esac

if [[ "$input_style" != "$fasta_style" ]]; then
  if [[ -z "$CHR_REMAP" ]]; then
    echo "ERROR: chrom-naming mismatch -- input is $input_style ('$first_chrom'), fasta is $fasta_style ('$fasta_first'). Pass --chr-remap to fix." >&2
    exit 1
  fi
fi

mkdir -p "$(dirname "$OUTPUT")"
tmp=$(mktemp -d)
# On failure, preserve the tmp dir + its norm.log so the caller can see the
# actual bcftools error message. On success, clean up as usual.
cleanup() {
  local status=$?
  if [[ $status -ne 0 ]]; then
    echo "[stage0] FAILED (exit $status) -- preserving bcftools log at $tmp/norm.log" >&2
    if [[ -s "$tmp/norm.log" ]]; then
      echo "----- $tmp/norm.log (last 20 lines) -----" >&2
      tail -20 "$tmp/norm.log" >&2
      echo "----- end -----" >&2
    fi
  else
    rm -rf "$tmp"
  fi
}
trap cleanup EXIT

# Pipeline:
#   1. bcftools norm -m- -f FASTA : split multi-allelic, left-anchor indels
#   2. bcftools annotate --rename-chrs (optional): rename chroms if requested
#   3. bgzip + tabix index
log() { printf '[stage0] %s\n' "$*" >&2; }

log "input: $INPUT  ($input_style chroms)"
log "fasta: $FASTA  ($fasta_style chroms)"
[[ -n "$CHR_REMAP" ]] && log "chr remap: $CHR_REMAP"

if [[ -n "$CHR_REMAP" ]]; then
  # Rename FIRST, then normalize. bcftools norm -f FASTA needs the contig
  # names in the VCF to match the FASTA index -- if we run norm before
  # renaming, it looks up "chr1" against an Ensembl-style FASTA and fails
  # with "faidx_fetch_seq failed at chr1:NNN". Correct order is
  # annotate --rename-chrs -> norm.
  bcftools annotate --rename-chrs "$CHR_REMAP" -Ou "$INPUT" 2> "$tmp/norm.log" \
    | bcftools norm -m- -f "$FASTA" -Oz -o "$OUTPUT" - 2>> "$tmp/norm.log"
else
  bcftools norm -m- -f "$FASTA" -Oz -o "$OUTPUT" "$INPUT" 2> "$tmp/norm.log"
fi

tabix -p vcf -f "$OUTPUT"

# Quick sanity stats from the bcftools norm log
log "$(grep -E 'Lines|Total' "$tmp/norm.log" | tr '\n' ' | ')"
log "wrote: $OUTPUT (indexed)"
