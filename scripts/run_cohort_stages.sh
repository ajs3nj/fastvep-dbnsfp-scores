#!/usr/bin/env bash
# run_cohort_stages.sh -- runs Stages 2-5 on an already-populated $OUT_DIR.
#
# Expects: $OUT_DIR/per_sample/*.annotated.tab.gz  (cumulative across batches)
#          $OUT_DIR/per_sample/*.genotypes.tsv     (cumulative across batches)
#
# Thin wrapper around the v1 cohort_pipeline.sh invoked with --skip-annotate.
# Keeps the battle-tested Stage 2-5 logic without duplicating it; v2's batched
# orchestrator calls this once after all batches have finished annotation.
#
# Two responsibilities:
#   1. Force --skip-annotate so a stray Stage 1 retry can't blow up the
#      cumulative per_sample/ output.
#   2. Bridge the v2 manifest format (sample_id<TAB>batch_id<TAB>vcf_url<...>)
#      to v1's expected format (sample_id<TAB>vcf_path). We build a temporary
#      v1-compatible manifest with dummy vcf_paths -- with --skip-annotate,
#      cohort_pipeline.sh only needs the sample_id column, but its upfront
#      validation still hard-checks the two-column header.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Extract --manifest and --out-dir values from args. We need both for the
# format bridge -- --manifest to translate, --out-dir to build the dummy
# vcf_path values that point into the actual per_sample/ dir.
manifest_in=""
out_dir=""
have_skip=0
have_manifest_flag=0
new_args=()
i=0
args=("$@")
while [[ $i -lt ${#args[@]} ]]; do
  case "${args[$i]}" in
    --manifest)
      manifest_in="${args[$((i+1))]}"
      have_manifest_flag=1
      i=$((i+2))
      ;;
    --out-dir)
      out_dir="${args[$((i+1))]}"
      new_args+=("${args[$i]}" "${args[$((i+1))]}")
      i=$((i+2))
      ;;
    --skip-annotate)
      have_skip=1
      new_args+=("${args[$i]}")
      i=$((i+1))
      ;;
    *)
      new_args+=("${args[$i]}")
      i=$((i+1))
      ;;
  esac
done

# If the caller supplied a manifest, translate it. Otherwise skip the bridge.
tmp_manifest=""
if [[ -n "$manifest_in" ]]; then
  [[ -e "$manifest_in" ]] || { echo "ERROR: manifest not found: $manifest_in" >&2; exit 1; }

  # Detect whether this is already a v1 manifest (has vcf_path in header).
  # If so, pass through unchanged; if it's v2, translate.
  header=$(head -1 "$manifest_in")
  if grep -qE $'\t''vcf_path(\t|$)' <(printf '%s\n' "$header"); then
    # Already v1 format -- pass through
    new_args+=(--manifest "$manifest_in")
  else
    # v2 format -- translate to a temp v1 manifest with dummy vcf_path
    tmp_manifest=$(mktemp --tmpdir cohort_manifest_v1.XXXXXX.tsv)
    trap 'rm -f "$tmp_manifest"' EXIT
    awk -F'\t' -v OFS='\t' -v out_dir="$out_dir" '
      NR == 1 {
        for (i=1;i<=NF;i++) c[$i] = i
        if (!("sample_id" in c)) {
          print "ERROR: manifest missing sample_id column" > "/dev/stderr"
          exit 1
        }
        print "sample_id","vcf_path"
        next
      }
      # Emit sample_id + a dummy path. With --skip-annotate the path is
      # never opened for annotation, but cohort_pipeline.sh does an upfront
      # chromosome-naming check that reads $2 of the first data row and
      # zcat-checks the first non-header chrom in it. Point at the already-
      # produced per_sample wide tab -- gzipped, first data row starts with
      # an Ensembl-style chrom (1, 2, ...), never gets opened for anything
      # else because Stage 1 is skipped.
      { print $c["sample_id"], out_dir "/per_sample/" $c["sample_id"] ".annotated.tab.gz" }
    ' "$manifest_in" > "$tmp_manifest"
    new_args+=(--manifest "$tmp_manifest")
    echo "[run_cohort_stages] translated v2 manifest -> $tmp_manifest" >&2
  fi
fi

# Force --skip-annotate if not already present
if [[ "$have_skip" -eq 0 ]]; then
  new_args+=(--skip-annotate)
fi

exec "$SCRIPT_DIR/cohort_pipeline.sh" "${new_args[@]}"
