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
# Args are passed straight through to cohort_pipeline.sh; this wrapper only
# enforces --skip-annotate so a stray Stage 1 retry can't blow up the cumulative
# per_sample/ output.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Detect whether --skip-annotate was already passed; if not, inject it.
have_skip=0
for arg in "$@"; do
  [[ "$arg" == "--skip-annotate" ]] && { have_skip=1; break; }
done

if [[ "$have_skip" -eq 0 ]]; then
  set -- "$@" --skip-annotate
fi

exec "$SCRIPT_DIR/cohort_pipeline.sh" "$@"
