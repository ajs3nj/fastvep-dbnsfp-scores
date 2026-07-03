#!/usr/bin/env bash
# scripts/lib/batch_state.sh -- sentinel-based state tracking for batched runs.
#
# Source this from any script that needs to read or write batch state:
#   source "$(dirname "${BASH_SOURCE[0]}")/lib/batch_state.sh"
#
# Sentinel files live under $STATE_DIR (typically $OUT_DIR/state/) and are
# written atomically (write-to-temp + rename). The convention is one sentinel
# per (batch, stage) pair:
#
#   state/batch_<N>.downloaded
#   state/batch_<N>.normalized
#   state/batch_<N>.annotated
#   state/batch_<N>.cleaned
#
# Each sentinel's content is a one-line summary (timestamp, sample count, etc.)
# for human inspection. The rename-into-place pattern means a half-written
# sentinel never confuses a resumed run.

set -u  # don't set -e here -- sourced into different shells

# Caller must set STATE_DIR before sourcing.
: "${STATE_DIR:?STATE_DIR must be set before sourcing batch_state.sh}"

mkdir -p "$STATE_DIR"

# state_path BATCH STAGE -> prints the sentinel path
state_path() {
  local batch="$1" stage="$2"
  printf '%s/batch_%s.%s\n' "$STATE_DIR" "$batch" "$stage"
}

# state_mark BATCH STAGE [extra...] -- atomically write the sentinel
state_mark() {
  local batch="$1" stage="$2"; shift 2
  local target tmp
  target=$(state_path "$batch" "$stage")
  tmp="$target.tmp.$$"
  {
    printf 'timestamp:    %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'batch:        %s\n' "$batch"
    printf 'stage:        %s\n' "$stage"
    [[ $# -gt 0 ]] && printf '%s\n' "$@"
  } > "$tmp"
  mv "$tmp" "$target"
}

# state_done BATCH STAGE -> returns 0 if the sentinel exists, 1 otherwise
state_done() {
  local batch="$1" stage="$2"
  [[ -e "$(state_path "$batch" "$stage")" ]]
}

# state_clear BATCH STAGE -- remove a sentinel (used by --retry-batch)
state_clear() {
  local batch="$1" stage="$2"
  rm -f "$(state_path "$batch" "$stage")"
}

# state_clear_batch BATCH -- remove ALL sentinels for a batch (full retry)
state_clear_batch() {
  local batch="$1"
  rm -f "$STATE_DIR"/batch_"$batch".*
}

# state_list_batches -- prints all distinct batch IDs that have any sentinel
state_list_batches() {
  ls "$STATE_DIR" 2>/dev/null \
    | sed -n 's/^batch_\([^.]*\)\..*/\1/p' \
    | sort -u
}

# state_summary -- pretty-print the status of each batch (which stages done)
state_summary() {
  local b stages stages_str
  for b in $(state_list_batches); do
    stages=$(ls "$STATE_DIR" 2>/dev/null \
      | sed -n "s/^batch_${b}\.//p" \
      | tr '\n' ',' | sed 's/,$//')
    printf '  batch %-8s  %s\n' "$b" "$stages"
  done
}
