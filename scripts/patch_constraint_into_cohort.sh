#!/usr/bin/env bash
# patch_constraint_into_cohort.sh
#
# One-off recovery script: joins gnomAD v4.1 gene constraint metrics
# (LOEUF, pLI, mis_z) into an existing cohort.variants.tsv whose loeuf/pli
# columns are empty because csq_to_wide_tab.py (pre-fix) mis-keyed the
# gene-level FV_GNOMAD_GENE projection.
#
# After this fix lands in csq_to_wide_tab.py (commit handling SYMBOL-keyed
# projections), fresh cohort runs do NOT need this script -- it's kept in
# the repo so collaborators with already-built cohort tables can recover
# without re-running Stage 1 across all samples.
#
# Per-transcript constraint rows are deduped by filtering to
# `mane_select == "true"` so each gene contributes one canonical row.
#
# Usage:
#   patch_constraint_into_cohort.sh \
#       <constraint.tsv> <cohort.variants.tsv> [<out.tsv>]
#
# If <out.tsv> is omitted the input is replaced in place (a .bak copy is kept).
#
# Source for the constraint file:
#   https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv

set -euo pipefail

CONSTRAINT="${1:?usage: patch_constraint_into_cohort.sh <constraint.tsv> <cohort.variants.tsv> [<out.tsv>]}"
IN_TAB="${2:?missing cohort.variants.tsv}"
OUT_TAB="${3:-}"

if [[ ! -f "$CONSTRAINT" ]]; then
  echo "constraint file not found: $CONSTRAINT" >&2
  exit 1
fi
if [[ ! -f "$IN_TAB" ]]; then
  echo "cohort variants file not found: $IN_TAB" >&2
  exit 1
fi

IN_PLACE=0
if [[ -z "$OUT_TAB" ]]; then
  IN_PLACE=1
  OUT_TAB="${IN_TAB}.patched.$$"
fi

# Streaming awk two-pass join.  Phase 1 reads the constraint file (small,
# ~20k rows after MANE filter).  Phase 2 streams the cohort table and
# overwrites the loeuf / pli / mis_z columns where the gene matches.
awk -F'\t' 'BEGIN { OFS="\t" }
  NR == FNR {
    if (FNR == 1) {
      for (i = 1; i <= NF; i++) gc[$i] = i
      next
    }
    # Dedupe per-transcript constraint rows to one canonical row per gene.
    if ($gc["mane_select"] != "true") next
    sym = $gc["gene"]
    loeuf = $gc["lof.oe_ci.upper"]
    pli   = $gc["lof.pLI"]
    mis_z = $gc["mis.z_score"]
    constraint[sym] = loeuf "\t" pli "\t" mis_z
    next
  }
  FNR == 1 {
    for (i = 1; i <= NF; i++) c[$i] = i
    if (!("gene"  in c)) { print "ERROR: input lacks `gene` column"  | "cat 1>&2"; exit 1 }
    if (!("loeuf" in c)) { print "ERROR: input lacks `loeuf` column" | "cat 1>&2"; exit 1 }
    if (!("pli"   in c)) { print "ERROR: input lacks `pli` column"   | "cat 1>&2"; exit 1 }
    print
    next
  }
  {
    sym = $c["gene"]
    if (sym in constraint) {
      split(constraint[sym], v, "\t")
      $c["loeuf"] = v[1]
      $c["pli"]   = v[2]
      if ("mis_z" in c) $c["mis_z"] = v[3]
      n_patched++
    } else {
      n_unmatched++
    }
    n_total++
    print
  }
  END {
    printf "[patch] total rows:     %d\n", n_total      | "cat 1>&2"
    printf "[patch] patched:        %d\n", n_patched+0  | "cat 1>&2"
    printf "[patch] unmatched gene: %d\n", n_unmatched+0 | "cat 1>&2"
  }
' "$CONSTRAINT" "$IN_TAB" > "$OUT_TAB"

if [[ "$IN_PLACE" -eq 1 ]]; then
  cp -p "$IN_TAB" "${IN_TAB}.preconstraintpatch.bak"
  mv "$OUT_TAB" "$IN_TAB"
  echo "[patch] in-place; original backed up at ${IN_TAB}.preconstraintpatch.bak"
else
  echo "[patch] wrote $OUT_TAB"
fi
