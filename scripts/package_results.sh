#!/usr/bin/env bash
# package_results.sh
# -----------------------------------------------------------------------------
# Bundle cohort tiering outputs into a self-contained tarball for handoff to
# collaborators. Includes the results doc, main tables, figures, validation
# logs, provenance (git commit, manifest, config), and a top-level README
# that tells them where to start.
#
# Explicitly excludes:
#   - cohort.genotypes.tsv        (typically 30-100 GB; kept on server)
#   - cohort.filtered.with_cohort.tab (intermediate; not a deliverable)
#   - per_sample/*                (per-sample fastVEP outputs; intermediate)
#   - figures/qc_ck_*             (streaming checkpoints; not deliverables)
#   - figures/pre-adaptive-afmax/ (superseded figures)
#
# Usage:
#   scripts/package_results.sh \
#       --out-dir     /data/nf1/outputs \
#       --repo        /home/rstudio/src/fastvep-dbnsfp-scores \
#       --manifest    /home/rstudio/src/fastvep-dbnsfp-scores/scripts/nf1_gwas_manifest.v2.tsv \
#       --results-doc docs/nf1_gwas_647_run.md \
#       --cohort-name "NF1 GWAS 647-sample cohort" \
#       --tarball     /data/nf1/deliverables/nf1_gwas_647_results.tar.gz
#
# Only --out-dir and --tarball are required; the rest have sensible defaults.
# -----------------------------------------------------------------------------

set -euo pipefail

# ------------------------------ defaults + args ------------------------------
OUT_DIR=""
REPO=""
MANIFEST=""
RESULTS_DOC=""
COHORT_NAME="cohort"
TARBALL=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir)     OUT_DIR="$2";     shift 2 ;;
    --repo)        REPO="$2";        shift 2 ;;
    --manifest)    MANIFEST="$2";    shift 2 ;;
    --results-doc) RESULTS_DOC="$2"; shift 2 ;;
    --cohort-name) COHORT_NAME="$2"; shift 2 ;;
    --tarball)     TARBALL="$2";     shift 2 ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | sed 's/^# //; s/^#//'
      exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$OUT_DIR" ]] && { echo "--out-dir required" >&2; exit 1; }
[[ -z "$TARBALL" ]] && { echo "--tarball required" >&2; exit 1; }
[[ -d "$OUT_DIR" ]] || { echo "OUT_DIR not found: $OUT_DIR" >&2; exit 1; }

# Default repo: infer from this script's location if not given
if [[ -z "$REPO" ]]; then
  REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

# Default results doc under repo
if [[ -z "$RESULTS_DOC" ]]; then
  # First match under docs/ ending in _run.md; fall back to results_guide.md
  RESULTS_DOC="$(ls "$REPO"/docs/*_run.md 2>/dev/null | head -1)"
  [[ -z "$RESULTS_DOC" ]] && RESULTS_DOC="$REPO/docs/results_guide.md"
fi

# ------------------------------ helpers --------------------------------------
log() { printf '[package] %s\n' "$*" >&2; }
die() { printf '[package] ERROR: %s\n' "$*" >&2; exit 1; }

# ------------------------------ staging area ---------------------------------
STAGE="$(mktemp -d --tmpdir cohort_pkg.XXXXXX)"
trap 'rm -rf "$STAGE"' EXIT

# The tarball's top-level directory name -- derived from the tarball basename
# so extraction produces "<pkgname>/..." rather than dumping into cwd.
PKG_NAME="$(basename "$TARBALL" .tar.gz)"
PKG_NAME="${PKG_NAME%.tgz}"
PKG_DIR="$STAGE/$PKG_NAME"
mkdir -p "$PKG_DIR"/{tables,figures,validation,provenance}

log "staging in $PKG_DIR"

# ------------------------------ tables ---------------------------------------
log "copying tables..."
for t in cohort.variants.tsv cohort.genes.tsv; do
  src="$OUT_DIR/$t"
  [[ -s "$src" ]] || die "missing required table: $src"
  cp "$src" "$PKG_DIR/tables/"
  ( cd "$PKG_DIR/tables" && md5sum "$t" > "${t}.md5" )
done

# ------------------------------ figures --------------------------------------
log "copying figures..."
if [[ -d "$OUT_DIR/figures" ]]; then
  # Cohort tier figures (01_..08_)
  find "$OUT_DIR/figures" -maxdepth 1 -type f \
       \( -name '0[1-8]_*.png' \) -exec cp -t "$PKG_DIR/figures/" {} + 2>/dev/null || true

  # AF source figure + summary
  for f in af_source_distribution.png af_source_summary.tsv; do
    [[ -f "$OUT_DIR/figures/$f" ]] && cp "$OUT_DIR/figures/$f" "$PKG_DIR/figures/"
  done

  # Batch QC figures + tables
  for f in qc_batch_burden_t1.png qc_batch_burden_high_priority.png \
           qc_batch_burden_summary.tsv qc_batch_stats.txt qc_batch_outliers.tsv; do
    [[ -f "$OUT_DIR/figures/$f" ]] && cp "$OUT_DIR/figures/$f" "$PKG_DIR/figures/"
  done
else
  log "WARNING: no figures/ subdirectory under $OUT_DIR -- packaging without figures"
fi

# ------------------------------ validation -----------------------------------
log "copying validation logs..."
for f in validate_cohort.log validate_cohort2.log tier_rerun.log tier_rerun2.log \
         cohort_figures.log af_source.log qc_batch_effects.log; do
  [[ -f "$OUT_DIR/$f" ]] && cp "$OUT_DIR/$f" "$PKG_DIR/validation/"
done

# Re-run the validator on the packaged tables and capture a fresh log. This
# is the definitive record of the tier-table state as shipped -- older logs
# might reflect a previous tier table if the R script was re-run without
# rewriting them.
if [[ -x "$REPO/scripts/validate_cohort.sh" ]]; then
  log "re-running validator against packaged tables..."
  bash "$REPO/scripts/validate_cohort.sh" --out-dir "$OUT_DIR" \
       > "$PKG_DIR/validation/validate_cohort.packaged.log" 2>&1 || true
fi

# ------------------------------ provenance -----------------------------------
log "capturing provenance..."

# Git commit + status of the repo used for the run
{
  echo "# repository"
  echo "path:    $REPO"
  if (cd "$REPO" && git rev-parse HEAD >/dev/null 2>&1); then
    echo "commit:  $(cd "$REPO" && git rev-parse HEAD)"
    echo "branch:  $(cd "$REPO" && git rev-parse --abbrev-ref HEAD)"
    echo "clean:   $(cd "$REPO" && git diff-index --quiet HEAD -- && echo yes || echo NO)"
    echo
    echo "# most recent commits (last 5)"
    (cd "$REPO" && git log --oneline -5)
    if ! (cd "$REPO" && git diff-index --quiet HEAD --); then
      echo
      echo "# UNCOMMITTED CHANGES (repo was not clean at package time)"
      (cd "$REPO" && git status --short)
    fi
  else
    echo "commit:  (repo is not a git working tree)"
  fi
} > "$PKG_DIR/provenance/git_commit.txt"

# Manifest (if provided)
if [[ -n "$MANIFEST" && -f "$MANIFEST" ]]; then
  cp "$MANIFEST" "$PKG_DIR/provenance/manifest.tsv"
fi

# Tiering config -- extract the CFG block from tier_variants.R so the exact
# thresholds are preserved with the deliverable, not just a reference.
if [[ -f "$REPO/R/tier_variants.R" ]]; then
  awk '
    /^CFG <- list\(/ { in_cfg = 1 }
    in_cfg           { print }
    /^\)/ && in_cfg  { in_cfg = 0; exit }
  ' "$REPO/R/tier_variants.R" > "$PKG_DIR/provenance/tier_config.txt"
fi

# Run summary (row counts, tier counts, modifier count) computed on the
# packaged tables -- this is the numbers that MUST match what the doc says.
if [[ -s "$PKG_DIR/tables/cohort.variants.tsv" ]]; then
  awk -F'\t' '
    NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
    { t[$c["tier"]]++; total++
      if ($c["modifier_candidate"] == "TRUE") mod++
    }
    END {
      print "# quick summary of packaged cohort.variants.tsv"
      printf "total variants:    %d\n", total
      for (k=1; k<=5; k++) printf "tier %d:            %d  (%5.2f%%)\n", k, t[k]+0, 100 * (t[k]+0) / total
      printf "modifier_candidate: %d\n", mod+0
    }
  ' "$PKG_DIR/tables/cohort.variants.tsv" > "$PKG_DIR/provenance/summary_counts.txt"
fi

# ------------------------------ results doc + README -------------------------
log "copying results doc + writing README..."
if [[ -f "$RESULTS_DOC" ]]; then
  cp "$RESULTS_DOC" "$PKG_DIR/$(basename "$RESULTS_DOC")"
fi

DOC_BASENAME="$(basename "$RESULTS_DOC")"
DATE="$(date '+%Y-%m-%d')"

cat > "$PKG_DIR/README.md" <<EOF
# ${COHORT_NAME} — results package

Packaged $DATE.

## Where to start

1. \`${DOC_BASENAME}\` — the full results notes. Explains the cohort, the
   tier system, both data-quality fixes applied during this run, known
   limitations, and recommended follow-up. **Read this first.**

2. \`tables/cohort.variants.tsv\` — one row per unique tiered variant. This
   is the primary result. Every row has tier (1-5), modifier_candidate
   flag, af_used, af_source, and full annotation (consequence, gene,
   AlphaMissense, ESM1b, LOEUF, pLI, ClinVar, cohort AC/AN, etc.). See
   the results doc for column-by-column interpretation.

3. \`tables/cohort.genes.tsv\` — per-gene rollup with best_tier and
   sample-burden columns. Use this for gene-level candidate ranking.

4. \`figures/\` — presentation-ready PNGs. Start with:
   - \`01_tier_distribution.png\`
   - \`04_modifier_landscape.png\`  (top genes by modifier candidate count)
   - \`05_nf_anchor_landscape.png\` (NF1/NF2/SMARCB1/LZTR1 positive control)
   - \`qc_batch_burden_high_priority.png\` (batch effect QC — passes)

## Directory contents

\`\`\`
$(cd "$PKG_DIR" && find . -maxdepth 2 -not -path '.' -not -path './.*' | sort)
\`\`\`

## Provenance

- \`provenance/git_commit.txt\` — repo commit SHA + branch used for this run
- \`provenance/tier_config.txt\` — thresholds baked into the tier logic
- \`provenance/manifest.tsv\` — samples included in the run
- \`provenance/summary_counts.txt\` — tier + modifier counts (must match
  what the results doc says)
- \`validation/validate_cohort.packaged.log\` — fresh validator output
  against these exact tables

## Not included

The following are on the server but excluded from this package due to size
or being intermediate outputs:

- \`cohort.genotypes.tsv\` (~30–100 GB long-format genotype table) — request
  separately if the collaborator needs per-sample genotypes at each
  Tier 1+2 variant
- \`cohort.filtered.with_cohort.tab\` (intermediate; not a deliverable)
- \`per_sample/*\` (per-sample fastVEP outputs; intermediate)
EOF

# ------------------------------ tarball --------------------------------------
mkdir -p "$(dirname "$TARBALL")"
log "creating tarball..."
tar -czf "$TARBALL" -C "$STAGE" "$PKG_NAME"

size="$(du -h "$TARBALL" | cut -f1)"
sha="$(sha256sum "$TARBALL" | awk '{print $1}')"

log "wrote $TARBALL ($size)"
log "sha256: $sha"

# Print the tar contents so the caller can verify what shipped without
# needing to untar first.
echo
echo "=== tarball contents ==="
tar -tzf "$TARBALL" | sort
