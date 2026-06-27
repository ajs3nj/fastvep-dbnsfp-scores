# fastvep-dbnsfp-scores

A variant annotation and tiering pipeline for rare-disease cohorts, built and validated on a
209-sample neurofibromatosis type 1 (NF1) WGS study. Annotates variants with calibrated
missense pathogenicity predictors (AlphaMissense, ESM1b, REVEL), assigns class-conditional
priority tiers, flags modifier candidates as a parallel track, and produces collaborator-ready
per-variant and per-gene cohort summaries with clinical-presentation figures.

If you received this repo as a handoff (tarball or git clone), this README is the place to
start for a quick orientation. For the **column-by-column dictionary of every output file
and the complete tiering rules table**, see `docs/results_guide.md`. For comprehensive
methods, see `docs/pipeline_methods.md`.

---

## What the outputs mean

The cohort pipeline produces two main tables and a set of figures.

### `cohort.variants.tsv` — one row per unique variant

Every variant observed in any sample, deduplicated. Carries its full annotation set plus
a tier, modifier-candidate flag, and cohort-level carrier counts. **`tier_reason` is the most
important column for downstream interpretation** — sort by it.

Tier scale:

| Tier | Meaning |
|------|---------|
| **1** | Most likely pathogenic. Calibrated strong evidence appropriate for the variant's class. |
| **2** | Likely pathogenic. One strong line of class-appropriate evidence. |
| **3** | Uncertain, leaning damaging. Sub-threshold class-appropriate signal. Manual review warranted. |
| **4** | Uncertain — no information. Rare but no class-appropriate predictor is available or neutral. |
| **5** | Very likely benign. Predictor calls benign, OR fails rarity gate, OR ClinVar B/LB. |

`tier_reason` values that appear in this cohort, in rough order of evidence strength:

- `missense: AM likely_pathogenic + ESM1b damaging (orthogonal agreement) in constrained gene`
  — both an MSA-based model (AlphaMissense) and an MSA-free protein language model (ESM1b)
  independently call the missense damaging, in a gene the cohort treats as LoF-intolerant.
  Highest-confidence missense bucket.
- `pLoF in constrained gene; not NMD-escape` — predicted loss-of-function (stop-gain,
  frameshift, splice ±1/±2) in a constrained gene, outside the last exon. Standard
  high-confidence LoF.
- `missense: AM>=0.85 in constrained gene` — AlphaMissense at or above the strong-pathogenic
  threshold in a constrained gene.
- `pLoF in NF-spectrum tumor suppressor (NF1/NF2/SMARCB1/LZTR1)` — manual override for the
  four NF anchor genes; rare pLoF reaches Tier 1 regardless of cohort constraint metric. See
  `docs/pipeline_methods.md` §7.9 for the rationale (gnomAD v4 reports LZTR1 as unconstrained,
  contradicting its established schwannomatosis role).
- `missense: AM>=0.85 in NF-spectrum tumor suppressor` — same override, missense branch.
- `ClinVar P/LP (star filter disabled; see methods §7.10)` — ClinVar Pathogenic or
  Likely_pathogenic. **In this cohort the star filter is disabled** due to an upstream parser
  bug (now fixed in repo for next runs); this bucket therefore includes 0-star
  "no-assertion-criteria" submissions and should be down-weighted when clinical-grade
  confidence is required. Methods §7.10 has full detail.
- `missense: ESM1b damaging in constrained gene (AM missing)` (Tier 2) — ESM1b-only fallback
  for variants outside AlphaMissense coverage.
- `pLoF in non-constrained gene / NMD-escape` (Tier 2) — pLoF that doesn't meet Tier 1 bar.
- `missense: AM likely_pathogenic` (Tier 2) — AlphaMissense 0.564–0.85 without orthogonal
  ESM1b confirmation.
- `non-canonical splice + SpliceAI ...` (Tiers 1/2/3 by score band) — when SpliceAI is
  available; this cohort runs without SpliceAI (dropped from Stage 1 for throughput; see
  methods §7.2 and `docs/noncoding_v2_plan.md`).

In parallel with the tier, each variant carries `modifier_candidate` (boolean) plus
`modifier_evidence`. Modifier candidates are not Tier 1 candidates — they're the parallel
track for variants of moderate effect that might modulate phenotype severity rather than cause
disease. See methods §3.

### `cohort.genes.tsv` — one row per gene

Per-gene rollups: variant counts at each tier, modifier-candidate counts, descriptions, plus
sample-burden columns. Useful as a starting point for gene-level prioritization or for joining
with external gene-set analyses.

### `figures/` — 8 cohort-summary figures

Slide-ready PNGs at 300 dpi (10–13 in wide). All caption text is in-figure so they're
shareable without context.

| # | Figure | Purpose |
|---|--------|---------|
| 01 | `01_tier_distribution.png` | Overall variant counts per tier across cohort |
| 02 | `02_variant_class_by_tier.png` | Variant-class composition (pLoF / missense / etc.) within each tier |
| 03 | `03_per_sample_burden.png` | Distribution of Tier 1+2 variants per sample. Outliers (very high or zero) warrant QC review |
| 04 | `04_modifier_landscape.png` | Top 20 genes by modifier-candidate count |
| 05 | `05_nf_anchor_landscape.png` | **Positive control.** Tier breakdown for NF1, NF2, SMARCB1, LZTR1. Confirms the pipeline detects expected signal in established NF tumor-suppressor genes |
| 06 | `06_variant_recurrence.png` | Variant-count vs n_carriers (log-log). Shows what fraction of high-priority variants are private (carriers=1) vs recurrent |
| 07 | `07_top_recurrent_t12.png` | Top 30 most carrier-rich Tier 1+2 variants. Founder mutations vs QC artifacts in low-complexity regions |
| 08 | `08_sample_gene_hotspots.png` | (sample, gene) pairs with ≥2 distinct rare Tier 1+2 variants. Compound-het / multi-hit candidates. NF anchor genes always promoted to top. QC-filtered for mappability-poor families (ZNF*, KRT*, MUC*, OR*, HLA-*, etc.) |

---

## How to prioritize results

For a typical "find the strongest candidates" pass through `cohort.variants.tsv`:

1. **Start with NF anchor genes.** Filter to `gene ∈ {NF1, NF2, SMARCB1, LZTR1}` and `tier == 1`.
   These are your highest-confidence diagnostic-grade candidates. Fig 05 summarizes the same
   set graphically.
2. **Then the orthogonal-agreement missense bucket.** Filter to
   `tier_reason` containing `"orthogonal agreement"`. Two independent predictors agreed plus
   the gene is constrained — strong evidence even outside the anchors.
3. **Then pLoF in constrained genes.** `tier_reason` containing `"pLoF in constrained gene"`.
4. **Down-weight ClinVar-only Tier 1.** The ClinVar P/LP bucket is broader in this cohort than
   you'd usually want (star filter disabled). Manually re-check star ratings via ClinVar
   web or restrict to variants that ALSO hit a non-ClinVar Tier 1 path.
5. **Check compound-het candidates.** Fig 08 + `cohort.genes.tsv` `n_compound_het_candidates`.
   NF1 patients with two distinct rare Tier 1+2 NF1 variants are particularly interesting
   (and should be phase-checked via trios or long-read).
6. **For modifier-of-severity questions**, switch to the parallel track: filter
   `modifier_candidate == TRUE` and intersect with your phenotype-stratified subset.

---

## Known limitations to flag in interpretation

Brief summary; full detail in `docs/pipeline_methods.md` §7.

- **Non-coding variants are under-characterized** in this v1 run. SpliceAI was dropped from
  Stage 1 for throughput; conservation (PhyloP/GERP) and regulatory scoring aren't yet
  integrated. About 137,972 rare non-coding variants land in Tier 4 by default. v2 plan in
  `docs/noncoding_v2_plan.md`.
- **LOFTEE is not currently integrated** — pLoF Tier 1 calls have a known ~5–15% over-call
  rate from missing LOFTEE LC/HC filtering. NMD-escape heuristic partially compensates.
  Methods §7.3.
- **ClinVar star filter currently disabled.** The bucket includes 0-star submissions; the
  parser fix is in the repo and next runs will restore star filtering. Methods §7.10.
- **gnomAD v4 LOEUF / pLI shift.** v4's stricter LoF curation makes some classically
  constrained genes (notably LZTR1) appear unconstrained. We use an OR rule
  (LOEUF < 0.35 OR pLI ≥ 0.9) and an NF-anchor override to compensate. Methods §7.9.
- **Mappability-poor gene families are silently filtered** out of the sample-gene
  hotspot figure (ZNF*, KRT*, MUC*, OR*, HLA-*, etc.) because short-read alignments
  routinely mis-map between paralogs. They are NOT filtered out of the main tier
  table — be cautious about variants in those families in the variant table. See
  `QC_ARTIFACT_GENE_PATTERNS` in `R/cohort_figures.R`.
- **AF filter is single-population popmax**, applied at 1×10⁻⁴ for tier rarity and
  variants without a gnomAD entry pass through. Methods §7.5 / 7.8.

---

## Re-running the pipeline

If you received only the outputs and want to re-run on a different cohort, or modify
parameters and re-tier:

- **Just want to re-tier with different thresholds?** Edit the `CFG` block at the top of
  `R/tier_variants.R` and re-run `Rscript R/tier_variants.R --input <annotated.tab>
  --out-prefix <prefix>`. No fastVEP run needed.
- **Want to re-run from raw VCFs?** See `INSTALL.md` for dependencies, fastVEP setup, and
  cohort manifest format.
- **Want to re-render figures only?** `Rscript R/cohort_figures.R --variants <variants.tsv>
  --genes <genes.tsv> --genotypes <genotypes.tsv> --out-dir <dir>`.
- **Want to validate that the outputs are internally consistent?**
  `scripts/validate_cohort.sh --out-dir <dir>` runs 9 sections of structural and rule checks.

---

## Document map

| Document | What it covers |
|----------|----------------|
| **`README.md`** *(this file)* | Quick orientation: what the outputs are, how to prioritize results, limitations summary, doc map |
| **`docs/results_guide.md`** | **Reading reference.** Full column dictionary for every output file, complete tiering rules table (all 23 rules with conditions), reading recipes for common questions. Start here when reading the data |
| **`docs/pipeline_methods.md`** | Comprehensive methods: data sources, every tier rule with rationale, validation, known gotchas. The reference document for collaborators publishing results from this pipeline |
| **`docs/tiering.md`** | Original tiering-scheme design spec. Higher-level conceptual reference for the class-conditional logic |
| **`docs/noncoding_v2_plan.md`** | Roadmap for v2 non-coding annotation (SpliceAI restore, PhyloP/GERP, regulatory, LOFTEE) |
| **`docs/dbnsfp_score_columns.md`** | Technical reference for the dbNSFP column extraction used by the fastVEP fork |
| **`INSTALL.md`** | Dependencies, fastVEP fork setup, cohort manifest format, environment requirements, integration steps. Read this if you're running the pipeline rather than interpreting its outputs |
| **`scripts/validate_cohort.sh`** | Self-documenting validation: read the file header for the 9 check sections |

---

## Questions / issues

This repo is the working code for a specific NF1 study; if you find a bug, ambiguity, or
have a suggestion, please open an issue or reach out directly. For interpretation questions
specifically about Tier 1 candidates in your samples, the right starting point is usually the
`tier_reason` column plus the methods doc section that defines that tier_reason.
