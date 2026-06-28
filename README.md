# fastvep-dbnsfp-scores

A variant annotation and tiering pipeline for rare-disease cohorts, built and validated on
a 209-sample neurofibromatosis type 1 (NF1) WGS study. The pipeline annotates variants
with calibrated missense pathogenicity predictors (AlphaMissense, ESM1b, REVEL) using a
patched fastVEP fork, assigns class-conditional priority tiers, flags phenotype-modifier
candidates as a parallel track, and produces per-variant + per-gene cohort summaries with
slide-ready figures.

---

## What it does

Given a manifest of per-sample VCFs:

1. Runs fastVEP per sample to attach VEP consequences plus dbNSFP scores (AlphaMissense,
   ESM1b, REVEL alongside SIFT/PolyPhen) in a single annotation pass.
2. Merges per-sample annotations into a deduplicated cohort table with cohort allele
   counts and carrier statistics.
3. Pre-filters out common and no-signal variants.
4. Assigns each surviving variant a tier (1–5) using class-conditional rules — pLoF,
   missense, splice, in-frame, non-coding, and regulatory each get an evidence path
   appropriate to that class. An independent `modifier_candidate` flag runs in parallel.
5. Rolls up per-gene burden counts and emits clinical-presentation figures.

---

## What you get out

| Output | Contents |
|--------|----------|
| `cohort.variants.tsv` | One row per unique variant. Annotation + tier + modifier flag + cohort carrier statistics |
| `cohort.genes.tsv` | One row per gene. Variant counts by tier, modifier counts, sample-level burden |
| `cohort.genotypes.tsv` | Long-format per-(variant, sample) carrier table |
| `figures/01–08.png` | Tier distribution, per-sample burden, modifier landscape, NF anchor positive control, variant recurrence, top recurrent variants, sample-gene hotspots |
| `_provenance.log` | Git SHA + parameter values + reference data paths for the run |

For the **complete column-by-column dictionary**, the **tiering rules table**, and
**reading recipes**, open `docs/results_guide.md` alongside the output.

---

## Running the pipeline

See **`INSTALL.md`** for dependencies, the fastVEP fork build, reference data, manifest
format, full `cohort_pipeline.sh` invocation with all flags, per-stage runtime expectations,
common re-run / re-figure / validate workflows, and troubleshooting.

---

## Known limitations

Brief; full detail in `docs/pipeline_methods.md` §7.

- Non-coding variants under-characterized in v1 (SpliceAI dropped from Stage 1 for
  throughput; conservation and regulatory scoring not yet integrated). v2 plan in
  `docs/noncoding_v2_plan.md`.
- LOFTEE not integrated — pLoF Tier 1 calls have a known ~5–15% over-call rate.
- ClinVar star filter currently disabled in the v1 NF1 cohort (parser bug, now fixed in
  repo); ClinVar P/LP Tier 1 bucket includes 0-star submissions.
- gnomAD v4 LOEUF / pLI shift mis-classifies LZTR1 as unconstrained; an NF-anchor
  override compensates for the four anchor genes.
- Mappability-poor gene families (ZNF*, KRT*, MUC*, OR*, HLA-*) filtered from the
  hotspot figure but NOT from the main tier table.

---

## Document map

| Document | What it covers |
|----------|----------------|
| **`README.md`** *(this file)* | What the pipeline is and what it produces |
| **`INSTALL.md`** | Everything operational: dependencies, build, manifest format, running the pipeline, troubleshooting |
| **`docs/results_guide.md`** | Reading reference: tiering rules table, column dictionary for every output file, reading recipes |
| **`docs/pipeline_methods.md`** | Comprehensive methods: data sources, every rule with rationale, validation, known gotchas |
| **`docs/tiering.md`** | Original tiering-scheme design spec |
| **`docs/noncoding_v2_plan.md`** | Roadmap for v2 non-coding annotation |
| **`docs/dbnsfp_score_columns.md`** | Technical reference for the dbNSFP column extraction |
| **`INTEGRATION.md`** | fastVEP fork integration steps (Rust-side) |

---

## Questions / issues

This repo is the working code for a specific NF1 study; if you find a bug, ambiguity, or
have a suggestion, please open an issue or reach out directly.
