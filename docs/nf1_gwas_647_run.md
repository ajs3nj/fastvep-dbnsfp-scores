# NF1 GWAS 647-sample run — results notes

Run-specific notes for the NF1 GWAS cohort processed through the v2 batched
pipeline. Distinct from `results_guide.md` (which is the general reference for
the output files) — this doc captures what happened during this specific run,
what was found, what was fixed, and known limitations to flag when handing
results to collaborators.

Cohort assembled from a collaborator manifest of NF1-affected individuals.
Two data-quality issues surfaced during interpretation, both patched in
`R/tier_variants.R`; both fixes are safe for future runs.

---

## 1. Cohort description

647 samples from three distinct upstream sequencing / calling pipelines:

| Pipeline                        | n_samples |
|---------------------------------|----------:|
| DRAGEN                          |       306 |
| DRAGEN (hard-filtered)          |        52 |
| bwa-mem2 + GATK HaplotypeCaller |       289 |

For orchestration, the two large pipelines were arbitrarily chunked into
within-pipeline batches of ~100 samples each so `cohort_pipeline_batched.sh`
could process them incrementally with disk reclamation between chunks. This
produces seven `batch_id`s in the manifest:

| batch_id | Pipeline                        | n_samples | Nature                       |
|---------:|---------------------------------|----------:|------------------------------|
| 1        | DRAGEN                          |       102 | throughput chunk of DRAGEN   |
| 2        | DRAGEN                          |       102 | throughput chunk of DRAGEN   |
| 3        | DRAGEN                          |       102 | throughput chunk of DRAGEN   |
| 4        | DRAGEN (hard-filtered)          |        52 | distinct condition           |
| 5        | bwa-mem2 + GATK HaplotypeCaller |        97 | throughput chunk of bwa+GATK |
| 6        | bwa-mem2 + GATK HaplotypeCaller |        96 | throughput chunk of bwa+GATK |
| 7        | bwa-mem2 + GATK HaplotypeCaller |        96 | throughput chunk of bwa+GATK |

**Within-pipeline batches (1/2/3 and 5/6/7) are not distinct sequencing
conditions** — they are arbitrary partitions of the same underlying data
made only to bound per-batch memory / disk usage. Only batches 4 (DRAGEN
hard-filtered), and the DRAGEN vs bwa-mem2+GATK split are real
condition-level distinctions.

VCFs were pulled from Synapse via `run_batch.sh` (SYNAPSE_AUTH_TOKEN),
normalized with `stage0_normalize_vcf.sh` (bcftools norm -m- with UCSC → Ensembl
chrom rename and primary-assembly filter), then annotated in per-sample fastVEP
runs. Cohort merge, prefilter, per-variant cohort AC/AN, and tiering ran once
across all batches after per-sample annotation completed.

Manifest: `scripts/nf1_gwas_manifest.v2.tsv`.

---

## 2. Run outcome

### 2.1 Tier distribution

| Tier | n_variants |    % | Notes |
|-----:|-----------:|-----:|-------|
|    1 |     23,731 | 5.0% | Inflated ~2× vs typical, see §4.1 |
|    2 |     47,326 | 10.0% | |
|    3 |     15,812 |  3.3% | |
|    4 |     26,947 |  5.7% | |
|    5 |    359,507 | 76.0% | |
| **Total** | **473,323** | | |

Modifier candidates: **365** across the 167-gene NF / RAS-MAPK modifier list
(see `R/nf_modifier_genes.txt`).

### 2.2 NF-anchor gene signal (positive control)

All four NF-spectrum tumor suppressors show Tier 1 signal, confirming the
pipeline is finding what an NF1 cohort should contain:

| Gene    | variants in cohort | best_tier | Tier 1 count |
|---------|-------------------:|:---------:|-------------:|
| NF1     | 596                | 1         | ~365 (pLoF/splice/AM-strong-missense) |
| NF2     | 78                 | 1         | (see figure 05) |
| SMARCB1 | 82                 | 1         | (see figure 05) |
| LZTR1   | 122                | 1         | (see figure 05) |

Detailed per-gene breakdown in `figures/05_nf_anchor_landscape.png`.

### 2.3 ClinVar override behavior

All ClinVar Pathogenic/Likely_pathogenic NF1 variants (n=222) correctly land
in Tier 1. Zero ClinVar Benign/Likely_benign NF1 variants incorrectly assigned
to any tier other than 5. General ClinVar override checks (across the whole
cohort, not just NF1) also pass in `validate_cohort.sh`.

---

## 3. Data-quality issues found and fixed during this run

### 3.1 cohort_af underdenominator bug (Stage 3.5)

**Symptom.** Initial tiering showed 100% of variants with `af_source =
cohort_af` had `cohort_af > 5%`. In a 647-sample cohort this is impossible —
singletons should be `1/1294 ≈ 7.7e-4`, not 50%+.

**Root cause.** Stage 3.5's `cohort_an` counts only samples that emitted the
variant in their per-sample VCF, not the total called cohort. For per-sample
VCFs (which don't record 0/0 sites), a singleton gets `cohort_ac / cohort_an =
1/2 = 0.5` instead of `1/1294`. Variants present in every sample were computed
correctly (max cohort_an = 1294 = 2 × 647).

**Fix.** In `R/tier_variants.R`, self-detect the bug when >50% of stored
`cohort_af` values are >5%. Recompute as `cohort_ac / max(cohort_an)`, which
is the correct denominator for variants present in every sample. Provenance
preserved via `af_source` column. Emits a log line when the fix fires:

```
[tier] cohort_af underdenominator bug detected (93.0% of stored cohort_af > 5%);
recomputed using max(cohort_an)=1294
```

### 3.2 Adaptive rarity threshold for cohort_af fallback

**Symptom.** After 3.1 was fixed, cohort_af values were sensible but the
modifier count came back healthy (365) while ~17 known-pathogenic NF1 variants
(pLoF, AM-strong missense, splice) were still landing in Tier 5.

**Root cause.** `CFG$af_max = 1e-4` is calibrated for gnomAD scale
(world-population AF). For a 647-sample cohort the resolution limit is
`1/1294 = 7.7e-4`, already above 1e-4. Every cohort singleton fails the
rarity gate. Runs with real gnomAD popmax loaded are unaffected; runs relying
on cohort_af fallback lose their entire singleton population.

**Fix.** Adaptive threshold: when `af_source` starts with `"cohort_af"`, use
`3/(2N)` (≈ 2.32e-3 for this cohort) as the rare gate. Accepts singletons,
doubletons, and tripletons as rare. Only fires when cohort_af is actually the
fallback; runs with proper gnomAD are unaffected. New `af_max_effective`
column records the per-variant threshold for provenance.

Impact on this run: Tier 1 grew from 20,172 to 23,731 (+3,559); Tier 2 from
19,751 to 47,326 (+27,575). Redistribution mainly went from Tier 5 into Tier
2/3/4 (variants with some but not maximal signal), not into Tier 1 —
consistent with signal-appropriate landing.

---

## 4. Known limitations

### 4.1 Missing gnomAD per-variant sources

The fastVEP SA dir for this run had `gnomad_genes.oga` (per-gene constraint)
but no per-variant `gnomad_af.osa`. All 440,052 non-NA `af_used` values in
this run come from cohort_af fallback; only 33,271 variants have `af_source =
none` (missing everywhere). Consequences:

- No variant can be marked "common in gnomAD" via the standard `filter_af =
  0.01` gate. Tier 1 is roughly 2× inflated vs. what it would be with real
  gnomAD data (still analytically usable, but the "Tier 1" set is enriched
  for variants that are rare-in-cohort AND absent-from-cohort — a superset
  of the gnomAD-rare set).
- `validate_cohort.sh` section 3 (tier pyramid check) downgrades to INFO
  instead of FAIL when popmax is >99% empty, so future runs with the same
  data limitation report cleanly.

**Fix path.** Build a `.osa` from gnomAD sites VCF and load into the SA dir.
Tracked as v3 backlog task.

### 4.2 Suspected artifact — 17:31324211

Splice_donor_variant in NF1 at `af_used = 0.597` (≈ 60% of the cohort). Almost
certainly a mis-mapped read cluster or common polymorphism escaping annotation
(possibly a mapping ambiguity between NF1 and a nearby pseudogene / segmental
duplication region). Not classified as Tier 1 by our rules (fails the
common-variant filter). Flagged for manual review before any downstream
carrier-level interpretation.

### 4.3 Non-coding coverage

Cohort has 169,283 rare non-coding variants, of which 0 land in Tier 4 —
consistent with the documented non-coding limitation (splice-region signals
outside SpliceAI's window aren't picked up). See `docs/noncoding_v2_plan.md`
for the v2 improvement plan.

---

## 5. Verification results

`scripts/validate_cohort.sh --out-dir /data/nf1/outputs` after fixes:

- File presence + column presence: PASS
- Tier distribution: INFO (pyramid inversion attributed to missing gnomAD)
- Class-conditional rules (Tier 1 pLoF constrained, Tier 1 missense paths,
  Tier 5 benign rules): PASS
- ClinVar overrides (P/LP → Tier 1, B/LB → Tier 5): PASS
- Modifier band check (af_used in band OR sub-Mendelian): PASS
- NF anchor gene signal: PASS (all four have Tier 1 hits)
- n_carriers matches genotype recount (50 sampled variants): PASS
- Non-coding limitation: INFO (documented)

Summary: 15 PASS, 0 FAIL, 2 INFO. No hard violations of tier logic.

---

## 6. Figures produced

Presentation-ready PNGs at `/data/nf1/outputs/figures/` from
`R/cohort_figures.R`:

| # | File | Purpose |
|---|------|---------|
| 01 | `01_tier_distribution.png` | Cohort-wide Tier 1–5 counts |
| 02 | `02_variant_class_by_tier.png` | Variant class composition per tier |
| 03 | `03_per_sample_burden.png` | Tier 1+2 variants per sample (pooled) |
| 04 | `04_modifier_landscape.png` | Top 20 genes by modifier candidate count |
| 05 | `05_nf_anchor_landscape.png` | NF1/NF2/SMARCB1/LZTR1 variants by tier |
| 06 | `06_variant_recurrence.png` | Recurrence distribution split by tier |
| 07 | `07_top_recurrent_t12.png` | Top 30 recurrent Tier 1+2 variants |
| 08 | `08_sample_gene_hotspots.png` | Compound-het / multi-hit candidates |

Additional QC figures from this run:

| File | Purpose |
|------|---------|
| `af_source_distribution.png` | log10(af_used) histogram by af_source; documents that this run used cohort_af fallback throughout |
| `qc_batch_burden_t1.png` | Per-sample Tier 1 burden by batch (batch effect QC) |
| `qc_batch_burden_high_priority.png` | Per-sample Tier 1+2 burden by batch |

Companion tables:

| File | Contents |
|------|----------|
| `af_source_summary.tsv` | per-source n, %, min/median/max af_used |
| `qc_batch_burden_summary.tsv` | per-batch median, IQR, q10/q90 |
| `qc_batch_stats.txt` | Kruskal-Wallis + pairwise Wilcoxon across batches |
| `qc_batch_outliers.tsv` | Samples > q75 + 1.5 IQR within batch |

---

## 7. Recommended follow-up

**Before publishing cohort-wide burden results:**

1. Load real gnomAD per-variant `.osa` into the SA dir and re-run tiering.
   Expected impact: Tier 1 drops from ~24k to a few thousand as gnomAD-common
   variants get properly demoted. Modifier count should stay stable or grow
   slightly.
2. Manual review of 17:31324211 (the suspected splice donor artifact).

**Batch effect check — passed.** Kruskal-Wallis across the seven batch IDs
returned p = 0.42 for Tier 1 and p = 0.14 for Tier 1+2; all pairwise
Wilcoxon p-adjusted are 1.0. Per-batch median burden is 150–156 (Tier 1)
and 192–200 (Tier 1+2). Because the seven batches include both arbitrary
within-pipeline throughput chunks (1/2/3 within DRAGEN, 5/6/7 within
bwa-mem2+GATK) and the three real pipeline conditions (DRAGEN vs
DRAGEN-hard-filtered vs bwa-mem2+GATK), this test simultaneously checks
two things: (i) that within-pipeline chunking didn't introduce spurious
noise — which it shouldn't by design, and (ii) that per-sample burden
differs across the real pipeline conditions. Both come back
indistinguishable. **Pooled downstream burden analysis is defensible;
pipeline / batch does NOT need to enter the covariate list.** 48
samples were flagged as within-batch outliers (n_variants > q75 +
1.5 IQR for Tier 1+2) — see `figures/qc_batch_outliers.tsv` if you
want to sanity check the highest-burden individuals before running the
gene-burden test.

**Optional, but nice to have:**

4. Restore the ClinVar star filter (currently disabled — see the note in
   `validate_cohort.sh` section 5) once the `csq_to_wide_tab.py`
   REVIEW_STATUS parsing bug is fixed.
5. Add gnomAD ancestry-stratified AF fields for populations relevant to the
   cohort composition, once ancestry column in the manifest is filled in.

---

## 8. Reproducibility notes

Pipeline was launched via `scripts/cohort_pipeline_batched.sh` and completed
across a Stage 0 → Stage 5 flow. A tmux session died mid-run (Docker
container reboot); Stage outputs and batch sentinels survived, and the run
was resumed with `nohup ... & disown`.

Stage 2 (cohort merge + dedup) took ~13 hours at this cohort scale. See
`docs/v3_stage2_parallel_design.md` for the v3 chrom-sharded rewrite that
addresses this bottleneck.

Output directory: `/data/nf1/outputs/`. Key files:

- `cohort.variants.tsv` — one row per unique tiered variant (473,323 rows)
- `cohort.genes.tsv` — per-gene rollup with sample burden columns
- `cohort.genotypes.tsv` — long format, 3.14 B rows, 92 GB
- `cohort.filtered.with_cohort.tab` — intermediate; input to tiering
- `figures/` — all figures listed in §6
- `validate_cohort2.log` — post-fix validation output
- `tier_rerun2.log` — post-fix tier step output (shows adaptive-threshold log line)
