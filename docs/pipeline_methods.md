# Cohort variant annotation and tiering pipeline

## Methods documentation for collaborators

This document describes the methods used to annotate, filter, and prioritize
genomic variants in cohort-scale rare-disease studies, with an immediate focus
on neurofibromatosis (NF) cohort modifier discovery. It is intended for
clinical and research collaborators who need to understand and interpret the
pipeline outputs, whether or not they run the pipeline themselves.

The pipeline is open-source and lives at
`github.com/<owner>/fastvep-dbnsfp-scores`. The reproducible artifact for any
given run is the commit hash of that repo plus the manifest of input VCFs,
reference data versions, and command-line flags used.

---

## 1. Overview

### 1.1 Purpose

For each variant carried by any sample in a cohort, the pipeline produces a
**tier** (1 = most likely pathogenic, 5 = very likely benign) and a parallel
**modifier candidate flag** (TRUE if the variant plausibly modulates phenotype
rather than causes disease). Per-gene rollups aggregate variant-level signal
to support both Mendelian candidate hunting and modifier-discovery hypotheses.

The pipeline is designed for **germline rare-disease research prioritization
in small to mid-size cohorts** (~50–500 samples). It is not a clinical-grade
variant classification system (it does not produce ACMG-compliant
pathogenicity assertions by default — though an `--mode acmg` knob exists in
the tiering R script for that use case), and it is not optimized for somatic
or large-cohort GWAS workflows.

### 1.2 Inputs

- One per-sample VCF per cohort sample (GRCh38 reference, single-sample or
  joint-called and split). Chromosome naming should be Ensembl-style (no `chr`
  prefix) to match the transcript cache; the pipeline includes a chromosome
  rename utility for VCFs using UCSC-style naming.
- A two-column TSV **manifest** listing `sample_id` and `vcf_path`. Additional
  metadata columns (phenotype, severity, batch, sex) are passed through to a
  side file but do not affect tiering.

### 1.3 Outputs

- `cohort.variants.tsv` — one row per unique variant in the cohort, with
  annotation, tier, modifier flag, and cohort-level carrier counts.
- `cohort.genes.tsv` — one row per gene with rollup statistics (variant
  counts by tier, modifier candidate count, sample carrier counts, gene
  description).
- `cohort.genotypes.tsv` — long-format genotype table (one row per
  variant×sample). Used internally for cohort-level burden; retained for
  downstream burden testing.

### 1.4 Pipeline stages

```
manifest          ┌──────────────────┐
   │              │ Stage 1:         │   per-sample
   │              │ fastvep annotate │   .annotated.tab.gz +
   ▼              │ + CSQ converter  │   .genotypes.tsv
per-sample        │ + bcftools query │
  VCFs            └────────┬─────────┘
                           │
                  ┌────────▼─────────┐
                  │ Stage 2:         │   cohort.annotated.tab
                  │ cohort merge     │   cohort.genotypes.tsv
                  │ + n_samples join │   _variant_counts.tsv
                  └────────┬─────────┘
                           │
                  ┌────────▼─────────┐
                  │ Stage 3:         │   cohort.filtered.tab
                  │ pre-filter       │
                  └────────┬─────────┘
                           │
                  ┌────────▼─────────┐
                  │ Stage 3.5:       │   cohort.variant_summary.tsv
                  │ aggregate gts    │   cohort.filtered.with_cohort.tab
                  │ + join cohort    │
                  └────────┬─────────┘
                           │
                  ┌────────▼─────────┐
                  │ Stage 4:         │   cohort.variants.tsv
                  │ tier_variants.R  │   cohort.genes.tsv
                  └────────┬─────────┘
                           │
                  ┌────────▼─────────┐
                  │ Stage 5:         │   cohort.gene_burden.tsv
                  │ per-gene burden  │   cohort.genes.tsv (merged)
                  └──────────────────┘
```

---

## 2. Per-stage description

### 2.1 Stage 1 — per-sample annotation

For each sample VCF, three things happen in parallel across cores:

1. **fastvep annotate**, producing a VCF with the standard 48-field VEP CSQ
   INFO field plus per-allele `FV_*` INFO fields for each loaded supplementary
   annotation database (`FV_DBNSFP`, `FV_CLINVAR`, `FV_GNOMAD`, `FV_GNOMAD_GENE`,
   etc.). Reference data: GRCh38 Ensembl release 115 transcript cache,
   GRCh38 primary-assembly FASTA with `.fai` index.

2. **CSQ-to-wide-tab conversion** (`scripts/csq_to_wide_tab.py`), which parses
   the VCF, splits each CSQ entry and each `FV_*` pipe-delimited projection
   into named columns, and writes a tab-separated file with one row per
   `(variant, transcript)` consequence. Output is gzipped (~80% smaller than
   plain TSV) and includes the columns the downstream tiering expects.

3. **Genotype extraction** via `bcftools query` into a long-format TSV
   (`chrom, pos, ref, alt, sample, gt`), suitable for cohort-level aggregation
   downstream.

Per-sample work is embarrassingly parallel via `xargs -P` (or GNU parallel if
available); per-sample failures don't block other samples.

### 2.2 Stage 2 — cohort merge

The 209-sample (or N-sample) per-sample tabs are merged into:

- `cohort.annotated.tab` — one row per unique `(chrom, pos, ref, alt, transcript)`
  across the cohort, with an added `n_samples_annotated` column reporting how
  many cohort samples have this variant in their annotation table.

- `cohort.genotypes.tsv` — concatenation of all per-sample genotype TSVs.

The merge is implemented as a gzip-aware streaming awk pass using a sentinel
marker between files; per-sample uniqueness is tracked in a hash that resets
between files, while cohort-wide counts persist. Memory peaks at ~5–10 GB on
a 200-sample WGS cohort.

### 2.3 Stage 3 — pre-filter

The merged cohort table contains tens to hundreds of millions of rows. The
pre-filter drops rows that have no chance of being meaningful candidates,
keeping the input to downstream stages tractable. A variant survives if **any
one of these is true**:

1. `IMPACT == HIGH` or `MODERATE` (coding/splice consequence)
2. `spliceai_ds_max >= 0.2` (cryptic splice signal — when SpliceAI is loaded)
3. `clin_sig` non-empty (ClinVar opinion exists)
4. `gnomad_popmax_af <= filter_af` (default 0.01 — explicitly rare)

By default, **missing gnomAD popmax does not bypass the rarity gate** — a
non-coding variant with no functional signal and no explicit rarity evidence
gets dropped. A `--filter-allow-missing-popmax` escape hatch restores the loose
behavior. Earlier pipeline versions used the loose rule and produced cohort
tables 10–20× larger that exceeded R's memory budget.

By default, chromosomes X and Y are also dropped at this stage. A
`--include-sex-chroms` flag overrides this for sex-aware analyses.

### 2.4 Stage 3.5 — cohort genotype aggregation

The 30+ GB long-format `cohort.genotypes.tsv` is condensed via streaming awk
into `cohort.variant_summary.tsv` (~5 GB), which has one row per unique
variant with these aggregate columns:

- `cohort_ac` — alt allele count across cohort
- `cohort_an` — allele number (2 × number of called samples)
- `cohort_af` — internal cohort allele frequency (`cohort_ac / cohort_an`)
- `n_het` — number of heterozygous carriers
- `n_hom` — number of homozygous-alt carriers
- `n_carriers` — `n_het + n_hom` (samples with at least one alt allele)

These six columns are then joined onto the pre-filtered table to produce
`cohort.filtered.with_cohort.tab`. This avoids the R tiering script having to
fread the full long-format genotype table (which OOMs on a typical cohort).

### 2.5 Stage 4 — variant classification (tiering)

`R/tier_variants.R` reads `cohort.filtered.with_cohort.tab`, assigns each
variant a tier (1–5) and a modifier candidate flag, computes a per-gene
summary, and writes `cohort.variants.tsv` + `cohort.genes.tsv`. The tiering
logic is described in detail in §4 below, and the design rationale is in
`docs/tiering.md`.

### 2.6 Stage 5 — per-gene sample burden

The per-gene rollup gets supplemented with sample carrier counts via awk over
the long-format genotype table joined against the tiered variants:

- `n_samples_high_impact` — distinct samples carrying ≥1 HIGH-impact variant
  in this gene
- `n_samples_tier1` — distinct samples carrying ≥1 Tier 1 variant
- `n_samples_tier12` — Tier 1 or 2
- `n_samples_qualifying` — Tier 1 or 2, OR any pLoF at Tier ≤3 (catches
  NMD-escape / LOFTEE-LC pLoFs demoted from Tier 1)
- `n_samples_with_modifier_candidate` — distinct samples carrying ≥1
  modifier-flag variant

These columns are merged into `cohort.genes.tsv` in place.

---

## 3. Data sources

| Source | Version | Use | Notes |
|---|---|---|---|
| Reference FASTA | GRCh38 Ensembl primary assembly | Reference allele check, HGVS | `.fai` index required for fast memory-mapped access |
| Gene/transcript models | Ensembl release 115 | Consequence prediction, transcript-aware HGVS | Built into a binary fastVEP transcript cache once; cached for all subsequent runs |
| dbNSFP | v5.3.1a (Nov 2024) | AlphaMissense, ESM1b, REVEL per-variant missense scores | Pulled via our extended `--source dbnsfp` parser; provides scores via a single SA lookup. SIFT and PolyPhen also extracted. |
| SpliceAI | Ensembl MANE-restricted v1.4 | Splice-altering signal for SNVs in/near MANE Select exons | **Currently dropped from Stage 1** due to throughput issues (40 GB `.osa`, 5+ min per chr22-scale annotation). Planned post-tier reattachment per `docs/noncoding_v2_plan.md`. |
| ClinVar | 2025 release | Pathogenicity assertions, review status | Used as a hard override on tier assignment (P/LP at ≥1 star → Tier 1; B/LB → Tier 5) |
| gnomAD v4.1 | Constraint metrics | LOEUF, pLI, mis_z, syn_z per gene | Drives the "constrained gene" criterion in tiering rules |
| gnomAD v4.1 | Per-variant AF | popmax, per-population AF, AC, AN | Drives rarity filtering at all stages |
| CADD | not used | (would provide genome-wide deleteriousness) | Deliberately excluded for v1 (delivery and integration cost); see `docs/noncoding_v2_plan.md` for status |

Reference data is downloaded once and built into local fastSA databases
(`.osa` / `.oga`) once per release. The reproducible setup commands are in
the project README.

---

## 4. Filtering and tiering rules

### 4.1 Pre-filter (Stage 3)

A variant is kept for tiering if any of:

| Rule | Reason |
|---|---|
| IMPACT in {HIGH, MODERATE} | Coding/splice change — always potentially relevant |
| SpliceAI ds_max ≥ 0.2 | Cryptic splice signal (when loaded) |
| ClinVar significance non-empty | Existing clinical opinion warrants review |
| gnomAD popmax ≤ 1e-2 | Explicitly rare in the population |

The pre-filter is conservative on the "rare" path: missing popmax alone is
not enough. Variants without functional signal AND without explicit rarity
evidence are dropped at this stage.

### 4.2 Tier definitions

Tiers reflect pathogenicity likelihood. The class-conditional logic in §4.3
operationalizes them per variant class; the *meaning* of each tier is fixed:

| Tier | Meaning |
|---|---|
| **1** | Most likely pathogenic. Strong calibrated evidence appropriate for the class; ClinVar P/LP at ≥1 star is always Tier 1. |
| **2** | Likely pathogenic. One strong calibrated line of class-appropriate evidence. |
| **3** | Uncertain, leaning damaging. Suggestive sub-threshold signal. Manual review warranted. |
| **4** | Uncertain — no information. Rare, but the class-appropriate predictor is missing, not applicable, or neutral. The honest "we don't know" bucket. |
| **5** | Very likely benign. Class-appropriate evidence calls benign, or fails rarity gate, or ClinVar B/LB. |

### 4.3 Class-conditional rules

Variants are classified by their *most severe consequence* on the
representative transcript (MANE Select → canonical → most-severe SO term).
One primary predictor per class — no cross-class consensus stacking, which
would systematically penalize variant classes for lacking predictors that
physically don't apply to them.

| Class | Tier 1 | Tier 2 | Tier 3 | Tier 4 | Tier 5 |
|---|---|---|---|---|---|
| **pLoF** (stop_gained, frameshift, splice_acceptor, splice_donor, start_lost, transcript_ablation) | rare + constrained gene (LOEUF<0.35 or pLI≥0.9), not NMD-escape, LOFTEE HC if present | rare LoF in unconstrained gene; or NMD-escape; or LOFTEE LC | — | — | common, or ClinVar B/LB |
| **non-canonical splice** | rare + SpliceAI ≥ 0.8 | rare + SpliceAI 0.5–0.8 | rare + SpliceAI 0.2–0.5 | rare, no SpliceAI signal | common |
| **missense / protein_altering** | ClinVar P/LP ≥1 star; or rare + AlphaMissense ≥ 0.85 in constrained gene | rare + AlphaMissense likely_pathogenic (>0.564) | rare + AlphaMissense ambiguous (0.34–0.564); or rare missense in constrained gene with AM missing | rare, AM neutral or missing | AM likely_benign (<0.34), or common, or ClinVar B/LB |
| **in-frame indel** | — | rare + LOEUF<0.35; or SpliceAI ≥0.5 | rare in moderately constrained gene; or SpliceAI 0.2–0.5 | rare, no constraint | common |
| **synonymous / UTR / intronic / non-coding tx** | — | rare + SpliceAI ≥0.5 (cryptic splice) | rare + SpliceAI 0.2–0.5 | rare, no SpliceAI signal | common |
| **regulatory** (motif/TFBS, regulatory_region_variant) | — | rare + motif disruption + strong SpliceAI | motif/TFBS overlap + rare | rare, no signal | common |

Two override rules applied after class-conditional assignment:
- ClinVar P/LP at ≥1 star → Tier 1 (regardless of computational signal)
- ClinVar B/LB → Tier 5 (regardless of computational signal)

The rarity gate (default popmax ≤ 1e-4) is the prerequisite for any
non-Tier-5 assignment. A variant failing rarity is always Tier 5.

### 4.4 Modifier candidate flag

Parallel to the tier, each variant gets a modifier candidate flag (TRUE /
FALSE). The flag asks a different question: "given that this variant likely
isn't causing disease on its own, could it modulate phenotype?" Modifiers
tend to be common-to-low-frequency variants of moderate effect, often in
downstream pathway genes.

A variant is `modifier_candidate = TRUE` if **all** hold:

1. **AF in modifier band** — either gnomAD popmax `1e-4 < AF ≤ 5e-2`, or
   rare-but-sub-Mendelian (rare + Tier 3 or 4).
2. **Functional signal appropriate for class** — at least one of:
   - missense with `0.34 < AlphaMissense ≤ 0.85` (moderate)
   - non-canonical-splice or non-coding with `0.2 ≤ SpliceAI < 0.5`
   - any pLoF (LoF carriers are biologically plausible for dose-sensitive
     pathway components)
   - regulatory_region_variant or TFBS overlap
3. **Not failed-common** (popmax > 5%). Above 5% it's a polymorphism.

The same variant can be `tier = 4` AND `modifier_candidate = TRUE` — that
combination is exactly what a real modifier looks like (rare or low-freq
with no Mendelian-strength signal, but with sub-Mendelian functional
evidence in a relevant gene).

A starter gene list for the modifier flag is in `R/nf_modifier_genes.txt` —
NF1, NF2, SMARCB1, LZTR1 as anchors, then their direct interactors,
canonical RAS/MAPK and Hippo pathway members, and published NF modifier
loci (CDKN2A/B, ATM, TP53, EED, SUZ12). Pass via `--modifier-genes` to
restrict the modifier flag to a curated set; default is genome-wide.

---

## 5. Output schemas

### 5.1 cohort.variants.tsv (one row per unique variant)

Identity:
- `chrom`, `pos`, `ref`, `alt` — variant key (GRCh38, Ensembl-style chrom naming)
- `gene`, `gene_id` — HGNC symbol and Ensembl gene ID of the representative transcript
- `transcript` — representative transcript ID (MANE Select / canonical / most severe consequence)
- `hgvsc`, `hgvsp` — HGVS coding and protein notation on the representative transcript
- `existing_variation` — RS / ClinVar IDs if present
- `mane_select`, `canonical`, `mane_plus_clinical` — transcript-selection flags

Consequence and impact:
- `consequence` — most-severe SO term on the representative transcript
- `impact` — VEP impact rating (HIGH/MODERATE/LOW/MODIFIER)
- `variant_class` — pipeline-assigned class (pLoF / splice_noncanonical / missense / inframe / noncoding / regulatory / other)
- `exon`, `intron` — exon/intron position (`i/N` format) on representative transcript
- `biotype` — protein_coding, lncRNA, etc.

Functional scores:
- `alphamissense`, `am_class` — AlphaMissense score and class (likely_pathogenic / ambiguous / likely_benign)
- `esm1b`, `esm1b_class` — ESM1b log-likelihood ratio (negative = damaging) and class
- `revel` — REVEL ensemble missense score
- `sift`, `polyphen` — SIFT and PolyPhen2-HDIV predictions
- `spliceai_ds_max` — max of SpliceAI delta scores (acceptor/donor × gain/loss). May be empty if SpliceAI not loaded for this run.

Clinical and population:
- `clin_sig`, `clin_stars` — ClinVar significance and review-status star count (0–4)
- `gnomad_af`, `gnomad_faf`, `gnomad_popmax_af` — gnomAD allele frequency, filtering allele frequency, and population-max AF

Gene-level (constant within a gene):
- `loeuf`, `pli`, `mis_z`, `syn_z` — gnomAD constraint metrics
- `omim_phenotype` — OMIM phenotype text for the gene, if available

Tiering:
- `tier` — 1 to 5 (see §4.2)
- `tier_reason` — short text explaining why this tier was assigned
- `modifier_candidate` — TRUE / FALSE (see §4.4)
- `modifier_evidence` — short text explaining the modifier flag

Cohort (from Stage 3.5):
- `cohort_ac`, `cohort_an`, `cohort_af` — internal cohort allele count/number/freq
- `n_het`, `n_hom`, `n_carriers` — cohort genotype distribution
- `n_samples_annotated` — distinct samples for whom this variant appeared in the annotated tab (slightly different from n_carriers; useful for tracking annotation coverage)

Secondary worst-across-any-tx (kept for review of isoform-specific hits):
- `max_alphamissense_any_tx`, `min_esm1b_any_tx` — worst score across all transcripts that overlap this variant, not just the representative transcript

### 5.2 cohort.genes.tsv (one row per gene)

Identity and description:
- `gene` — HGNC symbol
- `Gene` — Ensembl gene ID
- `BIOTYPE` — protein_coding, lncRNA, etc.
- `OMIM_phenotype` — OMIM phenotype text
- `gene_description` — concatenated text: SYMBOL [BIOTYPE] Gene_ID | LOEUF=X | OMIM: phenotype
- `LOEUF`, `pLI`, `mis_z`, `syn_z` — gnomAD constraint

Variant-level counts:
- `n_variants` — total tiered variants in this gene
- `n_lof` — pLoF variants
- `n_splice_ge05` — variants with SpliceAI ≥ 0.5
- `n_clinvar_plp` — variants with ClinVar P/LP
- `n_modifier_candidates` — variants with modifier_candidate=TRUE
- `n_tier1`, `n_tier2`, `n_tier3`, `n_tier4`, `n_tier5` — variant counts per tier
- `best_tier` — minimum tier across the gene (1 = strongest signal)

Score extremes (across all variants in the gene):
- `max_alphamissense`, `min_esm1b`, `max_revel`, `max_spliceai`, `min_gnomad_af`

Cohort sample burden (from Stage 5):
- `n_samples_high_impact` — samples with ≥1 HIGH-impact variant in this gene
- `n_samples_tier1` — samples with ≥1 Tier 1 variant
- `n_samples_tier12` — Tier 1 or 2
- `n_samples_qualifying` — Tier 1 or 2 OR any pLoF at Tier ≤3
- `n_samples_with_modifier_candidate` — samples with ≥1 modifier-flag variant

Non-coding counts (from interim noncoding addendum):
- `n_noncoding_total`, `n_noncoding_rare` — total and rare non-coding variants per gene
- `n_intronic`, `n_utr`, `n_regulatory`, `n_noncoding_tx` — per-class breakdown

Rows are sorted by `best_tier` ascending, then `n_tier1` descending, then
`n_modifier_candidates` descending — top of the file is genes with the
strongest Mendelian or modifier signal.

---

## 6. Interpretation guide

### 6.1 Reading a Tier 1 candidate

A row with `tier=1` is the most likely pathogenic variants in the cohort.
Standard follow-up:

1. Check `tier_reason` to understand which rule fired.
2. Check `n_carriers` and `n_samples_with_variant` — recurrent variants
   across the cohort are stronger candidates than singletons.
3. For pLoF: confirm the gene is constrained (`loeuf` < 0.35), confirm
   `nmd_escape` is FALSE, and verify the consequence call against the
   VCF (fastVEP's call should match Ensembl VEP).
4. For missense: check `alphamissense` and `am_class` — Tier 1 missense
   should have AM ≥ 0.85. Cross-check `revel` and `esm1b` as orthogonal
   confirmation (they should agree, but the tier rule does *not* require
   their agreement to avoid double-counting correlated predictors).
5. For ClinVar P/LP overrides: confirm `clin_stars ≥ 1` and review the
   variant in ClinVar directly.

### 6.2 Reading a modifier candidate

A row with `modifier_candidate=TRUE` is a hypothesis for phenotype
modification, not for direct causation. Standard follow-up:

1. Check `modifier_evidence` for the rationale.
2. Check `gene` against the curated modifier list — variants in NF1/NF2
   anchor genes or RAS-MAPK pathway members are the highest-priority hits.
3. Check `cohort_af` and `n_carriers` — modifier signal requires either
   recurrence (multiple carriers) or extreme effect (homozygous).
4. For phenotype-stratified analyses: join `cohort.variants.tsv` against
   `sample_metadata.tsv` on `sample_id` to enrich for variants
   overrepresented in cases vs controls / severe vs mild.

### 6.3 Reading per-gene rollup

`cohort.genes.tsv` is sorted by `best_tier` ascending, so the top rows are
the strongest Mendelian-candidate genes. For modifier-discovery views,
re-sort by `n_modifier_candidates` desc or `n_samples_with_modifier_candidate`
desc:

```bash
# Top modifier-discovery genes
sort -t$'\t' -k<col>,<col>nr cohort.genes.tsv | head -20
```

The `n_samples_*` columns are the key signals — a gene with 1 Tier 1
variant in 50 of 209 samples (recurrent founder hit) tells a different
story than 50 Tier 1 variants each in 1 sample (private-variant burden).

### 6.4 Reading non-coding counts (interim)

The `n_noncoding_*` columns are **counts only** — no functional score is
attached. They surface that a gene has a high non-coding-variant burden so
you can flag it for follow-up under the v2 plan, but they don't say which
of those variants are functional. A gene with 100 rare non-coding variants
and 0 tiered coding variants is currently invisible to the tiering scheme
beyond these counts.

---

## 7. Known gotchas and limitations

### 7.1 Non-coding variants under-characterized

The current pipeline largely punts on non-coding variants — most rare
non-coding variants land in Tier 4 ("we don't know") because there's no
functional score attached. For NF rare-disease work this is a meaningful
gap: deep-intronic NF1 splice variants are a documented disease cause,
and modifier hypotheses likely involve regulatory variants in pathway
genes. See `docs/noncoding_v2_plan.md` for the planned addition path
(SpliceAI restoration → conservation scores → regulatory element overlap).

### 7.2 SpliceAI currently dropped from Stage 1

For throughput reasons. SpliceAI's 40 GB `.osa` added ~5 min per chr22-scale
annotation, making the cohort run infeasible at scale. Planned remedy:
rebuild the `.osa` cleanly (the current size suggests a build error) AND/OR
run SpliceAI as a post-tier pass against the much smaller tiered variant
set. Both paths described in `docs/noncoding_v2_plan.md`.

### 7.3 Pre-filter coupling to tiering

Any new tiering criterion must also become a pre-filter keep-rule, or
variants with the new signal silently get dropped before tiering can see
them. The pre-filter and tiering rules grow together — see the checklist
in `docs/noncoding_v2_plan.md §7`.

### 7.4 AF filter trade-offs

- `--filter-af 1e-4` (very strict) drops the entire modifier band
  (`1e-4 < popmax ≤ 5e-2`) and produces tier output that omits common-band
  modifier candidates. Useful as a memory hedge when a cohort table is too
  large for R, harmful for modifier discovery.
- `--filter-af 0.01` (default) keeps modifier band, captures most signal.
- `--filter-af 0.05` keeps the full modifier band plus a buffer.

The default 0.01 is the right setting for routine cohort runs.

### 7.5 Single-cohort tiering, not joint-called

The pipeline operates on per-sample VCFs as given. If the VCFs are from
single-sample calling (e.g., GATK HaplotypeCaller on each sample
independently), reference-genotype information may be incomplete — variants
present in one sample's VCF will be absent (not "0/0") from others where
they weren't called. This affects `cohort_an` and per-variant AF
calculations slightly. For more rigorous cohort genotypes, joint-call the
cohort upstream and re-derive per-sample VCFs from the joint VCF.

### 7.6 Sample-level QC not included

The pipeline assumes per-sample VCFs are quality-controlled upstream
(GQ/DP filtering, sample-level PCA, contamination checks, etc.). Outlier
samples carrying unusually many Tier 1+2 variants may indicate QC issues
rather than biological signal — review the per-sample burden figure
(`03_per_sample_burden.png`) for outliers.

### 7.7 Reliance on gnomAD popmax annotation

Variants not in gnomAD have no popmax attached, and under the default
strict pre-filter they fail rarity gating unless they have HIGH/MOD/splice/
ClinVar signal. This is appropriate for the dominant rare-disease use case
(if it's truly rare and in a constrained gene, it'll have impact > MODIFIER
and survive), but it could miss truly novel population-private regulatory
variants in non-anchor genes. The v2 conservation and regulatory scoring
mitigate this.

---

## 8. Validation

`scripts/validate_cohort.sh` runs against finished cohort outputs and
asserts the filtering and tiering rules hold. Checks include:

- File presence + non-empty
- Required column presence
- Tier distribution is pyramid-shaped (T4 > T3 > T2 > T1)
- Class-conditional rules hold (every Tier 1 pLoF in constrained gene; every
  Tier 1 missense has AM ≥ 0.85 or ClinVar P/LP; every Tier 5 is common or
  ClinVar B/LB or AM-benign)
- ClinVar overrides applied (P/LP at ≥1 star always Tier 1; B/LB always Tier 5)
- Modifier flag matches its rule (AF in band or sub-Mendelian, AND class-appropriate signal)
- `n_carriers` matches recount from genotype table (on a 50-variant sample)
- Anchor genes have signal (NF1/NF2/SMARCB1/LZTR1 each have rows)
- Non-coding limitation visible (fraction of rare non-coding variants stuck
  in Tier 4 — documented gap, INFO not FAIL)

Hard violations (real bugs) exit non-zero; INFO items (documented
limitations) don't. Run after every cohort run:

```bash
scripts/validate_cohort.sh --out-dir /data/nf1/outputs/batch1
```

---

## 9. Reproducibility

A given cohort run is reproducible from:

1. Pipeline commit hash (`git rev-parse HEAD` from the repo at run time)
2. Manifest TSV (sample IDs + VCF paths)
3. Reference data versions:
   - Ensembl release number for transcript cache + FASTA
   - dbNSFP version
   - SpliceAI release (Ensembl MANE-restricted file date)
   - ClinVar release date
   - gnomAD version
4. Command-line flags (especially `--filter-af`, `--af-max`, `--include-sex-chroms`)

For publication-ready reporting, record all of the above in your methods
section. The validation script's output is also worth archiving alongside
the cohort outputs.

---

## 10. Future work

See `docs/noncoding_v2_plan.md` for the planned next iteration:

- **Tier A**: SpliceAI restoration (post-tier pass, then rebuild `.osa`)
- **Tier B**: Conservation scoring (PhyloP, GERP)
- **Tier C**: Regulatory-element overlap (ENCODE cCREs, Roadmap chromHMM)

Each independently shippable, prioritized by ROI / effort ratio.

Beyond v2:

- Move Stage 2 aggregations from awk to DuckDB for cohort scaling beyond 1000
  samples (awk hashes start hitting memory limits)
- Convert per-sample annotated tabs to Parquet for ~10× storage reduction
- Add CADD as a genome-wide non-coding signal once we have a fast SA build
  pipeline that handles 80+ GB sources
- Wire up the modifier candidate scoring to include direction-of-effect
  estimates from eQTL databases where available

---

## 11. Citations

The pipeline integrates and depends on:

- **fastVEP** — Huang Lab, GitHub:Huang-lab/fastVEP. Rust reimplementation of
  Ensembl VEP.
- **AlphaMissense** — Cheng, J. et al. (2023). Accurate proteome-wide missense
  variant effect prediction with AlphaMissense. *Science* 381, eadg7492.
- **ESM1b** — Brandes, N. et al. (2023). Genome-wide prediction of disease
  variants with a deep protein language model. *Nat Genet* 55, 1512–1522.
- **REVEL** — Ioannidis, N.M. et al. (2016). REVEL: An Ensemble Method for
  Predicting the Pathogenicity of Rare Missense Variants. *Am J Hum Genet*
  99, 877–885. ClinGen SVI calibration: Pejaver, V. et al. (2022). *Am J Hum
  Genet* 109, 2163–2177.
- **SpliceAI** — Jaganathan, K. et al. (2019). Predicting Splicing from
  Primary Sequence with Deep Learning. *Cell* 176, 535–548.e24.
- **dbNSFP** — Liu, X. et al. (2020). dbNSFP v4: a comprehensive database of
  transcript-specific functional predictions and annotations for human
  nonsynonymous and splice-site SNVs. *Genome Med* 12, 103.
- **ClinVar** — Landrum, M.J. et al. (2020). ClinVar: improvements to
  accessing data. *Nucleic Acids Res* 48, D835–D844.
- **gnomAD** — Karczewski, K.J. et al. (2020). The mutational constraint
  spectrum quantified from variation in 141,456 humans. *Nature* 581,
  434–443. v4: Chen, S. et al. (2024). *Nature* 625, 92–100.
- **ACMG/AMP framework** — Richards, S. et al. (2015). Standards and
  guidelines for the interpretation of sequence variants. *Genet Med* 17,
  405–424.
- **LOFTEE** — Karczewski, K.J. et al. (2020). gnomAD constraint metrics.

For NF1 modifier biology and the curated modifier-gene list, see refs in
`R/nf_modifier_genes.txt`.
