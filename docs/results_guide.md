# Results guide

This document explains what's in each output file and how the tier assignments are made.
It's intended as the first reference for anyone reading the cohort outputs — open this
alongside `cohort.variants.tsv` and you should be able to interpret every column.

For the comprehensive methods (data sources, rationale, validation, gotchas), see
`pipeline_methods.md`. For the original tiering design spec, see `tiering.md`. For pipeline
operation, see `INSTALL.md`.

---

## 1. Tiering — current rules

Every variant in `cohort.variants.tsv` gets exactly one `tier` (1–5) and one `tier_reason`
string explaining which rule fired.

### 1.1 Tier definitions

| Tier | Meaning | Typical evidence |
|------|---------|------------------|
| **1** | Most likely pathogenic. Strong calibrated evidence appropriate for the class. | pLoF in constrained gene, AM ≥ 0.85 in constrained gene, AM + ESM1b orthogonal agreement in constrained gene, NF-anchor pLoF / AM-strong missense, ClinVar P/LP |
| **2** | Likely pathogenic. One strong line of class-appropriate evidence. | pLoF in non-constrained gene or NMD-escape, AM likely_pathogenic without orthogonal confirmation, ESM1b damaging without AM, non-canonical splice with SpliceAI 0.5–0.8 |
| **3** | Uncertain, leaning damaging. Sub-threshold class-appropriate signal. | AM ambiguous (0.34–0.564), in-frame in moderately constrained gene, non-canonical splice with SpliceAI 0.2–0.5 |
| **4** | Uncertain — no information. Rare, but no class-appropriate predictor is available or neutral. | Rare missense outside AM coverage with ESM1b not damaging, rare non-coding outside splice regions, rare regulatory without motif evidence |
| **5** | Very likely benign. Predictor calls benign OR fails rarity gate OR ClinVar B/LB. | AM likely_benign, common (popmax > 1e-4), ClinVar B/LB |

### 1.2 Rule table (ordered as they apply)

The R script applies rules in the order below. Each rule's `WHEN` is a logical filter; if a
row matches, its `tier` and `tier_reason` are overwritten. Later rules win when their
conditions overlap with earlier ones — the late-override block at the bottom always wins.

| # | Variant class | WHEN | THEN tier | tier_reason |
|---|---------------|------|-----------|-------------|
| 1 | pLoF | rare ∧ constrained ∧ ¬NMD-escape ∧ (LOFTEE-HC ∨ LOFTEE absent) | 1 | `pLoF in constrained gene; not NMD-escape` |
| 2 | pLoF | rare ∧ (¬constrained ∨ NMD-escape ∨ LOFTEE-LC) | 2 | `pLoF in non-constrained gene / NMD-escape` |
| 3 | splice_noncanonical | rare ∧ SpliceAI ≥ 0.8 | 1 | `non-canonical splice + SpliceAI>=0.8` |
| 4 | splice_noncanonical | rare ∧ 0.5 ≤ SpliceAI < 0.8 | 2 | `non-canonical splice + SpliceAI 0.5-0.8` |
| 5 | splice_noncanonical | rare ∧ 0.2 ≤ SpliceAI < 0.5 | 3 | `non-canonical splice + SpliceAI 0.2-0.5` |
| 6 | missense | rare ∧ AM ≥ 0.85 ∧ constrained | 1 | `missense: AM>=0.85 in constrained gene` |
| 7 | missense | rare ∧ AM > 0.564 ∧ ¬(AM ≥ 0.85 ∧ constrained) | 2 | `missense: AM likely_pathogenic` |
| 8 | missense | rare ∧ constrained ∧ AM > 0.564 ∧ ESM1b ≤ −7.5 | 1 | `missense: AM likely_pathogenic + ESM1b damaging (orthogonal agreement) in constrained gene` |
| 9 | missense | rare ∧ constrained ∧ AM missing ∧ ESM1b ≤ −7.5 | 2 | `missense: ESM1b damaging in constrained gene (AM missing)` |
| 10 | missense | rare ∧ 0.34 < AM ≤ 0.564 ∧ AM ≤ 0.564 | 3 | `missense: AM ambiguous (0.34-0.564)` |
| 11 | missense | rare ∧ AM missing ∧ constrained ∧ ESM1b not damaging | 3 | `rare missense in constrained gene; AM missing, ESM1b not damaging` |
| 12 | inframe | rare ∧ (LOEUF < 0.35 ∨ SpliceAI ≥ 0.5) | 2 | `in-frame in constrained gene or SpliceAI>=0.5` |
| 13 | inframe | rare ∧ (constrained ∨ SpliceAI ≥ 0.2) ∧ tier > 2 | 3 | `in-frame in moderately constrained or SpliceAI 0.2-0.5` |
| 14 | noncoding | rare ∧ SpliceAI ≥ 0.5 | 2 | `non-coding + SpliceAI>=0.5 (cryptic splice)` |
| 15 | noncoding | rare ∧ 0.2 ≤ SpliceAI < 0.5 | 3 | `non-coding + SpliceAI 0.2-0.5` |
| 16 | regulatory | rare ∧ motif present ∧ SpliceAI ≥ 0.5 | 2 | `regulatory: motif + SpliceAI>=0.5` |
| 17 | regulatory | rare ∧ (motif ∨ SpliceAI ≥ 0.2) | 3 | `regulatory: motif or SpliceAI 0.2-0.5` |
| — | **late-override block (these always win):** | | | |
| 18 | missense | AM likely_benign | 5 | `missense: AM likely_benign` |
| 19 | any | ¬rare (popmax > 1e-4) | 5 | `common (failed rarity gate)` |
| 20 | any | ClinVar B/LB | 5 | `ClinVar B/LB` |
| 21 | pLoF | rare ∧ gene ∈ {NF1, NF2, SMARCB1, LZTR1} ∧ ¬NMD-escape | 1 | `pLoF in NF-spectrum tumor suppressor (NF1/NF2/SMARCB1/LZTR1); not NMD-escape` |
| 22 | missense | rare ∧ gene ∈ {NF1, NF2, SMARCB1, LZTR1} ∧ AM ≥ 0.85 | 1 | `missense: AM>=0.85 in NF-spectrum tumor suppressor (NF1/NF2/SMARCB1/LZTR1)` |
| 23 | any | ClinVar Pathogenic or Likely_pathogenic | 1 | `ClinVar P/LP (star filter disabled; see methods §7.10)` |

**Note on rule 23:** The star filter (require ClinVar review status ≥ 1 star) is disabled for
this cohort due to an upstream parser bug. The bucket includes 0-star "no-assertion-criteria"
submissions and should be down-weighted when clinical-grade confidence is required. See
`pipeline_methods.md` §7.10.

### 1.3 Definitions used in the rules

| Term | Definition |
|------|------------|
| **rare** | `popmax_af` is missing OR `popmax_af ≤ 1e-4` |
| **constrained** | `LOEUF < 0.35` OR `pLI ≥ 0.9` (gnomAD v4.1 gene constraint) |
| **NMD-escape** | Frameshift/nonsense in the last exon (parsed from VEP EXON i/N field) |
| **LOFTEE-HC/LC** | Loss-of-function transcript effect estimator confidence call. Currently NOT integrated; treated as absent (all pLoF assumed HC). See `pipeline_methods.md` §7.3. |
| **SpliceAI** | Maximum delta-score across {AG, AL, DG, DL} acceptors/donors. Currently dropped from Stage 1 for throughput; non-coding tiering limited. See §7.2. |
| **AM** | AlphaMissense (DeepMind, 2023). Calibrated thresholds: 0.34 (benign), 0.564 (likely-pathogenic), 0.85 (strong). |
| **AM_class** | AlphaMissense categorical call: `likely_benign` / `ambiguous` / `likely_pathogenic`. |
| **ESM1b** | Protein language model log-likelihood ratio (Brandes et al, 2023). Calibrated damaging threshold: ≤ −7.5. |
| **REVEL** | Ensemble missense predictor. Not currently used in tier promotion (correlated with AM); informational only. |
| **variant_class** | Internal grouping: `pLoF` (HIGH-impact LoF terms), `splice_noncanonical`, `missense`, `inframe`, `noncoding`, `regulatory`, `other`. |
| **NF-spectrum tumor suppressor** | Hardcoded set `{NF1, NF2, SMARCB1, LZTR1}`. See §7.9 of methods for rationale (gnomAD v4 mis-classifies LZTR1 as unconstrained). |

### 1.4 Modifier candidate flag (parallel track)

In addition to `tier`, every variant gets a boolean `modifier_candidate` flag computed
independently. A modifier candidate is a variant that plausibly modulates phenotype severity
rather than causing disease — typically a moderate-effect variant in a moderate-frequency band
(`1e-4 < popmax ≤ 5e-2`) in an NF/RAS-MAPK pathway gene. The same variant can be `tier = 4`
(Mendelian-uninformative) AND `modifier_candidate = TRUE`; that combination is exactly what a
real modifier looks like.

The `modifier_evidence` column carries a short string describing why it was flagged
(e.g., `"missense, AM=0.42, popmax=2.3e-3"`).

---

## 2. `cohort.variants.tsv` column dictionary

One row per unique variant `(chrom, pos, ref, alt)` observed in any sample, deduplicated.

### Variant identity

| Column | Source | Description |
|--------|--------|-------------|
| `chrom` | VCF | Chromosome (Ensembl-style: `1`–`22`, `X`, `Y`, `MT`) |
| `pos` | VCF | 1-based position |
| `ref` | VCF | Reference allele |
| `alt` | VCF | Alternate allele |

### Transcript / gene context

| Column | Source | Description |
|--------|--------|-------------|
| `gene` | VEP `SYMBOL` | HGNC gene symbol of the chosen consequence |
| `gene_id` | VEP `Gene` | Ensembl gene ID |
| `transcript` | VEP `Feature` | Ensembl transcript ID for the chosen consequence |
| `mane_select` | VEP | `TRUE` if this is the MANE Select transcript |
| `canonical` | VEP | `TRUE` if the transcript is flagged canonical |
| `mane_plus_clinical` | VEP | `TRUE` if MANE Plus Clinical |
| `biotype` | VEP `BIOTYPE` | e.g., `protein_coding`, `lncRNA`, `processed_pseudogene` |
| `consequence` | VEP `Consequence` | Comma-separated SO terms (e.g., `missense_variant`, `splice_donor_variant`) |
| `impact` | VEP `IMPACT` | `HIGH` / `MODERATE` / `LOW` / `MODIFIER` |
| `variant_class` | Internal | One of: `pLoF` / `splice_noncanonical` / `missense` / `inframe` / `noncoding` / `regulatory` / `other`. Derived from VEP consequence terms |
| `exon`, `intron` | VEP | Exon or intron number as `i/N` |
| `hgvsg`, `hgvsc`, `hgvsp` | VEP | HGVS nomenclature at genome / cDNA / protein level |
| `existing_variation` | VEP | dbSNP/COSMIC IDs known for this position |
| `domains` | VEP | Protein domain overlap (Pfam etc.) |
| `pheno` | VEP `PHENO` | Variant-level phenotype links from VEP |
| `gene_pheno` | VEP | Gene-level phenotype association flag |
| `motif_name` | VEP `MOTIF_NAME` | TF binding motif name if applicable |
| `nmd_escape` | Derived | `TRUE` if the variant is predicted to escape nonsense-mediated decay (last-exon heuristic from EXON i/N) |
| `loftee_lof` | LOFTEE / VEP | LOFTEE LoF confidence. **Currently always empty** — LOFTEE not integrated in v1. See methods §7.3 |

### Population / clinical annotation

| Column | Source | Description |
|--------|--------|-------------|
| `gnomad_af` | gnomAD v4.1 | Overall gnomAD allele frequency |
| `gnomad_faf` | gnomAD v4.1 | Filtering allele frequency (FAF95) — most conservative; used for rarity if present |
| `gnomad_popmax_af` | gnomAD v4.1 | Maximum AF across ancestry groups |
| `af_used` | Derived | The single AF actually used for the `rare` flag: prefers `faf`, then `popmax`, then `af` |
| `clin_sig` | ClinVar | Clinical significance string (`Pathogenic`, `Likely_pathogenic`, `Benign`, etc.) |
| `clin_stars` | ClinVar | Review status mapped to 0–4 stars. **Silently 0 for the whole cohort** in v1 due to a parser bug; star filter is disabled in the R rule. See methods §7.10 |

### Gene-level constraint and phenotype

| Column | Source | Description |
|--------|--------|-------------|
| `loeuf` | gnomAD v4.1 | LoF Observed/Expected upper bound (smaller = more constrained). Constrained threshold: `< 0.35` |
| `pli` | gnomAD v4.1 | Probability of LoF intolerance (0–1). Constrained threshold: `≥ 0.9` |
| `mis_z` | gnomAD v4.1 | Missense Z-score (positive = constrained against missense) |
| `syn_z` | gnomAD v4.1 | Synonymous Z-score (sanity check; should be ≈ 0 for well-modeled genes) |
| `omim_phenotype` | OMIM via fastvep | OMIM phenotype text associated with the gene, if any |

### Functional predictors

| Column | Source | Description |
|--------|--------|-------------|
| `spliceai_ds_max` | SpliceAI | Max delta-score across acceptors/donors. Currently empty for most variants (SpliceAI dropped from Stage 1) |
| `alphamissense` | dbNSFP v4.5+ → AlphaMissense | 0–1 missense pathogenicity. Calibrated thresholds: 0.34 (benign) / 0.564 (likely_pathogenic) / 0.85 (strong) |
| `am_class` | dbNSFP → AlphaMissense | Categorical call: `likely_benign` / `ambiguous` / `likely_pathogenic` |
| `esm1b` | dbNSFP → ESM1b | Log-likelihood ratio (negative = damaging). Calibrated damaging threshold: `≤ −7.5` |
| `esm1b_class` | dbNSFP → ESM1b | Categorical call: `damaging` / `tolerated` |
| `revel` | dbNSFP → REVEL | 0–1 ensemble missense score (informational; not used in tier promotion to avoid double-counting with AM) |
| `sift` | VEP → SIFT | Categorical: `deleterious` / `tolerated` (informational) |
| `polyphen` | VEP → PolyPhen-2 | Categorical: `probably_damaging` / `possibly_damaging` / `benign` (informational) |

### Cohort-level (added at Stage 3.5)

| Column | Source | Description |
|--------|--------|-------------|
| `n_samples_annotated` | Stage 2 | Number of distinct samples whose VCF contained this variant (irrespective of genotype) |
| `cohort_ac` | Stage 3.5 | Cohort allele count summed across called samples |
| `cohort_an` | Stage 3.5 | Cohort total allele number across called samples |
| `cohort_af` | Stage 3.5 | `cohort_ac / cohort_an` |
| `n_het` | Stage 3.5 | Number of samples heterozygous (0/1) for the alt allele |
| `n_hom` | Stage 3.5 | Number of samples homozygous (1/1) for the alt allele |
| `n_carriers` | Stage 3.5 | Number of samples carrying at least one alt allele (`n_het + n_hom`) |

### Tier + modifier flag (added at Stage 4)

| Column | Source | Description |
|--------|--------|-------------|
| `tier` | tier_variants.R | Integer 1–5 (see §1.1) |
| `tier_reason` | tier_variants.R | Short string naming the specific rule that set the tier (see §1.2). **Use this column for filtering / prioritization** |
| `modifier_candidate` | tier_variants.R | Boolean: `TRUE` if the variant is a candidate phenotype modifier (parallel track to tier — see §1.4) |
| `modifier_evidence` | tier_variants.R | Short string explaining why it was flagged (variant class, AM score, AF) |
| `rare` | tier_variants.R | Boolean: `af_used` is missing OR `≤ 1e-4` |
| `constrained` | tier_variants.R | Boolean: gene is LoF-intolerant by either LOEUF or pLI threshold |

---

## 3. `cohort.genes.tsv` column dictionary

One row per gene observed in the cohort.

### Gene identity + constraint

| Column | Description |
|--------|-------------|
| `gene_key` | HGNC symbol (primary key) |
| `SYMBOL`, `Gene`, `BIOTYPE` | HGNC symbol, Ensembl gene ID, biotype |
| `OMIM_phenotype` | OMIM phenotype text, if any |
| `LOEUF`, `pLI`, `mis_z`, `syn_z` | gnomAD v4.1 gene constraint metrics |
| `gene_description` | Pre-formatted human-readable summary string combining the above — useful for plotting / reports |

### Variant counts by tier

| Column | Description |
|--------|-------------|
| `n_variants` | Total distinct variants observed in this gene |
| `n_tier1`, `n_tier2`, `n_tier3`, `n_tier4`, `n_tier5` | Variant count per tier |
| `n_modifier_candidates` | Distinct variants flagged as modifier candidates |

### Sample-level burden (added at Stage 5)

These count distinct **samples** carrying at least one qualifying variant in the gene.

| Column | Description |
|--------|-------------|
| `n_samples_high_impact` | Samples carrying any VEP `HIGH`-impact variant in this gene |
| `n_samples_tier1` | Samples carrying any Tier 1 variant |
| `n_samples_tier12` | Samples carrying any Tier 1 or Tier 2 variant |
| `n_samples_qualifying` | Samples carrying any "qualifying" variant: Tier 1–2 in any class, OR any pLoF at Tier ≤ 3 (catches NMD-escape / LOFTEE-LC LoF that get demoted) |
| `n_samples_with_modifier_candidate` | Samples carrying any modifier candidate |

---

## 4. `cohort.genotypes.tsv` schema

Long-format per-sample carrier table written by Stage 3.5. Tab-separated, gzipped. Six
columns, no header in the file:

| Position | Column | Description |
|----------|--------|-------------|
| 1 | `chrom` | Chromosome |
| 2 | `pos` | 1-based position |
| 3 | `ref` | Reference allele |
| 4 | `alt` | Alternate allele |
| 5 | `sample` | Sample ID (matches the cohort manifest) |
| 6 | `gt` | Genotype string (`0/1`, `1/0`, `0|1`, `1|0`, `1/1`, `1|1`). Reference-only rows (`0/0`) are excluded |

The file is ~32 GB for a 209-sample WGS cohort. It exists so downstream tools (the figures
script, custom analyses) can recover per-sample carrier information without re-running
variant calling. Read it via `awk` streaming for any aggregation — `fread`-ing it into R will
OOM.

---

## 5. Variant ↔ sample relationship

The three output files are designed around one architectural decision: **variant
annotations are per-variant, sample carriership is per-(variant, sample)**, and the two are
kept in separate files joined on `(chrom, pos, ref, alt)`.

### 5.1 Each variant gets exactly one row in `cohort.variants.tsv`

Regardless of how many samples carry it. The pipeline deduplicates twice:

1. **Stage 2** (cohort merge) dedupes on `(chrom, pos, ref, alt, transcript)`. If three
   samples have the same variant on the same transcript, only one annotation row survives.
   Annotations (consequence, AlphaMissense score, gene constraint, etc.) are properties of
   the variant itself, not the carrier, so duplication would just waste storage.
2. **Stage 4** (`tier_variants.R`) collapses across transcripts to one row per
   `(chrom, pos, ref, alt)`. The chosen transcript is MANE Select where available,
   otherwise canonical.

For a 209-sample cohort the final file has ~624k rows even though the source per-sample
VCFs contain millions of variant–sample observations in aggregate.

### 5.2 What's per-variant vs per-sample vs aggregated

| Information | Where it lives | Per what |
|-------------|----------------|----------|
| Variant annotation (consequence, AM, ESM1b, LOEUF, pLI, ClinVar, ...) | `cohort.variants.tsv` | per variant |
| Tier + tier_reason + modifier flag | `cohort.variants.tsv` | per variant |
| Carrier counts (`n_carriers`, `n_het`, `n_hom`) | `cohort.variants.tsv` | per variant (aggregated across samples) |
| Cohort allele counts (`cohort_ac`, `cohort_an`, `cohort_af`) | `cohort.variants.tsv` | per variant (aggregated) |
| Which specific samples carry which variant + their genotype | `cohort.genotypes.tsv` | per (variant, sample) |
| Per-gene rollups (n_tier1, n_samples_tier1, ...) | `cohort.genes.tsv` | per gene (aggregated) |

So a variant carried by 5 samples (3 het, 2 hom alt) gets:
- **1 row** in `cohort.variants.tsv` with `n_carriers=5, n_het=3, n_hom=2, cohort_ac=7, cohort_an=418, cohort_af≈0.0167`
- **5 rows** in `cohort.genotypes.tsv`, one per carrier with that carrier's sample ID and genotype
- contribution to the gene's `n_tier1` / `n_samples_tier1` / etc. counts in `cohort.genes.tsv`

### 5.3 The `sample` column in `cohort.variants.tsv` is NOT meaningful

This column is inherited from the per-sample annotation tables produced by
`csq_to_wide_tab.py`. After the Stage 2 dedupe it holds whichever single sample's row
happened to win the deduplication (first occurrence) — it does **not** identify the
carrier(s) of a multi-sample variant. Use `n_carriers` for "how many samples have this
variant" and join to `cohort.genotypes.tsv` for "which samples specifically".

### 5.4 Joining variants to samples

The join key is `(chrom, pos, ref, alt)` in both files. A few common patterns:

**Get all samples carrying a specific variant:**

```bash
awk -F'\t' '$1=="17" && $2=="31196120" && $3=="C" && $4=="T" {print $5, $6}' \
    cohort.genotypes.tsv
```

**For every Tier 1 NF1 variant, list the carrier samples and their genotypes:**

```bash
# Step 1: emit (key, tier_reason) pairs for Tier 1 NF1 variants
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
  $c["gene"]=="NF1" && $c["tier"]==1 {
    print $c["chrom"]"_"$c["pos"]"_"$c["ref"]"_"$c["alt"] "\t" $c["tier_reason"]
  }
' cohort.variants.tsv | sort > tier1_nf1.keys

# Step 2: emit (key, sample, gt) from genotypes
awk -F'\t' '{ print $1"_"$2"_"$3"_"$4 "\t" $5 "\t" $6 }' \
    cohort.genotypes.tsv | sort > genotypes.keys

# Step 3: join
join -t $'\t' tier1_nf1.keys genotypes.keys
```

The same pattern (variant key + sort + join) works for any tier subset.

**Get the sample list for every Tier 1 variant in a single command** (no sorting needed if
the genotypes file is small enough to hash in memory — for 209 samples this is fine on
modest RAM if you first filter the variants table):

```r
library(data.table)
v  <- fread("cohort.variants.tsv")[tier == 1L]
gt <- fread("cohort.genotypes.tsv",
            col.names = c("chrom","pos","ref","alt","sample","gt"))
v_keys <- v[, .(chrom = as.character(chrom), pos = as.integer(pos),
                ref, alt, gene, tier_reason)]
gt[, chrom := as.character(chrom)]
joined <- merge(v_keys, gt, by = c("chrom","pos","ref","alt"))
# joined now has one row per (Tier 1 variant, carrier sample) pair
```

For larger filters (e.g., all Tier 1+2 ≈ 130k variants), the awk-stream pattern is
faster and avoids loading the 32 GB genotype file into R.

---

## 6. Reading recipes

### "Show me the strongest NF1 candidates"

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
            $c["gene"]=="NF1" && $c["tier"]==1' \
    cohort.variants.tsv | column -t -s $'\t' | less -S
```

### "How many samples carry a Tier 1 variant in each NF anchor gene?"

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
            $c["gene_key"] ~ /^(NF1|NF2|SMARCB1|LZTR1)$/
            { print $c["gene_key"], $c["n_samples_tier1"] }' \
    cohort.genes.tsv
```

### "List Tier 1 missense by tier_reason bucket"

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
            $c["tier"]==1 && $c["variant_class"]=="missense"
            { n[$c["tier_reason"]]++ }
            END { for (r in n) print n[r], r }' \
    cohort.variants.tsv | sort -rn
```

### "Filter to high-confidence Tier 1 only (exclude ClinVar-relaxed bucket)"

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
            $c["tier"]==1 && $c["tier_reason"] !~ /ClinVar P\/LP/' \
    cohort.variants.tsv > tier1_strict.tsv
```

This drops the ClinVar relaxed bucket (which currently includes 0-star submissions) and
keeps only variants reaching Tier 1 via constraint-based or orthogonal-predictor evidence.

### "Find compound-het candidates: samples with ≥2 distinct Tier 1+2 variants in one gene"

The pre-computed answer is in figure `08_sample_gene_hotspots.png`. To recompute manually:

```bash
# Build (chrom,pos,ref,alt) -> gene table for Tier 1+2 variants
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
            $c["tier"] <= 2 { print $c["chrom"]"_"$c["pos"]"_"$c["ref"]"_"$c["alt"], $c["gene"] }' \
    cohort.variants.tsv > t12_variants_by_gene.txt

# Stream genotypes, count (sample, gene) pairs with ≥2 distinct T1+T2 variants
awk -F'\t' '
    NR==FNR { vg[$1] = $2; next }
    $6 ~ /^(0\/1|1\/0|0\|1|1\|0|1\/1|1\|1)$/ {
      k = $1"_"$2"_"$3"_"$4
      if (k in vg) cnt[$5"|"vg[k]]++
    }
    END { for (sg in cnt) if (cnt[sg] >= 2) print cnt[sg], sg }' \
    t12_variants_by_gene.txt cohort.genotypes.tsv | sort -rn
```

---

## 7. Where to read more

| Question | See |
|----------|-----|
| Why was rule X chosen? | `pipeline_methods.md` §4 (tiering rules) and §5 (modifier flag) |
| What are the known gaps in v1? | `pipeline_methods.md` §7 (Known gotchas and limitations) |
| How do I re-run the pipeline? | `INSTALL.md` |
| What's planned for v2? | `noncoding_v2_plan.md` |
| What's the original tiering design? | `tiering.md` |
