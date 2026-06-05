# Variant classification: causation tier + parallel modifier track

Two parallel classifications per variant, computed in one pass over the merged VEP + dbNSFP table:

- **Causation tier** (`tier`, values 1–5): pathogenicity likelihood for the variant acting on its own.
  Tier 1 = most likely pathogenic; Tier 5 = very likely benign. **Class-conditional logic** —
  the *evidence used* depends on the variant class; the *meaning of each tier* does not.
- **Modifier candidate** (`modifier_candidate`, bool + reason): independent flag for variants that
  plausibly modulate phenotype rather than cause disease — typically common or low-frequency variants
  of moderate effect that would never appear in Tier 1.

The same variant can be `tier = 4` (Mendelian-uninformative) **and** `modifier_candidate = TRUE`.
That combination is exactly what a real modifier looks like.

Implemented in `R/tier_variants.R`. Assumes germline data; ACMG-strict mode is a knob (`--mode acmg`).

---

## 1. Causation tier — principled definition

Tiers are a pathogenicity-likelihood gradient. The class-specific table below operationalizes them;
the *meaning* of each tier is fixed:

| Tier | Meaning |
|------|---------|
| **1** | Most likely pathogenic. Strong calibrated evidence appropriate for the class; ClinVar P/LP at ≥1 star is always Tier 1. |
| **2** | Likely pathogenic. One strong calibrated line of class-appropriate evidence, without disqualifiers. |
| **3** | Uncertain, leaning damaging. Suggestive but sub-threshold class-appropriate signal. Manual review warranted. |
| **4** | Uncertain — no information. Rare, but the class-appropriate predictor isn't applicable, is missing, or is neutral. The honest "we don't know" bucket. |
| **5** | Very likely benign. Class-appropriate evidence calls benign, OR fails rarity gate, OR ClinVar B/LB. |

### 1.1 Class-conditional evidence

Variant class is the *most severe consequence* on the representative transcript
(MANE Select → canonical → most-severe SO term). One primary predictor per class. No cross-class
"consensus" — that systematically punishes non-missense variants for lacking predictors that
physically don't apply (the original concern). REVEL, ESM1b, SIFT and PolyPhen ride along on every
row as **review columns**, never as tier inputs.

| Class | Tier 1 | Tier 2 | Tier 3 | Tier 4 | Tier 5 |
|-------|--------|--------|--------|--------|--------|
| **pLoF** (stop_gained, frameshift, splice_acceptor, splice_donor, start_lost, transcript_ablation) | rare + constrained gene (LOEUF<0.35 or pLI≥0.9), not NMD-escape, LOFTEE HC if present | rare LoF in unconstrained gene; or NMD-escape; or LOFTEE LC | — | — | common, or ClinVar B/LB |
| **non-canonical splice** (splice_region, donor_5th_base, polypyrimidine, donor_region) | rare + SpliceAI ds_max ≥ 0.8 | rare + SpliceAI 0.5–0.8 | rare + SpliceAI 0.2–0.5 | rare, no SpliceAI signal | common |
| **missense / protein_altering** | ClinVar P/LP ≥1 star; or rare + AM ≥0.85 in constrained gene | rare + AM likely_pathogenic (>0.564) | rare + AM 0.34–0.564 (ambiguous); or rare missense in constrained gene with AM missing | rare, AM neutral or missing | AM likely_benign (<0.34), or common, or ClinVar B/LB |
| **in-frame indel** | — | rare + LOEUF<0.35; or SpliceAI ≥0.5 | rare in moderately constrained gene; or SpliceAI 0.2–0.5 | rare, no constraint | common |
| **synonymous / UTR / intronic / non-coding tx** | — | rare + SpliceAI ≥0.5 (cryptic splice) | rare + SpliceAI 0.2–0.5 | rare, no SpliceAI signal | common |
| **regulatory** (motif/TFBS, regulatory_region_variant) | — | rare + motif disruption + strong SpliceAI | motif/TFBS overlap + rare | rare, no signal | common |

Wiggle room: a Tier 1 pLoF and a Tier 1 missense use different evidence but represent the same
probability claim. ClinVar P/LP (≥1 star) always overrides up to Tier 1; B/LB always overrides down
to Tier 5. The rarity gate (default gnomAD popmax / FAF ≤ 1e-4 for dominant) failing always sends a
variant to Tier 5.

### 1.2 Why CADD doesn't appear

We're not pulling CADD. CADD itself is genome-wide, but our delivery vehicle (dbNSFP) only carries
CADD for coding variants, and rather than add CADD as a separate full-genome `.osa2` source, we're
treating SpliceAI as the sole signal on non-coding classes. A consequence: rare deep-intronic /
regulatory variants land Tier 4 unless SpliceAI flags them — which is honest given we don't have a
genome-wide non-coding predictor we trust in this pipeline. For *NF1* deep-intronic cryptic-splice
hunting this matters; SpliceAI is your only chance.

### 1.3 Predictor non-independence (the ACMG-strict caveat)

The class-conditional design *also* solves the predictor non-independence problem: by using exactly
one primary predictor per class, no PP3 evidence is double-counted. `--mode acmg` further constrains
to ACMG-calibrated thresholds (AlphaMissense or REVEL as the single missense PP3/BP4 input,
SpliceAI as a separate splicing PP3) and emits P / LP / VUS / LB / B rather than T1–T5.

---

## 2. Modifier candidate track — parallel

Different question: *given that this variant likely isn't causing the disease on its own, could it
modulate phenotype severity / penetrance / age-of-onset?* Real modifiers tend to be:

- **Not Mendelian-rare.** Common (AF 1–5%) or low-frequency. The 1e-4 causation rarity gate would
  filter them out.
- **Moderate-effect.** Missense in the ambiguous-leaning-damaging band, sub-threshold splice nudges,
  LoF in pathway/network genes, regulatory variants.
- **Pathway-related.** Often downstream of or interacting with the disease gene (for NF1: RAS/MAPK
  network — *KRAS*, *NRAS*, *HRAS*, *BRAF*, *RAF1*, *MAP2K1/2*, *MAPK1/3*, *PTPN11*, *SPRED1/2*,
  *RASA1*; and known NF risk loci like *CDKN2A/B* / 9p21, *ATM*, *TP53* for MPNST).

### 2.1 Modifier flag rule (default: genome-wide)

A variant is `modifier_candidate = TRUE` if **all** hold:

1. **AF band**: either gnomAD popmax `1e-4 < AF ≤ 0.05` (the modifier band), or rare but in
   Mendelian tier 3 or 4 (rare + sub-threshold for causation).
2. **Functional signal appropriate for class**: at least one of
   - missense with `0.34 < AlphaMissense ≤ 0.85` (moderate, ambiguous-leaning-damaging — excludes
     both confident-benign and confident-pathogenic which belong to the causation track), OR
   - non-canonical-splice / non-coding with SpliceAI 0.2–0.5 (sub-Mendelian splice nudge), OR
   - any pLoF in any gene (LoF carriers in modifier roles are biologically plausible for
     dose-sensitive pathway components), OR
   - regulatory_region_variant or TFBS / motif overlap.
3. **Not failed-rarity-common** (popmax > 0.05). Above 5% it's a polymorphism, not a candidate.

The `modifier_evidence` column gives a short string ("missense AM=0.72 in modifier AF band"; "pLoF
in pathway gene PTPN11"). No numeric score — that would imply more precision than the rule has.

### 2.2 Gene scope

Default: **genome-wide** — flag any variant meeting the AF + class-signal criteria. This is noisier
but doesn't presume which genes matter.

If you pass `--modifier-genes path/to/genes.txt` (one HGNC symbol or Ensembl ID per line), the flag
is restricted to that list. A starter NF/RAS-MAPK gene list is provided at
`R/nf_modifier_genes.txt` — use it with `--modifier-genes R/nf_modifier_genes.txt`.

### 2.3 Causation × modifier matrix

Both flags are independent. The expected joint distribution:

|                       | modifier_candidate=FALSE | modifier_candidate=TRUE |
|-----------------------|--------------------------|-------------------------|
| **Causation Tier 1**  | classical Mendelian hit  | rare (compound / second-hit edge cases) |
| **Causation Tier 2**  | likely-causal candidate  | uncommon |
| **Causation Tier 3–4**| Mendelian uninformative  | **the modifier sweet spot** |
| **Causation Tier 5**  | common polymorphism      | low-effect common modifier candidate |

---

## 3. Output granularity

### 3.1 Per-variant table (one row per variant)

Identity + VEP: chrom, pos, ref, alt, gene SYMBOL, Gene (Ensembl ID), transcript (representative),
BIOTYPE, MANE_SELECT/CANONICAL flag, consequence, IMPACT, EXON/INTRON, HGVSg, HGVSc, HGVSp,
Existing_variation, DOMAINS, VARIANT_CLASS, PHENO / GENE_PHENO, MOTIF_NAME.

Scores: AlphaMissense + class, ESM1b, REVEL, SpliceAI ds_max + per-acceptor/donor deltas, SIFT,
PolyPhen, gnomAD popmax / FAF / subpopulation AFs, ClinVar CLIN_SIG + review stars, LOFTEE LOF
(HC/LC) if present, gnomAD LOEUF / pLI.

Secondary: `max_alphamissense_any_tx`, `min_esm1b_any_tx` (worst score across *all* transcripts —
catches isoform-only hits at loci like *NF1* exon 23a).

Tiering: `tier`, `tier_reason`, `variant_class`, evidence flags (am_path, splice_likely, is_lof,
constrained, nmd_escape, clinvar_plp, clinvar_blb), `modifier_candidate`, `modifier_evidence`.

Cohort (only if `--genotypes` provided): `cohort_ac`, `cohort_an`, `cohort_af`, `n_het`, `n_hom`,
`n_carriers`.

### 3.2 Per-gene summary (one row per gene)

Identity + description: SYMBOL, Gene, BIOTYPE, OMIM_phenotype (from fastVEP gene annotations),
LOEUF, pLI, mis_z, syn_z, `gene_description` (concatenated human-readable string).

Variant counts: `n_variants`, `n_lof`, `n_splice_ge05`, `n_clinvar_plp`, per-tier counts
`n_tier1..n_tier5`, `best_tier`, `n_modifier_candidates`.

Cohort burden (only if `--genotypes` provided): `n_samples_high_impact` (≥1 IMPACT=HIGH variant
carrier), `n_samples_tier1`, `n_samples_tier12`, `n_samples_with_modifier_candidate`,
`n_samples_qualifying` (configurable; default = pLoF or missense AM>0.564).

Ordering: by `best_tier` asc, then `n_tier1` desc, then `n_modifier_candidates` desc.

**Note on what "gene description" pulls.** SYMBOL + BIOTYPE come from VEP CSQ directly. OMIM
phenotype text and gene constraint come from fastVEP's gene-level annotations (`.oga`), which are
loaded by fastVEP at annotation time. **HGNC long names** (e.g., "neurofibromin 1") are *not* in
standard VEP output — if you want them, the cleanest add is a one-time join against HGNC's free
`hgnc_complete_set.txt`. A stub for that is left in the R script.

---

## 4. Sample-level summaries — input format

Optional `--genotypes path/to/genotypes.tsv` long-format table with columns
`chrom, pos, ref, alt, sample, gt`. `gt` follows VCF convention: `0/0`, `0/1`, `1/1`, `./.`, or
phased `0|1` etc. Anything containing the alt allele increments AC by one per copy.

A `bcftools query` recipe to produce it from your VCF:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input.vcf.gz | \
  awk 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt","sample","gt"} {
    for (i=5; i<=NF; i++) { split($i, a, "="); print $1,$2,$3,$4,a[1],a[2] }
  }' > genotypes.tsv
```

For multi-allelic sites this assumes you've already split with `bcftools norm -m-`. With 400 samples
genome-wide the file is large; consider streaming the VCF directly later if performance bites.

---

## 5. Caveats that matter for NF small-cohort work

- **Class-conditional logic eliminates** the predictor non-independence problem at the tier level,
  but for ACMG write-ups still use `--mode acmg`.
- **NMD escape** uses a position-based heuristic if `EXON`/`INTRON` columns are present (last exon
  or within 50bp of last exon–exon junction → demote Tier 1 to Tier 2). Best-effort; LOFTEE's
  judgement is better when available.
- **Cohort AF is descriptive only.** With ~400 samples, internal AF noise is high; rarity filtering
  for tiering still uses gnomAD/FAF. Per-gene sample-burden counts are exploratory, not powered
  association — do burden testing properly downstream (SKAT, ACAT, etc.).
- **fastVEP maturity** — validate consequence concordance vs Ensembl VEP on a chromosome before
  trusting tiers cohort-wide. IMPACT and consequence drive the LoF tiers.
- **Modifier flag is hypothesis-generating.** The genome-wide default will surface many candidates
  per sample; that's expected. Filter to your candidate-gene list or pathway for analysis.

## Expected output on the bundled example

Running:

```
Rscript R/tier_variants.R \
  --input tests/example_annotated_tx.tsv \
  --genotypes tests/example_genotypes.tsv \
  --out-prefix /tmp/demo
```

should reproduce these per-variant results (logic mirrored in Python and verified):

| variant | gene | class | tier | modifier | reason |
|---|---|---|---|---|---|
| 17:31226649 C>T | NF1 | pLoF | 1 |  | pLoF constrained, not NMD-escape, LOFTEE HC |
| 17:31229119 G>A | NF1 | pLoF | 1 |  | pLoF (canonical splice donor) constrained + SpliceAI 0.97 |
| 17:31259089 G>A | NF1 | missense | 1 |  | AM 0.98 (≥0.85) in constrained gene |
| 22:29687290 C>T | NF2 | noncoding | 2 |  | synonymous + SpliceAI 0.78 (cryptic splice) |
| 22:29683017 A>G | NF2 | missense | 2 |  | AM 0.61 likely_pathogenic |
| 17:31346200 T>delA | NF1 | pLoF | 2 |  | frameshift but LOFTEE LC + NMD-escape (last exon) → demoted |
| 17:31309500 A>G | NF1 | noncoding | 3 | **Y** | intron + SpliceAI 0.42 — Mendelian uncertain, *modifier candidate* |
| 12:25245350 C>A | KRAS | missense | 5 | **Y** | common (AF 0.012), AM 0.72 in modifier band — *modifier candidate* |
| 3:56865792 G>C | ARHGEF3 | missense | 5 | **Y** | common (AF 0.031), AM 0.45 ambiguous — *modifier candidate* |
| 9:21974827 C>T | CDKN2A | regulatory | 5 | **Y** | common (AF 0.018), CTCF motif — *modifier candidate* |
| 22:29664398 T>C | NF2 | missense | 5 |  | AM 0.21 likely_benign |
| 1:45331556 A>G  | PTEN | missense | 5 |  | ClinVar Benign, common |

Per-gene rollup ordered by `best_tier`: NF1 best=1 (n_tier1=3, n_tier2=1, n_modifier=1); NF2 best=2;
KRAS / ARHGEF3 / CDKN2A best=5 but each with n_modifier_candidates=1 — these dominate the modifier
track even though they're absent from the causation top tiers. With the bundled 4-sample genotype
file, NF1 has 4/4 qualifying carriers; KRAS / ARHGEF3 / CDKN2A each have 2–3 modifier carriers.

## References

- Richards et al. 2015, *Genet Med* — ACMG/AMP framework.
- Pejaver et al. 2022, *Am J Hum Genet* — ClinGen SVI in-silico calibration (REVEL thresholds).
- Cheng et al. 2023, *Science* — AlphaMissense (class thresholds 0.34 / 0.564).
- Brandes et al. 2023, *Nat Genet* — ESM1b genome-wide variant effects.
- Jaganathan et al. 2019, *Cell* — SpliceAI (0.2 / 0.5 / 0.8 thresholds).
- Karczewski et al. 2020 / gnomAD v4 — LOEUF, pLI, filtering allele frequency.
- Karczewski et al. 2017 (LOFTEE) — LoF confidence (HC/LC).
