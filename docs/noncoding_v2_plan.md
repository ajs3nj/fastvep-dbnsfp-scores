# Non-coding variants: gap and v2 plan

This document describes how the current cohort_pipeline.sh handles non-coding variants
(badly), the resulting gaps in tiered output, and the concrete additions that would let
us prioritize non-coding candidates the way we already prioritize coding ones.

Audience: collaborators running cohort_pipeline.sh on rare-disease WGS cohorts (e.g. the
NF1/NF2 work that motivated this) who care about regulatory, intronic, splice-altering,
and UTR variants in addition to the coding signal we currently score.

---

## 1. What "non-coding" means here

Non-coding consequence classes from fastvep / VEP, in rough order of importance for
rare-disease variant interpretation:

| class | example consequence terms | typical fraction of WGS variants |
|---|---|---|
| canonical splice | splice_acceptor_variant, splice_donor_variant | <0.1% (counted as pLoF, not non-coding here) |
| splice region | splice_region_variant, splice_donor_5th_base_variant, splice_polypyrimidine_tract_variant | ~0.5% |
| 5' UTR | 5_prime_UTR_variant | ~1% |
| 3' UTR | 3_prime_UTR_variant | ~2% |
| intronic | intron_variant | ~35–45% |
| upstream/downstream | upstream_gene_variant, downstream_gene_variant | ~5–10% |
| regulatory | regulatory_region_variant, TF_binding_site_variant, TFBS_ablation, TFBS_amplification | ~1% |
| intergenic | intergenic_variant | ~30–40% |
| non-coding transcript | non_coding_transcript_exon_variant, non_coding_transcript_variant, mature_miRNA_variant | ~5% |
| synonymous | synonymous_variant, start_retained, stop_retained | ~1% |

Roughly: **non-coding is ~95% of variants in a typical WGS sample.** Whatever signal
we can pull from this fraction is where the rare-disease modifier hypothesis lives —
and where most under-explored disease causality currently is.

## 2. Current state: we punt on non-coding

### Annotation phase

- **dbNSFP scores** (AlphaMissense, ESM1b, REVEL) are coding-only. They don't apply to
  any non-coding variant.
- **SpliceAI** *would* apply to splice-region and deep-intronic variants within ±50 nt
  of MANE exons, **but is currently dropped from Stage 1 for throughput reasons** (a 40 GB
  `.osa` that adds ~5 minutes per sample). So non-coding variants have no splice score
  attached.
- **gnomAD constraint** (LOEUF, pLI) is gene-level and attaches to any variant with a
  gene_id, including intronic and regulatory variants — but constraint reflects coding
  intolerance, not regulatory importance.
- **ClinVar** entries exist for some non-coding variants but coverage is sparse outside
  splice sites.
- **No conservation, regulatory-element, or chromatin-state annotation** is added.

### Pre-filter phase

The Stage 3 filter keeps a variant if any of:
1. impact == HIGH or MODERATE (coding hits)
2. SpliceAI ds_max >= 0.2 (currently always FALSE — SpliceAI not loaded)
3. clin_sig non-empty (mostly coding)
4. gnomAD popmax <= 1e-4 (rarity-only path)

A non-coding variant survives pre-filter **only** by being rare in gnomAD. Once it's
through, it has no functional score to drive tiering.

### Tiering phase

The non-coding class rules in `tier_variants.R`:

| tier | rule | currently met? |
|---|---|---|
| 1 | (none) | n/a |
| 2 | rare + SpliceAI >= 0.5 | almost never (SpliceAI off) |
| 3 | rare + SpliceAI 0.2–0.5 | almost never |
| 4 | rare, no signal | **everything else** |
| 5 | common | filter-out |

Modifier-flag rules for non-coding are similar — require SpliceAI 0.2–0.5 or membership
in the `regulatory` class. Neither path lights up without SpliceAI.

### Result

A 209-sample NF1 cohort produces ~7M rare variants after strict filtering. Of those,
roughly 5–6M are non-coding. **All of them land in Tier 4 by default** with no further
discrimination. A deep-intronic NF1 cryptic-splice variant looks identical in the
output to a random intronic SNP 200 kb from any gene.

For an NF1 cohort this is a real omission: deep-intronic NF1 variants are documented
in ~5% of clinically-confirmed cases, and modifier hypotheses almost certainly involve
regulatory variants in pathway genes.

## 3. v2 plan: three tiers of additions

In priority order. Each tier is independently useful; later tiers depend on earlier
ones only loosely.

### Tier A — Restore SpliceAI cleanly (highest ROI)

**Why it matters first.** SpliceAI is the only existing predictor that operates on
non-coding variants without requiring new infrastructure. It catches cryptic splice
in intronic and UTR variants within ±50 nt of an annotated exon — exactly the
deep-intronic NF1 / NF2 variants we're missing.

**Why the current SpliceAI build is broken.** The `.osa` is 40 GB (vs the typical
5–10 GB for a MANE-restricted build). Per-variant lookup is ~5 minutes per chr22-sized
slice, vs <1 second for any other SA source. The root cause is likely (a) v1 fastSA
format instead of v2, or (b) the sort step before sa-build wasn't applied / was
incomplete.

**Concrete steps.**

```bash
# 1. Re-download the canonical Ensembl MANE-restricted SpliceAI SNV file
wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz

# 2. Verify it's sorted (--reheader-required check)
bcftools view spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz | grep -v '^#' | \
  awk 'BEGIN{lc=""; lp=0} {if ($1 != lc) {lc=$1; lp=0} if ($2 < lp) {print "UNSORTED at " $1":"$2; exit 1} lp=$2}' && \
  echo "OK sorted"

# 3. Rebuild with explicit v2 (chunked, bloom-filtered) format
fastvep sa-build --source spliceai \
  --input spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz \
  --output $SA/spliceai \
  --assembly GRCh38
# expected output: $SA/spliceai.osa (5-10 GB) and $SA/spliceai.osa.idx (a few MB)
du -sh $SA/spliceai.osa

# 4. Re-time on chr22 — target is ~30 seconds, not 5 minutes
time fastvep annotate -i chr22_test.vcf --sa-dir $SA --output-format vcf -o /tmp/out.vcf
```

**Alternative path (post-tier SpliceAI pass) if rebuild doesn't help.** Even at slow
throughput, SpliceAI on the ~7M tiered variants (not the full per-sample 1B) is a
single fastvep invocation taking ~5 hours. Acceptable as an overnight post-tier
addition. Recipe:

```bash
# After cohort.variants.tsv exists, build a small VCF of just the tiered variant keys
tail -n+2 cohort.variants.tsv | awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) c[$i]=i; print "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; next }
  { print $c["chrom"]"\t"$c["pos"]"\t.\t"$c["ref"]"\t"$c["alt"]"\t.\t.\t." }
' > tiered.vcf

# Annotate just these with the full SA dir
fastvep annotate -i tiered.vcf --sa-dir $SA_with_spliceai \
  --output-format vcf -o tiered.spliceai.vcf

# Join spliceai_ds_max from that VCF back into cohort.variants.tsv
# (small Python/awk script extracts the SpliceAI INFO field per variant key
# and adds a spliceai_ds_max column).
```

**R tiering changes**: none, the existing non-coding class rules already use
spliceai_ds_max; they just stop being no-ops once it's populated.

**Expected impact**: ~5–15% of currently-Tier-4 non-coding variants get promoted to
Tier 2 or 3 based on cryptic-splice signal. Deep-intronic NF1 splice variants surface
as candidates.

### Tier B — Conservation scoring (PhyloP, GERP)

**Why it matters.** Conservation tells us "is this position evolutionarily preserved?"
For non-coding variants without splice signal, conservation is the next-best functional
prior — conserved non-coding positions are far more likely to be functional than
non-conserved ones. PhyloP-100way and GERP++ are the standard scores.

**Concrete steps.**

```bash
# PhyloP-100way (wigFix format)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
# (or the .wigFix variant; fastvep sa-build --source phylop accepts wigFix)
fastvep sa-build --source phylop \
  --input hg38.phyloP100way.wigFix \
  --output $SA/phylop \
  --assembly GRCh38
# expected: $SA/phylop.osa (~3-5 GB)

# GERP++ (TSV format)
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/All_hg38_RS.bw  
# or download via the standard GERP distribution; convert to TSV if needed
fastvep sa-build --source gerp \
  --input gerp_scores.tsv.gz \
  --output $SA/gerp \
  --assembly GRCh38
# expected: $SA/gerp.osa (~3-5 GB)
```

**Output columns to add** to the wide tab (csq_to_wide_tab.py): `phylop`, `gerp` (just
the numeric scores; one per variant position).

**R tiering changes** in `tier_variants.R`:

```r
# Add to CFG:
phylop_conserved   = 2.0,   # rough functional threshold
gerp_conserved     = 4.0,

# In add_evidence(), add:
v[, conserved := T0(num_col(v, COLS$phylop) >= CFG$phylop_conserved) |
                 T0(num_col(v, COLS$gerp)   >= CFG$gerp_conserved)]

# Non-coding tier rule additions:
v[is_nc & rare & splice_flag & !splice_likely,
  `:=`(tier = 3L, tier_reason = "non-coding + SpliceAI 0.2-0.5")]
v[is_nc & rare & !splice_flag & conserved,                       # NEW
  `:=`(tier = 3L, tier_reason = "non-coding + conservation (PhyloP>=2 or GERP>=4)")]
```

**Expected impact**: another ~10–20% of Tier-4 non-coding variants get promoted to
Tier 3 based on conservation. Reduces the "we don't know" bucket meaningfully.

### Tier C — Regulatory-element overlap (ENCODE cCREs / Roadmap chromHMM)

**Why it matters.** A variant overlapping an enhancer, promoter, or CTCF binding site
in a disease-relevant tissue is dramatically more likely to be functional than a random
non-coding variant. ENCODE candidate cis-Regulatory Elements (cCREs) and Roadmap
Epigenomics chromHMM tracks provide this.

**Concrete steps.**

```bash
# ENCODE cCREs (single BED, ~1M elements)
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb
bigBedToBed encodeCcreCombined.bb encodeCcreCombined.bed

# Roadmap chromHMM for specific tissues (per-tissue 15-state BEDs)
# For NF: neural crest, schwann cell, brain
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E082_15_coreMarks_segments.bed
# (E082 = fetal brain; pick relevant epigenomes for your phenotype)

# Build as custom BED annotations -- fastvep --custom supports BED
# (or build as a fastSA custom interval source)
fastvep sa-build --source custom_bed \
  --input encodeCcreCombined.bed \
  --output $SA/encode_ccre \
  --assembly GRCh38

fastvep sa-build --source custom_bed \
  --input E082_15_coreMarks_segments.bed \
  --output $SA/roadmap_brain \
  --assembly GRCh38
```

**Output columns to add**: `encode_ccre_class` (PLS / pELS / dELS / CA / TF / CTCF or
blank), `roadmap_state` (the 15-state chromHMM label).

**R tiering changes**:

```r
v[, in_active_regulatory := T0(grepl("PLS|pELS|dELS", v[[COLS$encode_ccre_class]])) |
                            T0(grepl("Enh|TssA|TssAFlnk|EnhBiv|TssBiv", v[[COLS$roadmap_state]]))]

# Tier 2 non-coding with strong regulatory + conservation evidence:
v[is_nc & rare & in_active_regulatory & conserved & splice_likely, ...]
v[is_nc & rare & in_active_regulatory & conserved,
  `:=`(tier = 2L, tier_reason = "non-coding + active regulatory + conserved")]
v[is_nc & rare & in_active_regulatory,
  `:=`(tier = 3L, tier_reason = "non-coding + active regulatory element")]
```

**Expected impact**: another ~5–10% of Tier-4 non-coding variants get reclassified.
For NF-specific work, Schwann-cell / neural-crest regulatory tracks would let us
specifically prioritize variants in NF1-relevant regulatory elements.

## 4. Tiering principle for the v2 non-coding rules

The class-conditional design from the existing tiering scheme extends naturally. A
non-coding variant's tier reflects the *strength of the strongest applicable
functional signal*:

| tier | non-coding criterion |
|---|---|
| 1 | rare + SpliceAI >= 0.8 + in pathway/anchor gene (matches existing pLoF Tier 1 standard) |
| 2 | rare + SpliceAI >= 0.5; OR rare + (active regulatory + conserved); OR rare + (motif disruption in TF binding site of anchor gene) |
| 3 | rare + SpliceAI 0.2–0.5; OR rare + conserved (PhyloP>=2 or GERP>=4); OR rare + active regulatory element |
| 4 | rare, no signal — "we don't know" |
| 5 | common, OR benign by all class-appropriate predictors |

The modifier flag for non-coding gains parallel paths:
- modifier_candidate = TRUE if non-coding + AF in modifier band + (any of: SpliceAI
  0.2–0.5, conservation, regulatory element, in anchor/pathway gene).

## 5. Implementation order

1. **Tier A (SpliceAI rebuild)** — single biggest improvement; days of effort.
   Add a post-tier SpliceAI join step now while the rebuild is pending.
2. **Tier B (PhyloP + GERP)** — broad gain across all non-coding; downloads + builds + R
   tiering rule additions; days of effort.
3. **Tier C (ENCODE + Roadmap regulatory)** — biggest interpretability gain for
   pathway/disease-specific work; bigger lift because custom BED integration needs
   fastvep `--custom` testing.

Each tier is independently shippable. Order is by ROI / effort ratio, not strict
dependency.

## 6. Interim noncoding reporting (already in place)

While the v2 work is pending, `cohort.genes.tsv` has been extended with simple
non-coding count columns:

- `n_noncoding_total` — all non-coding variants per gene
- `n_noncoding_rare` — non-coding variants with gnomAD popmax <= 1e-4
- `n_intronic`, `n_utr`, `n_regulatory`, `n_noncoding_tx` — per-category breakdowns

These are counts only, not prioritization — no functional score is attached to any of
these variants in the current output. They surface gene-level non-coding burden so
reports can include statements like "NF1 has N rare non-coding variants across 209
samples, M of them in 5' UTR" while we wait for the v2 functional scoring.

## 7. Pre-filter must grow alongside tiering — AF gate dependency

The Stage 3 pre-filter and the `tier_variants.R` tiering rules share a structural
constraint: **any variant the tiering might want to score has to survive the
pre-filter first**. The current pre-filter rule is:

```
keep if: HIGH/MOD impact   OR   SpliceAI >= 0.2   OR   ClinVar non-empty   OR   popmax <= filter_af
```

If we add a new tiering criterion (e.g., conservation, regulatory overlap) without
adding a matching keep-rule to the pre-filter, **variants with the new signal but no
existing signal AND no rare popmax get dropped silently**. Two concrete failure modes:

1. **Conserved non-coding with missing popmax.** A non-coding variant in a highly
   conserved region (PhyloP > 4) with no gnomAD popmax — exactly the profile of a
   functionally important rare variant in a low-complexity / repetitive region where
   gnomAD coverage is incomplete — would survive Tier B's tiering rule if it reached
   R, but is dropped by Stage 3 today.

2. **Regulatory non-coding with missing popmax.** A variant overlapping an ENCODE
   enhancer of an anchor gene, missing from gnomAD: dropped by Stage 3 even after
   Tier C ships.

**Required pre-filter additions** alongside each tier:

| Tier | new pre-filter keep-rule |
|---|---|
| A (SpliceAI restored) | already in place — `spliceai_ds_max >= 0.2` |
| B (conservation) | `phylop >= 2.0 OR gerp >= 4.0` |
| C (regulatory) | `encode_ccre_class in (PLS, pELS, dELS) OR roadmap_state matches /Enh|TssA/` |

```awk
# Stage 3 awk additions when v2 ships:
else if (col["phylop"] && $col["phylop"] != "" && $col["phylop"]+0 >= 2.0) keep = 1
else if (col["gerp"]   && $col["gerp"]   != "" && $col["gerp"]+0   >= 4.0) keep = 1
else if (col["encode_ccre_class"] && $col["encode_ccre_class"] ~ /PLS|pELS|dELS/) keep = 1
else if (col["roadmap_state"] && $col["roadmap_state"] ~ /Enh|TssA/) keep = 1
```

### Modifier-band variants and the AF cutoff

The `--filter-af` default of `0.01` keeps variants up to 1% popmax, which captures
most of the modifier band (`1e-4 < popmax <= 5e-2`). Some emergency runs (like the
NF1 recovery in late 2025) used `--filter-af 1e-4` to shrink the filtered cohort
table when R OOM'd on the full table — that drops the modifier band entirely and
prevents modifier candidates with no functional signal from reaching tier_variants.R.

**Going forward, use `--filter-af 0.05`** (or leave at `0.01`) so the modifier band
survives Stage 3. The pre-aggregation in Stage 3.5 keeps R's input small enough
without needing the strict rarity cutoff. Strict cutoffs at the pre-filter are a
last-resort memory hedge, not a default — they conflict with the modifier track.

### Future-proofing: when adding any new tier rule

The checklist for adding a new signal to the tiering scheme:

1. Add the column to `csq_to_wide_tab.py` and document it in COLS contract.
2. Add the corresponding keep-rule to Stage 3's awk so variants with this signal
   survive pre-filter even without other evidence.
3. Add the tiering rule to `tier_variants.R`'s `assign_tier_research()`.
4. Add the modifier-flag rule to `add_modifier()` if the signal is relevant in the
   modifier AF band.
5. Update docs/tiering.md and docs/noncoding_v2_plan.md with the new thresholds.

Skipping step 2 produces the same silent-drop bug — the rule is in tiering but
nothing ever satisfies it because the relevant variants don't survive Stage 3.

## 8. Acceptance criteria for v2

A v2 implementation is complete when, on the test cohort:

- The fraction of rare non-coding variants in Tier 4 drops from ~100% to <50%.
- Known deep-intronic NF1 cryptic-splice variants (from published case reports) land
  in Tier 1 or 2 if their SpliceAI score warrants it.
- The per-gene rollup includes non-coding signal alongside coding-tier counts:
  n_tier1_noncoding, n_tier2_noncoding, etc., so cohort-scale enrichment in regulatory
  regions of anchor / pathway genes is visible.
- Modifier flag captures known NF modifier loci even when their variants are purely
  non-coding (validating the "modifier" hypothesis the existing tier scheme was
  designed to support).
