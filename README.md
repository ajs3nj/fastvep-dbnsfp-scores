# fastvep-dbnsfp-scores

Add **AlphaMissense**, **ESM1b**, and **REVEL** to a [fastVEP](https://github.com/Huang-lab/fastVEP)
annotation run by extending fastVEP's existing `dbNSFP` supplementary-annotation source — the "dbNSFP route".
Plus an R variant-classification layer that produces per-variant tiers, a parallel modifier-candidate flag,
cohort sample summaries, and a per-gene rollup with descriptions.

The idea: dbNSFP v4.5+ already bundles all four scores, mapped to GRCh38 coordinates and aligned per Ensembl
transcript. fastVEP already has a `--source dbnsfp` parser, but today it only extracts SIFT/PolyPhen. This repo
gives you the column logic and encoding to pull the four extra scores out of the *same* dbNSFP build, so they come
out in the same VEP pass instead of being joined afterward.

## What's here

| Path | Purpose |
|------|---------|
| `docs/dbnsfp_score_columns.md` | The exact dbNSFP columns to read, their multi-value/transcript-alignment semantics, signedness, and recommended quantization multipliers. **Read this first.** |
| `src/dbnsfp_scores.rs` | Self-contained, `std`-only Rust module that does the hard part correctly: header-indexed column lookup (version-safe), per-transcript `;`-split alignment, missing-value handling, and `value*multiplier -> u32` quantization with **zigzag for signed scores** (ESM1b, CADD raw). Has unit tests. Compiles on its own — no dependency on fastVEP internals. |
| `tools/extract_dbnsfp_scores.py` | Standalone extractor (pure Python 3, stdlib only). Reads a bgzipped dbNSFP TSV and emits a tidy per-(variant, transcript) score table. Use it (a) to get scores today without touching Rust, and (b) as a **validation oracle** to check the Rust path against. |
| `tests/sample_dbnsfp.tsv` | Tiny synthetic dbNSFP fixture (header + a few rows) used by both the Rust tests and the Python tool. |
| `INTEGRATION.md` | Step-by-step: where to call `dbnsfp_scores.rs` from your fork's `dbnsfp` parser, what to add to the `sa-build --source` dispatch, the `SupplementaryAnnotation` fields, and the tab formatter. Lists the **3 source files I need from your fork** to turn this into exact diffs. |
| `docs/tiering.md` | Revised variant **tiering scheme**: per-variant tiers from VEP consequence/IMPACT + gnomAD AF/constraint + ClinVar + SpliceAI + AlphaMissense/ESM1b/REVEL/CADD. Germline research-prioritization default, ACMG-strict as a knob, calibrated thresholds, the predictor non-independence caveat, and NF-specific notes. |
| `R/tier_variants.R` | data.table implementation: class-conditional causation tiering, **parallel modifier flag**, optional cohort sample summaries from `--genotypes`, per-gene rollup with description + sample burden. Tunable thresholds at the top. Logic verified on the example; run on EC2 to confirm against `docs/tiering.md`. |
| `R/nf_modifier_genes.txt` | Starter NF/RAS-MAPK candidate-gene list for `--modifier-genes` restriction. |
| `tests/example_annotated_tx.tsv` | Synthetic merged input across multiple variant classes (NF1/NF2/PTEN/KRAS/ARHGEF3/CDKN2A) for the tiering script. |
| `tests/example_genotypes.tsv` | Synthetic 4-sample genotype file demonstrating cohort summaries. |
| `scripts/cohort_pipeline.sh` | **End-to-end cohort pipeline**: takes a manifest TSV of `(sample_id, vcf_path)` rows and produces a single cohort `variants.tsv` + `genes.tsv`. Per-sample annotation runs in parallel; merge step deduplicates by `(chrom,pos,ref,alt,transcript)` and tracks `n_samples_annotated` per variant; pre-filter drops common/no-signal rows; one tier_variants.R run on the cohort. Idempotent with `--skip-annotate` / `--skip-merge`. |
| `scripts/cohort.manifest.example.tsv` | Example manifest format. |

## Cohort workflow

For an N-sample cohort run (the main use case), use the pipeline script rather than running tier_variants.R
per sample. It collapses variants across samples *before* tiering — pathogenicity is a property of the
variant, not the sample carrying it, so tiering 400 times wastes effort.

```bash
# Create a manifest TSV (header required):
#   sample_id<TAB>vcf_path
#   S001<TAB>/data/vcfs/S001.vcf.gz
#   ...
# See scripts/cohort.manifest.example.tsv

scripts/cohort_pipeline.sh \
  --manifest    my.manifest.tsv \
  --out-dir     /path/to/results \
  --data-dir    /path/to/data \    # holds GFF3 + FASTA + transcript cache
  --sa-dir      /path/to/sa \      # holds built .osa/.osi/.oga files
  --threads     16 \
  --modifier-genes R/nf_modifier_genes.txt
```

Outputs in `<out-dir>/`: `cohort.variants.tsv` (one row per unique variant with tier + modifier flag +
`n_samples_annotated` + cohort AC/AN/AF/n_het/n_hom/n_carriers) and `cohort.genes.tsv` (per-gene
rollup with sample-burden columns: n_samples_high_impact / tier1 / tier12 / qualifying /
with_modifier_candidate).

The pipeline tracks how many samples each variant appears in via two columns:
- `n_samples_annotated` (added during merge): number of distinct samples whose VCF contained the variant
- `n_carriers` (added by tier_variants.R from `--genotypes`): number of samples actually carrying the alt allele

For per-sample called VCFs these should agree; for joint-called multi-sample VCFs `n_samples_annotated` will
equal the cohort size and `n_carriers` is the real signal.

## Why this design

I built the version-sensitive, easy-to-get-wrong logic (column selection, transcript alignment, signed-score
encoding) as a standalone module that compiles and is tested in isolation, so it can't be silently wrong against
guessed APIs. Wiring it into your fork is then a thin, documented adapter. See `INTEGRATION.md`.

## Important caveats

- **CADD is intentionally not pulled.** Non-coding classes rely on SpliceAI alone; see `docs/tiering.md §1.2` for
  the design rationale. If you ever want CADD back, the path is a dedicated `--source cadd` Rust parser
  (genome-wide tabix) — out of scope here.
- **AlphaMissense and ESM1b are missense-only by design.** dbNSFP coverage is complete for them; nothing is lost.
- **ESM1b is a signed score (negative LLR = damaging).** The encoding/decoding path handles this correctly; only
  matters if you ever extend the parser to other negative-valued scores.
- **dbNSFP column names/order change between versions.** The parser keys off the header line by name, never by
  fixed index. Always confirm the column names against the `.readme` shipped with your dbNSFP build.
- **fastVEP is young.** Before trusting it across all ~400 samples, validate consequence concordance against Ensembl
  VEP on a chromosome of your own data — the preprint reports 100% concordance on 2,340 transcript-allele pairs,
  which is reassuring but small.

## Getting this onto your EC2 box

This is a plain git repo. Push it to your GitHub, then clone on EC2:

```bash
# from this directory
git remote add origin git@github.com:<you>/fastvep-dbnsfp-scores.git
git push -u origin main
# on EC2
git clone git@github.com:<you>/fastvep-dbnsfp-scores.git
```
