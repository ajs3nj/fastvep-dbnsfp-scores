# fastvep-dbnsfp-scores

Add **AlphaMissense**, **ESM1b**, **REVEL**, and (coding) **CADD** to a [fastVEP](https://github.com/Huang-lab/fastVEP)
annotation run by extending fastVEP's existing `dbNSFP` supplementary-annotation source — the "dbNSFP route".

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

## Why this design

I built the version-sensitive, easy-to-get-wrong logic (column selection, transcript alignment, signed-score
encoding) as a standalone module that compiles and is tested in isolation, so it can't be silently wrong against
guessed APIs. Wiring it into your fork is then a thin, documented adapter. See `INTEGRATION.md`.

## Important caveats

- **CADD via dbNSFP is coding-only.** dbNSFP only carries CADD for nonsynonymous/splice SNVs. If you need
  genome-wide (intronic/regulatory) CADD for your NF work, that needs the full `whole_genome_SNVs.tsv.gz` as its own
  fastSA source — there's a stub note in `INTEGRATION.md`. AlphaMissense and ESM1b are missense-only *by design*, so
  dbNSFP coverage is complete for them.
- **ESM1b and CADD-raw are signed.** They must use the zigzag path (`quantize_signed`). Using the unsigned path will
  wrap negative values into huge positive integers — silently corrupting scores. The module enforces this.
- **dbNSFP column names/order change between versions.** The parser keys off the header line by name, never by fixed
  index. Always confirm the column names against the `.readme` shipped with your dbNSFP build.
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
