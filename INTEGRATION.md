# Wiring the dbNSFP scores into your fastVEP fork

This is the thin adapter layer. The hard, version-sensitive logic already lives in `src/dbnsfp_scores.rs`
(tested in isolation). Below is where it plugs in. Four touch points, matching the crates described in the
fastVEP manuscript (`fastvep-sa`, `fastvep-core`, `fastvep-io`).

## Touch point 1 — drop the module into `fastvep-sa`

Copy `src/dbnsfp_scores.rs` to `crates/fastvep-sa/src/dbnsfp_scores.rs` and add `mod dbnsfp_scores;` (or
`pub mod`) to that crate's `lib.rs`. It's `std`-only, so no new dependencies.

## Touch point 2 — extend the dbNSFP source parser

Today your `dbnsfp` parser extracts SIFT/PolyPhen. In its per-line handler, after you've split the row into
`fields: &[&str]` and built a `HeaderIndex` once from the header:

```rust
use crate::dbnsfp_scores::{HeaderIndex, parse_row, score_field, quantize};

// once, when the dbNSFP file header is read:
let hidx = HeaderIndex::from_header(header_line);

// per data row:
let vs = parse_row(&fields, &hidx);

// Per-variant scores -> encode for the fastSA value array:
if let Some(r) = vs.revel    { /* push quantize(score_field("revel").unwrap(), r) */ }
if let Some(c) = vs.cadd_phred { /* push quantize(score_field("cadd_phred").unwrap(), c) */ }
if let Some(c) = vs.cadd_raw   { /* push quantize(score_field("cadd_raw").unwrap(), c) */ } // SIGNED -> zigzag

// Per-transcript scores: choose the transcript that matches the consequence's transcript at annotation
// time. Two options:
//   (a) Simplest, matches how positional .osa2 works: collapse to the MANE/canonical transcript's value
//       (or the max |score|) and store one value per (chr,pos,ref,alt).
//   (b) Most faithful: store the per-transcript vector and resolve against the transcript_id during
//       output. This needs a richer SA value layout; see touch point 3.
```

**Design decision you need to make:** AlphaMissense and ESM1b are *per-transcript* (isoform-specific — ~2M
variants are damaging only in specific isoforms). The fastSA `.osa2` model keys on `(chr,pos,ref,alt)`, i.e.
one value per allele, not per transcript. So either:

- **(a) collapse** to one representative value per allele (canonical/MANE transcript, or max-magnitude). Loses
  isoform specificity but fits the existing format with zero schema change. Fine for a first pass.
- **(b) keep per-transcript** by extending the SA value record to hold a small transcript→score map and
  joining on `transcript_id` when emitting output. Faithful to the biology, more code.

For NF work where isoform choice can matter (e.g. *NF1* has alternatively-spliced exon 23a; *NF2/SMARCB1*
isoforms), (b) is the biologically correct target, but I'd ship (a) first, validate, then upgrade. Tell me which
and I'll write that half concretely.

## Touch point 3 — add the fields to `SupplementaryAnnotation`

In `fastvep-core`, the `SupplementaryAnnotation` struct already carries score fields (REVEL, SpliceAI, etc.).
Add:

```rust
pub alphamissense: Option<f32>,
pub alphamissense_class: Option<String>,
pub esm1b: Option<f32>,
pub revel: Option<f32>,          // if not already present
pub cadd_phred: Option<f32>,
pub cadd_raw: Option<f32>,
```

On read-back from the `.osa2` value array, decode with `dbnsfp_scores::dequantize(field, stored)` using the SAME
`ScoreField` so the multiplier/zigzag match. Never hand-divide at the call site.

## Touch point 4 — surface them in tab output (the one that bites)

The tab formatter in `fastvep-io` emits a fixed column set (the manuscript notes a compact tab layout; the full
field set lives in the 48-field CSQ). Your scores won't appear in `--output-format tab` unless you add columns.
Add `AlphaMissense`, `AlphaMissense_pred`, `ESM1b`, `REVEL`, `CADD_PHRED`, `CADD_RAW` to:

- the tab header writer, and
- the per-row writer, pulling from `SupplementaryAnnotation` (print `.` / empty for `None` to match VEP).

JSON output already serializes the full `SupplementaryAnnotation`, so if you'd rather not touch the tab writer,
emit `--output-format json` and flatten in R. But you asked for tab, so extending the writer is the clean answer.

## Touch point 5 — `sa-build --source dbnsfp`

No new `--source` value is needed (you're extending dbNSFP, not adding a source). Just rebuild your dbNSFP
`.osa2` after the parser change so the new value arrays are written:

```bash
fastvep sa-build --source dbnsfp -i dbNSFP4.7a_variant.complete.gz -o sa_databases/dbnsfp --assembly GRCh38
```

Then annotate as usual — scores now come out in the same pass:

```bash
fastvep annotate -i sample.vcf --gff3 Homo_sapiens.GRCh38.115.gff3 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ --output-format tab --hgvs > sample.annotated.tab
```

## Genome-wide CADD (later, optional)

dbNSFP's CADD is coding-only. If you later need intronic/regulatory CADD genome-wide, that's a *new* source
(`--source cadd`) reading `whole_genome_SNVs.tsv.gz` + the indel file, modeled on the REVEL parser (positional,
allele-specific, two values: `CADD_raw` signed, `CADD_phred` unsigned). Build it into `.osa2` once so per-variant
lookups are Bloom-gated/LRU-cached rather than live tabix seeks — that, not the plugin, is what fixes the CADD
slowdown. Say the word and I'll write that parser too.

---

## To finalize this into exact diffs, paste me these 3 files from your fork

1. `crates/fastvep-sa/src/` — the **dbNSFP source parser** (whatever file holds the `--source dbnsfp` logic;
   likely `dbnsfp.rs` or a `sources/` submodule).
2. The **`--source` enum / `sa-build` dispatch** (the `match` that routes `clinvar`, `revel`, `spliceai`,
   `dbnsfp`, …).
3. The **tab output formatter** in `crates/fastvep-io/` (header + row writer) and the
   `SupplementaryAnnotation` struct definition in `crates/fastvep-core/`.

With those I'll turn touch points 2–4 into copy-paste diffs against your actual function signatures.
