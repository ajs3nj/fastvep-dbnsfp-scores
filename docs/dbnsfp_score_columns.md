# dbNSFP columns for AlphaMissense / ESM1b / REVEL

Target: dbNSFP **v4.5+** (AlphaMissense, ESM1b, EVE were added in v4.5, Nov 2023). v4.7 / v4.9 / 5.x
keep these columns but may add others and shift positions — **always parse by header name, never by
index**, and confirm against the `dbNSFP*.readme.txt` that ships with your download.

**CADD is intentionally not pulled.** See `docs/tiering.md §1.2` for the design choice — non-coding
classes rely on SpliceAI alone.

## Coordinate / key columns (GRCh38)

| Column | Notes |
|--------|-------|
| `#chr` | chromosome (no `chr` prefix in dbNSFP) |
| `pos(1-based)` | GRCh38 position. (There are separate `hg19_chr` / `hg19_pos(1-based)` columns — do **not** use those for a GRCh38 build.) |
| `ref`, `alt` | alleles |
| `aaref`, `aaalt` | reference / alt amino acid (single letter) |
| `Ensembl_transcriptid` | `;`-separated list; defines the per-transcript ordering |
| `Ensembl_proteinid` | `;`-separated, aligned with `Ensembl_transcriptid` |
| `genename` | `;`-separated gene symbol(s) |

## Score columns

| Score | dbNSFP column(s) | Range / type | Per-transcript? | Sign | Quantization |
|-------|------------------|--------------|-----------------|------|--------------|
| **AlphaMissense** | `AlphaMissense_score` | 0.0–1.0 (pathogenicity) | **Yes** (`;`-split, aligned to `Ensembl_transcriptid`) | unsigned | `*1e6 -> u32` |
| AlphaMissense class | `AlphaMissense_pred` | categorical: `B`/`A`/`P` (benign / ambiguous / pathogenic) | Yes | categorical | string→index |
| **ESM1b** | `ESM1b_score` | log-likelihood ratio, typically **negative** (~ -25 … +2) | **Yes** (per Ensembl protein) | **SIGNED** | `*1e6 -> i32 -> zigzag -> u32` |
| ESM1b class | `ESM1b_pred` | categorical (`D`/`T` damaging/tolerated) | Yes | categorical | string→index |
| **REVEL** | `REVEL_score` | 0.0–1.0 | single value per variant in most versions (may be `;`-split in some) — treat as broadcast if cardinality 1 | unsigned | `*1e6 -> u32` |

Optional rankscore columns (`*_rankscore`, all 0–1, unsigned) exist for each and are convenient for
thresholding; skip unless you want them.

## Multi-value / alignment rules

1. Missing value token is `.` (a literal dot). Treat `.` **and** empty string as `None`.
2. Per-transcript columns are `;`-separated and **positionally aligned** with `Ensembl_transcriptid`.
   Split all of them on `;` and zip by index.
3. A per-transcript score cell can itself be `.` for an individual transcript even when others have
   values — keep the alignment; emit `None` for that transcript only.
4. If a score column has cardinality 1 but `Ensembl_transcriptid` has N>1 entries, **broadcast** the
   single value to all transcripts (this is how dbNSFP represents per-variant scores like REVEL).
   The module handles this.

## Why signedness is not optional

fastSA v2 stores values as `u32`. The documented encodings are: unsigned `value*multiplier -> u32`,
**zigzag** for signed values, and string→index for categoricals. ESM1b is an LLR centered below
zero. If you push a value of `-12.3` through the unsigned path you get `(-12.3 * 1e6) as u32` which
wraps to ~4.28e9 — a silently corrupt score that will look like an extreme pathogenic hit. Use
`quantize_signed` (scale → `i32` → zigzag → `u32`) for ESM1b. The Rust module makes this the only
way to encode that field.

## Decode side (so downstream R gets real numbers)

Whatever multiplier you encode with must be applied in reverse when the value is written to the
tab/CSQ/JSON output:
`f32 = (zigzag_decode(u32_or_unsigned) as f32) / multiplier`. Keep the multiplier per-field in one
place (`SCORE_FIELDS` in `dbnsfp_scores.rs`) so encode and decode never drift.
