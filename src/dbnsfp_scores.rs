//! Self-contained dbNSFP score extraction + fastSA-style quantization.
//!
//! This module is intentionally `std`-only and has NO dependency on fastVEP's internal types, so it
//! compiles and is unit-tested in isolation. Wire it into your fork's `dbnsfp` parser via the thin
//! adapter described in INTEGRATION.md.
//!
//! It extracts AlphaMissense, ESM1b, REVEL and (coding) CADD from a dbNSFP v4.5+ TSV row, handling:
//!   * version-safe column lookup by header name (never fixed index),
//!   * `;`-separated per-transcript values aligned to `Ensembl_transcriptid` (with broadcast of
//!     single-cardinality columns such as REVEL/CADD),
//!   * the `.` / empty missing-value token,
//!   * fastSA v2 quantization: `value * multiplier -> u32`, with zigzag for SIGNED scores
//!     (ESM1b, CADD raw) so negative values do not wrap.

use std::collections::HashMap;

/// Per-field encoding metadata. Keep encode and decode using the SAME multiplier from here so they
/// can never drift. `signed = true` routes through scale -> i32 -> zigzag -> u32.
#[derive(Clone, Copy, Debug)]
pub struct ScoreField {
    pub name: &'static str,
    pub multiplier: f32,
    pub signed: bool,
}

/// The four scores we add via the dbNSFP route, plus their encoding choices.
/// (See docs/dbnsfp_score_columns.md for rationale.)
pub const SCORE_FIELDS: &[ScoreField] = &[
    ScoreField { name: "alphamissense", multiplier: 1e6, signed: false },
    ScoreField { name: "esm1b",         multiplier: 1e6, signed: true  }, // LLR, negative
    ScoreField { name: "revel",         multiplier: 1e6, signed: false },
    ScoreField { name: "cadd_phred",    multiplier: 1e3, signed: false },
    ScoreField { name: "cadd_raw",      multiplier: 1e6, signed: true  }, // can be negative
];

pub fn score_field(name: &str) -> Option<&'static ScoreField> {
    SCORE_FIELDS.iter().find(|f| f.name == name)
}

// ---------------------------------------------------------------------------
// Quantization (matches fastSA v2: value*multiplier -> u32, zigzag for signed)
// ---------------------------------------------------------------------------

/// Zigzag-encode a signed 32-bit integer into u32 (same scheme as protobuf / echtvar).
#[inline]
pub fn zigzag_encode(v: i32) -> u32 {
    ((v << 1) ^ (v >> 31)) as u32
}

/// Inverse of `zigzag_encode`.
#[inline]
pub fn zigzag_decode(v: u32) -> i32 {
    ((v >> 1) as i32) ^ -((v & 1) as i32)
}

/// Encode an unsigned score (range >= 0) as `round(value * multiplier)` clamped into u32.
#[inline]
pub fn quantize_unsigned(value: f32, multiplier: f32) -> u32 {
    let scaled = (value * multiplier).round();
    if scaled < 0.0 {
        0
    } else if scaled > u32::MAX as f32 {
        u32::MAX
    } else {
        scaled as u32
    }
}

/// Encode a SIGNED score (may be negative) as zigzag(round(value * multiplier)).
#[inline]
pub fn quantize_signed(value: f32, multiplier: f32) -> u32 {
    let scaled = (value * multiplier).round();
    let clamped = scaled.clamp(i32::MIN as f32, i32::MAX as f32) as i32;
    zigzag_encode(clamped)
}

/// Encode using a named `ScoreField`'s rules. This is the function the adapter should call so that
/// signedness is never chosen by hand at the call site.
pub fn quantize(field: &ScoreField, value: f32) -> u32 {
    if field.signed {
        quantize_signed(value, field.multiplier)
    } else {
        quantize_unsigned(value, field.multiplier)
    }
}

/// Decode back to f32 for output (tab/CSQ/JSON). Mirrors `quantize`.
pub fn dequantize(field: &ScoreField, stored: u32) -> f32 {
    if field.signed {
        zigzag_decode(stored) as f32 / field.multiplier
    } else {
        stored as f32 / field.multiplier
    }
}

// ---------------------------------------------------------------------------
// Header indexing (version-safe: look up columns by name)
// ---------------------------------------------------------------------------

/// Maps dbNSFP column name -> field index, built once from the `#chr ... ` header line.
pub struct HeaderIndex {
    map: HashMap<String, usize>,
}

impl HeaderIndex {
    /// Build from the tab-separated header line (with or without a leading '#').
    pub fn from_header(header_line: &str) -> Self {
        let mut map = HashMap::new();
        for (i, col) in header_line.trim_end_matches(['\n', '\r']).split('\t').enumerate() {
            map.insert(col.to_string(), i);
        }
        HeaderIndex { map }
    }

    pub fn get<'a>(&self, fields: &'a [&'a str], col: &str) -> Option<&'a str> {
        self.map.get(col).and_then(|&i| fields.get(i)).copied()
    }

    pub fn has(&self, col: &str) -> bool {
        self.map.contains_key(col)
    }
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

const MISSING: &str = ".";

/// Parse a possibly-missing float cell. `.` and empty are None.
fn parse_opt_f32(cell: &str) -> Option<f32> {
    let c = cell.trim();
    if c.is_empty() || c == MISSING {
        None
    } else {
        c.parse::<f32>().ok()
    }
}

/// Split a `;`-separated dbNSFP cell into Option<&str> per element.
fn split_multi(cell: &str) -> Vec<Option<&str>> {
    let c = cell.trim();
    if c.is_empty() || c == MISSING {
        return vec![];
    }
    c.split(';')
        .map(|e| {
            let e = e.trim();
            if e.is_empty() || e == MISSING { None } else { Some(e) }
        })
        .collect()
}

/// One transcript's protein-level scores.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TranscriptScore {
    pub transcript_id: Option<String>,
    pub protein_id: Option<String>,
    pub alphamissense: Option<f32>,
    pub alphamissense_class: Option<String>,
    pub esm1b: Option<f32>,
    pub esm1b_class: Option<String>,
}

/// All scores parsed from one dbNSFP variant row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VariantScores {
    pub chrom: String,
    pub pos: u32,
    pub reference: String,
    pub alt: String,
    /// Per-variant (broadcast) scores.
    pub revel: Option<f32>,
    pub cadd_phred: Option<f32>,
    pub cadd_raw: Option<f32>,
    /// Per-transcript scores, aligned to `Ensembl_transcriptid`.
    pub transcripts: Vec<TranscriptScore>,
}

/// Helper: return the i-th element of a split column, or None if absent.
fn nth<'a>(v: &'a [Option<&'a str>], i: usize) -> Option<&'a str> {
    v.get(i).copied().flatten()
}

/// Parse a single dbNSFP data row (already split into &str fields) into `VariantScores`.
///
/// Columns are resolved by name via `idx`, so this is robust across dbNSFP versions. Any column that
/// is absent in this dbNSFP build is simply treated as missing.
pub fn parse_row(fields: &[&str], idx: &HeaderIndex) -> VariantScores {
    let get = |c: &str| idx.get(fields, c).unwrap_or(MISSING);

    let chrom = get("#chr").to_string();
    // Parse position as an integer directly. Do NOT route through f32: positions exceed 2^24
    // (e.g. chr1 is ~249 Mb) and f32 cannot represent integers that large exactly, so an f32
    // round-trip would silently corrupt the coordinate.
    let pos = {
        let c = get("pos(1-based)").trim();
        if c.is_empty() || c == MISSING { 0 } else { c.parse::<u32>().unwrap_or(0) }
    };
    let reference = get("ref").to_string();
    let alt = get("alt").to_string();

    // Per-variant (single-value) scores.
    let revel = parse_opt_f32(get("REVEL_score"));
    let cadd_phred = parse_opt_f32(get("CADD_phred"));
    let cadd_raw = parse_opt_f32(get("CADD_raw"));

    // Per-transcript columns.
    let tx_ids = split_multi(get("Ensembl_transcriptid"));
    let prot_ids = split_multi(get("Ensembl_proteinid"));
    let am_scores = split_multi(get("AlphaMissense_score"));
    let am_preds = split_multi(get("AlphaMissense_pred"));
    let esm_scores = split_multi(get("ESM1b_score"));
    let esm_preds = split_multi(get("ESM1b_pred"));

    // Determine how many transcripts this row describes.
    let n = [
        tx_ids.len(),
        prot_ids.len(),
        am_scores.len(),
        am_preds.len(),
        esm_scores.len(),
        esm_preds.len(),
    ]
    .into_iter()
    .max()
    .unwrap_or(0);

    // Broadcast helper: a length-1 per-transcript column applies to all transcripts.
    let pick = |v: &[Option<&str>], i: usize| -> Option<String> {
        if v.len() == 1 {
            v[0].map(|s| s.to_string())
        } else {
            nth(v, i).map(|s| s.to_string())
        }
    };
    let pick_f32 = |v: &[Option<&str>], i: usize| -> Option<f32> {
        let cell = if v.len() == 1 { v.first().copied().flatten() } else { nth(v, i) };
        cell.and_then(parse_opt_f32)
    };

    let mut transcripts = Vec::with_capacity(n);
    for i in 0..n {
        transcripts.push(TranscriptScore {
            transcript_id: pick(&tx_ids, i),
            protein_id: pick(&prot_ids, i),
            alphamissense: pick_f32(&am_scores, i),
            alphamissense_class: pick(&am_preds, i),
            esm1b: pick_f32(&esm_scores, i),
            esm1b_class: pick(&esm_preds, i),
        });
    }

    VariantScores { chrom, pos, reference, alt, revel, cadd_phred, cadd_raw, transcripts }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zigzag_roundtrips_including_negatives() {
        for v in [-2_000_000_i32, -123, -1, 0, 1, 42, 2_000_000] {
            assert_eq!(zigzag_decode(zigzag_encode(v)), v);
        }
    }

    #[test]
    fn signed_encoding_does_not_wrap_negative_esm1b() {
        let f = score_field("esm1b").unwrap();
        // ESM1b LLR of -12.3 must NOT become a huge positive integer.
        let enc = quantize(f, -12.3);
        let dec = dequantize(f, enc);
        assert!((dec - (-12.3)).abs() < 1e-3, "decoded {dec}");
        // The unsigned path would corrupt it:
        assert_ne!(enc, quantize_unsigned(-12.3, f.multiplier));
    }

    #[test]
    fn unsigned_alphamissense_roundtrips() {
        let f = score_field("alphamissense").unwrap();
        let dec = dequantize(f, quantize(f, 0.9871));
        assert!((dec - 0.9871).abs() < 1e-5, "decoded {dec}");
    }

    #[test]
    fn cadd_raw_negative_roundtrips() {
        let f = score_field("cadd_raw").unwrap();
        let dec = dequantize(f, quantize(f, -2.345));
        assert!((dec - (-2.345)).abs() < 1e-4, "decoded {dec}");
    }

    fn idx_and_row() -> (HeaderIndex, Vec<&'static str>) {
        // Minimal header subset in arbitrary order to prove name-based lookup.
        let header = "#chr\tpos(1-based)\tref\talt\taaref\taaalt\tEnsembl_transcriptid\tEnsembl_proteinid\tREVEL_score\tCADD_phred\tCADD_raw\tAlphaMissense_score\tAlphaMissense_pred\tESM1b_score\tESM1b_pred\tgenename";
        let idx = HeaderIndex::from_header(header);
        // Two transcripts; REVEL/CADD single (broadcast); ESM1b negative; second tx AM missing.
        let row = vec![
            "17", "43045712", "C", "T", "R", "Q",
            "ENST00000357654;ENST00000468300",
            "ENSP00000350283;ENSP00000417148",
            "0.842",          // REVEL (single)
            "24.7",           // CADD_phred (single)
            "3.21",           // CADD_raw (single, positive here)
            "0.9912;.",       // AlphaMissense per-tx, 2nd missing
            "P;.",            // AlphaMissense_pred
            "-8.4;-2.1",      // ESM1b per-tx (negative)
            "D;T",            // ESM1b_pred
            "BRCA1;BRCA1",
        ];
        (idx, row)
    }

    #[test]
    fn parses_per_transcript_alignment_and_broadcast() {
        let (idx, row) = idx_and_row();
        let vs = parse_row(&row, &idx);

        assert_eq!(vs.chrom, "17");
        assert_eq!(vs.pos, 43045712);
        assert_eq!(vs.revel, Some(0.842));
        assert_eq!(vs.cadd_phred, Some(24.7));
        assert_eq!(vs.transcripts.len(), 2);

        let t0 = &vs.transcripts[0];
        assert_eq!(t0.transcript_id.as_deref(), Some("ENST00000357654"));
        assert_eq!(t0.alphamissense, Some(0.9912));
        assert_eq!(t0.alphamissense_class.as_deref(), Some("P"));
        assert_eq!(t0.esm1b, Some(-8.4));

        let t1 = &vs.transcripts[1];
        assert_eq!(t1.transcript_id.as_deref(), Some("ENST00000468300"));
        assert_eq!(t1.alphamissense, None); // was "."
        assert_eq!(t1.esm1b, Some(-2.1));
    }

    #[test]
    fn large_position_is_not_corrupted_by_f32() {
        // 200_000_007 is > 2^24 and not a multiple of the f32 step at that magnitude, so an f32
        // round-trip would change it. Parsing as u32 must return it exactly.
        let header = "#chr\tpos(1-based)\tref\talt\tEnsembl_transcriptid";
        let idx = HeaderIndex::from_header(header);
        let row = vec!["1", "200000007", "A", "G", "ENST00000000001"];
        let vs = parse_row(&row, &idx);
        assert_eq!(vs.pos, 200_000_007);
    }

    #[test]
    fn missing_value_token_is_none() {
        assert_eq!(parse_opt_f32("."), None);
        assert_eq!(parse_opt_f32(""), None);
        assert_eq!(parse_opt_f32("  .  "), None);
        assert_eq!(parse_opt_f32("0.5"), Some(0.5));
    }
}
