// dedup_one_chrom
//
// v3 Stage 2 Phase B: per-chromosome cohort dedup.
//
// Reads gzipped per-sample wide-tab shards for one chromosome (produced by
// v3 Phase A `split_wide_tabs_by_chrom.sh`), emits a single deduplicated
// tab-separated file with n_samples_annotated attached.
//
// Semantics identical to v1/v2 Stage 2 (Stages 2ab + 2c + 2d combined):
//   - Dedup rows by (chrom, pos, ref, alt, transcript)  (keeps the first
//     row seen for each such tuple, matches v1 awk behavior)
//   - Count how many samples emitted each (chrom, pos, ref, alt), attach
//     that as a new n_samples_annotated column at the end
//
// Usage:
//     dedup_one_chrom <chrom> <shard_root> <output_dir>
//
//   chrom:        e.g. "chr1", "chrX", "chrMT"  (matches Phase A naming)
//   shard_root:   e.g. /data/nf1/outputs/per_sample_sharded
//                 must contain <shard_root>/<sample_id>/<chrom>.tab.gz
//   output_dir:   where to write cohort.annotated.<chrom>.tab
//
// Design notes:
//   - Single-pass over the shards. We hold two hashes:
//       count: (chrom, pos, ref, alt) -> u32          (Stage 2ab equivalent)
//       seen:  (chrom, pos, ref, alt, transcript) -> row_index
//     For each row we increment count if it's the first time we've seen the
//     4-tuple in the current sample, and store the row in a Vec if it's the
//     first time we've seen the 5-tuple overall.
//   - The `count` hash is small (~2M entries for chr1). The `seen` +
//     stored-rows footprint scales with unique (variant, transcript) pairs
//     — for chr1 across 647 samples we estimate ~5-15M rows at ~200-500
//     bytes each, so ~1-8 GB per chr1 worker. Well within budget on the
//     123 GB EC2 instance even at 25-way parallelism.
//   - We use rustc_hash for both hashes: FxHash is ~3x faster than the
//     default SipHash for our small-string keys.

use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};

// Compact interned key types. We store owned Strings for the variable-length
// parts (ref, alt, transcript) and pack chrom+pos into a compact form.
// For cohort scale we're mostly bottlenecked by hash lookups rather than
// bytes-per-key, so this simple representation is fine.
type Variant = (String, String, String, String); // chrom, pos, ref, alt
type VariantTranscript = (Variant, String);

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: {} <chrom> <shard_root> <output_dir>", args[0]);
        std::process::exit(1);
    }
    let chrom = &args[1];
    let shard_root = PathBuf::from(&args[2]);
    let output_dir = PathBuf::from(&args[3]);

    fs::create_dir_all(&output_dir)?;
    let output_path = output_dir.join(format!("cohort.annotated.{}.tab", chrom));

    // Idempotency check: skip if output already exists and is non-trivially large.
    if output_path.exists() {
        if let Ok(meta) = fs::metadata(&output_path) {
            if meta.len() > 512 {
                eprintln!("[dedup {}] output already present ({} bytes); skipping",
                          chrom, meta.len());
                return Ok(());
            }
        }
    }

    // Discover shards: <shard_root>/<sample_id>/<chrom>.tab.gz
    let shards = discover_shards(&shard_root, chrom)?;
    eprintln!("[dedup {}] found {} sample shards", chrom, shards.len());

    if shards.is_empty() {
        // Emit an empty-but-headed output so downstream concat doesn't
        // fail if e.g. chrY has no shards in an all-female cohort.
        emit_empty_with_header(&shard_root, chrom, &output_path)?;
        return Ok(());
    }

    // Read the header from the first shard to determine the transcript column.
    let header = read_header(&shards[0])?;
    let field_names: Vec<&str> = header.split('\t').collect();
    let transcript_col = field_names
        .iter()
        .position(|&f| f == "transcript")
        .ok_or_else(|| anyhow!("no 'transcript' column in header of {:?}", shards[0]))?;
    eprintln!("[dedup {}] header has {} columns; transcript at index {}",
              chrom, field_names.len(), transcript_col);

    // Data structures.
    let mut count: FxHashMap<Variant, u32> = FxHashMap::default();
    let mut seen: FxHashSet<VariantTranscript> = FxHashSet::default();
    let mut rows: Vec<(Variant, String)> = Vec::new(); // (variant_key, raw row) for output

    let mut total_in = 0u64;
    let mut per_shard_in = 0u64;
    for (shard_idx, shard) in shards.iter().enumerate() {
        let reader = open_gz(shard)
            .with_context(|| format!("opening shard {:?}", shard))?;
        let mut lines = reader.lines();
        let _ = lines.next(); // skip header

        let mut sample_seen: FxHashSet<Variant> = FxHashSet::default();
        per_shard_in = 0;

        for line in lines {
            let line = line?;
            per_shard_in += 1;

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() <= transcript_col {
                continue; // malformed row -- skip defensively
            }

            let v: Variant = (
                fields[0].to_string(),
                fields[1].to_string(),
                fields[2].to_string(),
                fields[3].to_string(),
            );

            // Per-sample count: only increment count[v] once per sample
            if sample_seen.insert(v.clone()) {
                *count.entry(v.clone()).or_insert(0) += 1;
            }

            // Global dedup on (variant, transcript)
            let vt: VariantTranscript = (v.clone(), fields[transcript_col].to_string());
            if seen.insert(vt) {
                rows.push((v, line));
            }
        }
        total_in += per_shard_in;

        if (shard_idx + 1) % 100 == 0 || shard_idx + 1 == shards.len() {
            eprintln!("[dedup {}] processed {}/{} shards ({} rows in, {} unique so far)",
                      chrom, shard_idx + 1, shards.len(), total_in, rows.len());
        }
    }

    // Emit output.
    let mut writer = BufWriter::new(File::create(&output_path)?);
    writeln!(writer, "{}\tn_samples_annotated", header)?;
    for (v, row) in &rows {
        let n = count.get(v).copied().unwrap_or(0);
        writeln!(writer, "{}\t{}", row, n)?;
    }
    writer.flush()?;

    eprintln!("[dedup {}] wrote {} unique variant-transcript rows to {:?}",
              chrom, rows.len(), output_path);
    Ok(())
}

fn discover_shards(shard_root: &Path, chrom: &str) -> Result<Vec<PathBuf>> {
    let mut shards = Vec::new();
    let filename = format!("{}.tab.gz", chrom);
    for entry in fs::read_dir(shard_root)
        .with_context(|| format!("reading shard root {:?}", shard_root))?
    {
        let entry = entry?;
        let sample_dir = entry.path();
        if sample_dir.is_dir() {
            let shard = sample_dir.join(&filename);
            if shard.exists() && fs::metadata(&shard)?.len() > 0 {
                shards.push(shard);
            }
        }
    }
    shards.sort(); // deterministic ordering across runs
    Ok(shards)
}

fn open_gz(path: &Path) -> Result<BufReader<MultiGzDecoder<File>>> {
    let file = File::open(path).with_context(|| format!("opening {:?}", path))?;
    // MultiGzDecoder handles concatenated gzip streams (bgzip). BufReader on
    // top gives us line iteration and a decent read buffer.
    Ok(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
}

fn read_header(shard: &Path) -> Result<String> {
    let mut reader = open_gz(shard)?;
    let mut header = String::new();
    reader.read_line(&mut header)?;
    Ok(header.trim_end_matches('\n').to_string())
}

fn emit_empty_with_header(shard_root: &Path, chrom: &str, output_path: &Path) -> Result<()> {
    // Find ANY shard (in any sample, any chrom) so we can steal a header row.
    // If no shards exist at all, just write a stub header we know from schema.
    for entry in fs::read_dir(shard_root)? {
        let entry = entry?;
        let sample_dir = entry.path();
        if !sample_dir.is_dir() {
            continue;
        }
        for entry2 in fs::read_dir(&sample_dir)? {
            let shard = entry2?.path();
            if shard.extension().and_then(|s| s.to_str()) == Some("gz") {
                if let Ok(header) = read_header(&shard) {
                    let mut writer = BufWriter::new(File::create(output_path)?);
                    writeln!(writer, "{}\tn_samples_annotated", header)?;
                    eprintln!("[dedup {}] no shards found; wrote empty output with header",
                              chrom);
                    return Ok(());
                }
            }
        }
    }
    // Fall back: create the output file empty
    File::create(output_path)?;
    eprintln!("[dedup {}] no shards + no header source; wrote empty output", chrom);
    Ok(())
}
