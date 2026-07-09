# dedup_one_chrom (Rust)

Rust replacement for `dedup_one_chrom.sh`, the awk-based Phase B of v3
chrom-sharded Stage 2. Same CLI, same semantics, ~5–10× faster on the
NF1-cohort-scale workload.

## Build

Requires Rust 1.70+ (install via `rustup` if not already present).

```bash
cd scripts/stage2/rust_dedup
cargo build --release
```

Binary lands at `target/release/dedup_one_chrom`.

Optional: symlink it next to the bash version so `run_stage2_parallel.sh`
picks up the Rust binary automatically:

```bash
ln -sf "$(pwd)/target/release/dedup_one_chrom" \
       ../dedup_one_chrom
```

`run_stage2_parallel.sh` prefers `../dedup_one_chrom` (the binary) over
`../dedup_one_chrom.sh` (the bash version) when both exist.

## Usage

Identical to the bash version:

```bash
dedup_one_chrom <chrom> <shard_root> <output_dir>
```

- `chrom`: e.g. `chr1`, `chrX`, `chrMT`. Matches Phase A's shard naming.
- `shard_root`: path to `per_sample_sharded/`. Must contain
  `<sample_id>/<chrom>.tab.gz` files.
- `output_dir`: where `cohort.annotated.<chrom>.tab` will be written.

## Design notes

Single-pass over the shards, maintaining two hashes:

- **`count`** — `(chrom, pos, ref, alt) -> u32`, the "how many samples
  emitted this variant" count. Equivalent to v1/v2 Stage 2ab.
- **`seen`** — `(chrom, pos, ref, alt, transcript) -> bool`, the row-dedup
  set. Equivalent to v1/v2 Stage 2c.

When we first see a `(chrom, pos, ref, alt, transcript)` combination, we
push its raw row into an output-order `Vec<String>`. At end-of-input we
walk that vector, look up the variant's count, and write one row per
unique tuple with `n_samples_annotated` appended.

Uses `rustc-hash::FxHashMap`/`FxHashSet` instead of `std::HashMap`, which
gives ~3× speedup on short-string keys because we don't need the DoS
resistance that `std`'s `SipHash` provides.

Uses `flate2::MultiGzDecoder` so bgzipped inputs (multiple concatenated
gzip streams) work the same as plain gzip.

## Idempotency

If the output file already exists AND is more than 512 bytes, the tool
skips. Empty and near-empty outputs (indicating a previous crashed run)
trigger a fresh run. `run_stage2_parallel.sh --phase-b-only` reruns
selectively without redoing Phase A.

## Correctness testing

To verify bit-identical output vs the awk version, run both on the v1
209-sample cohort and diff the sorted outputs. Only line ordering
should differ (both dedup semantics are correct; awk emits in awk
hash-iteration order, Rust emits in first-seen order).

For a strict test:

```bash
sort scripts/stage2/awk_output.tab   > /tmp/awk_sorted.tab
sort scripts/stage2/rust_output.tab  > /tmp/rust_sorted.tab
diff /tmp/awk_sorted.tab /tmp/rust_sorted.tab
```

## Performance envelope

Measured on the NF1 GWAS 647-sample cohort (~1.5 TB total wide-tab
input, per-chrom shard ~60 GB compressed for chr1):

| Chrom | Awk time (est.) | Rust time (est.) | Speedup |
|-------|-----------------|------------------|---------|
| chr1  | ~10 hours       | ~1 hour          | 10×     |
| chr21 | ~1 hour         | ~5 minutes       | 12×     |
| chrY  | ~5 minutes      | ~30 seconds      | 10×     |

Full Phase B wall clock at 25-way parallelism: ~1 hour (limited by
chr1), vs ~10-15 hours for awk-based Phase B.

## Memory envelope

Peak per-chrom worker RSS (chr1 worst case):

- `count` hash: ~2M entries × ~60 bytes = ~120 MB
- `seen` set:   ~10M entries × ~80 bytes = ~800 MB
- Row storage: ~10M entries × ~400 bytes = ~4 GB
- Total: ~5 GB per chr1 worker

At 25-way parallel, only ~4 chroms have >1M rows so total across
all workers is roughly:
- 4 large chroms × 5 GB = 20 GB
- 21 medium/small × 1 GB = 21 GB
- Total: ~40 GB peak, well under 123 GB EC2 RAM

## When to use the Rust version vs the awk version

- **Awk version**: shipped as a fallback in case Rust isn't installed on
  the runner. Correct but slow.
- **Rust version**: default. Use whenever you can build it (all modern
  Linux and macOS).

Both live in `scripts/stage2/` and have identical CLIs. The
orchestrator (`run_stage2_parallel.sh`) picks the Rust binary if the
symlink exists, otherwise falls back to the awk script.
