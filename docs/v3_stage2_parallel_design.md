# v3 Stage 2 parallelization — design

## Problem

v1/v2 Stage 2 (cohort merge) is a single-threaded awk hash-dedup that streams
all per-sample wide tabs sequentially. It has TWO sub-passes over the whole
input:

- **2ab** — builds `count[chrom,pos,ref,alt] = n_samples_annotated`, dumps at END
- **2c**  — builds `seen[chrom,pos,ref,alt,transcript]`, streams-writes unique rows

At v1 scale (209 samples, ~325 GB input) this ran in ~30 minutes total. At
v2/v3 scale (647 samples, ~1.5 TB input) the 209-sample-tuned single-threaded
implementation extrapolates to **~30 hours** — and that's what we're
observing in production right now on the NF1 GWAS run.

The awk pattern has three inherent limits:

1. **Single-threaded** — cannot use all 16 cores on the EC2 instance
2. **In-memory hash** — RSS grows monotonically as unique variants accumulate;
   at some point cache pressure slows insertion
3. **Two full passes** — 2ab and 2c both stream every wide tab from disk,
   doubling I/O

## Design

Shard the work by chromosome. Chromosomes are natural boundaries because:

- Every wide-tab row has a chromosome as `$1` (post-Stage 0 normalization)
- Variants on different chromosomes are independent for dedup purposes
- 25 primary contigs (1-22, X, Y, MT) parallelize nicely onto 16 cores
- Chromosome-local hashes are ~1/25th the size of the whole-cohort hash

Two-phase implementation:

### Phase A: split per-sample wide tabs by chromosome

**Parallel across samples** (16-way). For each sample, one pass over the wide
tab writes a family of chrom-specific gzipped shards:

```
per_sample/NG104C2PK8.annotated.tab.gz
  ->  per_sample_sharded/NG104C2PK8/chr1.tab.gz
      per_sample_sharded/NG104C2PK8/chr2.tab.gz
      ...
      per_sample_sharded/NG104C2PK8/chrMT.tab.gz
```

Header written to each shard. Total output size ~= input size (still gzipped,
same content just partitioned). Awk memory usage is bounded by the number of
open output pipes (25). Runtime is I/O-bound; expect ~5-10 min per sample
serial, ~30-60 min at 16-way parallel across all 647 samples.

### Phase B: per-chrom cohort dedup

**Parallel across chromosomes** (25-way, we can oversubscribe past 16 cores
because 2ab and 2c are largely I/O-bound). For each chromosome:

- Run the existing 2ab awk on that chrom's shards from all 647 samples
- Run the existing 2c awk on the same shards
- Run 2d (attach n_samples_annotated) inline
- Emit `cohort.annotated.chr<N>.tab`

Each per-chrom job reads only that chromosome's shards — much less I/O than
the whole cohort. Hash size is ~1/25th of whole-cohort hash, comfortably
fits per worker.

### Phase C: concatenate

Single-threaded but trivial. Header from `cohort.annotated.chr1.tab`, then
concat all chroms:

```
cohort.annotated.chr1.tab  +  chr2  +  ...  +  chrMT  =  cohort.annotated.tab
```

Chromosome-order sort is a nice bonus for downstream tools.

## Expected performance at 647-sample scale

| Phase | Time | Notes |
|-------|------|-------|
| A (split) | ~30-60 min | I/O-bound, 16-way parallel across samples |
| B (per-chrom dedup) | ~30-90 min | I/O-bound, 25-way parallel across chroms; longest chrom (chr1) determines wall-clock |
| C (concat) | ~10-20 min | Serial cat, disk-bound |
| **Total Stage 2** | **~1-3 hours** | vs ~30 hours single-threaded |

## Compatibility with existing pipeline

The output `cohort.annotated.tab` has identical schema and content to what
v2 produces. Downstream stages (3, 3.5, 4, 5) don't need any changes. The
new logic slots into `run_cohort_stages.sh` as a replacement for Stage 2
only.

## Implementation plan

1. `scripts/stage2/split_wide_tabs_by_chrom.sh` — Phase A. Per-sample, called
   in parallel by `run_stage2_parallel.sh`.
2. `scripts/stage2/dedup_one_chrom.sh` — Phase B. Per-chrom, called in
   parallel across chromosomes.
3. `scripts/stage2/run_stage2_parallel.sh` — orchestrator: calls Phase A
   across samples, then Phase B across chroms, then concat.
4. `scripts/run_cohort_stages_v3.sh` — new wrapper that invokes
   `run_stage2_parallel.sh` and then delegates to v1 `cohort_pipeline.sh`
   with `--skip-annotate --skip-merge` for stages 3+.
5. Update `scripts/cohort_pipeline_batched.sh` to optionally use
   `run_cohort_stages_v3.sh` via a new `--stage2-parallel` flag.

## Migration and testing plan

1. Implement on 209-sample v1 cohort first (already have expected outputs
   to compare against — validate identical `cohort.variants.tsv` at end)
2. Run in parallel with the currently-running v2 pipeline once it completes,
   using same inputs, to confirm identical output
3. Cut v3 cohort_pipeline_batched.sh over to parallel Stage 2 as default
4. Keep v2 (single-threaded) available via `--stage2-serial` for
   troubleshooting

## Open questions

- **File-descriptor limit for Phase A.** Splitting into 25 gzip pipes per
  sample = 25 open fds per worker × 16 workers = 400 fds. Well under default
  ulimit (1024) but worth documenting.
- **Genotype TSV handling.** Do we shard genotypes too? Probably not for v3.1;
  stage 2e (genotype concat) is already fast (~few hours) and doesn't
  benefit as much from chrom sharding.
- **Downstream file sharding.** Should `cohort.variants.tsv` also be sharded
  by chromosome? Not for v3.1. Add later if downstream analysis benefits.
