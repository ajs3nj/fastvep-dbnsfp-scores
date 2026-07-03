# v2 Roadmap — 2-week sprint

Planning document for the v2 NF1 cohort run (~500 samples). 2-week sprint
window; SpliceAI deferred to v3 to fit timeline; primary architectural ask
is batched VCF processing so disk can be reclaimed between batches.

---

## 1. Goals

### 1.1 What v2 must deliver in 2 weeks

In priority order:

1. **Batched per-sample annotation with disk reclamation.** Download a batch
   of VCFs, annotate, delete the VCFs from local disk, then move on to the
   next batch. This is the central architectural change; everything else is
   subordinate to it. Allows the 500-sample run to fit on a tractable disk
   budget (one batch at a time) and lets the user clean up as they go.

2. **Bake in the v1 fixes that need a fresh annotation run to take effect:**
   - ClinVar REVIEW_STATUS parser fix (`csq_to_wide_tab.py` — already in
     repo, but cohort tables built before it have silently-zero stars)
   - VEP/VCF indel normalization (`bcftools norm -m-` before Stage 1)
   - SYMBOL-keyed gene-level projection fix (already in repo)

3. **Tighten the pLoF Tier 1 bucket via LOFTEE.** Use a post-VEP LOFTEE
   pass joined back into the annotation rather than a full fastVEP Rust
   integration — faster to land, achieves the same biological filter.

4. **Scaling work for 500 samples.** Stage 2 + Stage 3.5 parallelization
   by chromosome, sharded `cohort.genotypes.tsv`. Without this the run is
   uncomfortably slow at 500 samples; with this it fits in a workday.

5. **Ancestry-stratified popmax if metadata is available.** Big modifier
   sensitivity recovery if it can land. If metadata isn't ready in time,
   defer to v3 — don't block the 2-week run.

### 1.2 Explicitly deferred to v3 (won't fit in 2 weeks)

- **SpliceAI restoration.** Per user direction. Non-coding tiering remains
  limited; rare non-coding variants will continue to land in Tier 4.
  Documented as a known v2 caveat.
- **Conservation scores (PhyloP / GERP) as primary tier evidence.** May
  surface as informational columns if dbNSFP extraction is quick, but no
  tier rule changes.
- **Regulatory annotations** (Ensembl Reg Build, ENCODE).
- **fastVEP-side LOFTEE Rust integration.** Use the post-VEP join-back
  fallback for v2.

### 1.3 Non-goals

- No new variant predictor classes beyond AM + ESM1b.
- No phasing / compound-het resolution.
- No clinical-grade ACMG classification (research-tier mode only).

---

## 2. Batched processing architecture (the headline change)

The driving constraint: at ~500 samples × ~30 GB per VCF = 15 TB of raw VCFs,
none of which fits on the EC2 instance simultaneously. We need to process
batches and reclaim disk after each batch.

### 2.1 Design

```
manifest.tsv = sample_id, batch_id, vcf_url, [metadata...]

For each batch:
  1. Download batch VCFs into $WORK/inputs/batch_<N>/
  2. bcftools norm -m- + chrom-rename (Stage 0 normalization)
  3. fastVEP annotate per sample (Stage 1a)
  4. csq_to_wide_tab.py per sample (Stage 1b)
  5. Compress annotated wide tabs to $WORK/annotated/batch_<N>/*.tab.gz
  6. (Optional) Upload annotated wide tabs to S3 for redundancy
  7. DELETE $WORK/inputs/batch_<N>/ (reclaim disk)
  8. Move to next batch

After all batches complete:
  9. Cohort merge across all batches (Stage 2)
 10. Pre-filter (Stage 3)
 11. Cohort summary + sharded genotypes (Stage 3.5)
 12. R tiering (Stage 4)
 13. Per-gene burden + figures (Stage 5)
```

### 2.2 Idempotency + resumability

The batch loop must be resumable — a single failed batch shouldn't force
restarting the whole pipeline. Mechanism:

- A `$WORK/state/` directory tracks which batches have completed each stage:
  - `state/batch_<N>.downloaded`
  - `state/batch_<N>.normalized`
  - `state/batch_<N>.annotated`
  - `state/batch_<N>.cleaned`
- The batch loop checks for the latest sentinel and skips completed work.
- A `--retry-batch <N>` flag forces re-processing of a specific batch (e.g.,
  if a download was corrupted).
- A `--keep-vcfs` flag disables the disk-reclamation step for debugging.

### 2.3 Batch sizing

Trade-off: bigger batches mean less per-batch overhead but higher peak disk.
Recommended default: **50 samples per batch** for the v2 run.

- 50 samples × ~30 GB VCF = ~1.5 TB peak per batch
- fastVEP-annotated VCFs are ~2 GB each (compressed) → ~100 GB per batch
- Wide tab gzip outputs are ~200 MB each → ~10 GB per batch
- After cleanup of batch <N>: only the ~10 GB wide tabs remain on disk

For a 500-sample run, that's 10 batches with peak disk per batch staying
under 2 TB.

### 2.4 Manifest format

Extends the v1 format:

```tsv
sample_id    batch_id    vcf_url                                     ancestry         severity    family_id    coverage    contamination    qc_pass
NG104C2PK8   1           s3://nf-cohort/raw/NG104C2PK8.vcf.gz        AFR              3.2         FAM-001      32.5        0.012            TRUE
NG1JSQFDDW   1           s3://nf-cohort/raw/NG1JSQFDDW.vcf.gz        NFE              4.1         FAM-002      31.8        0.008            TRUE
...
```

`batch_id` is mandatory in v2 (was implicit in v1). `ancestry`, `severity`,
`family_id`, `coverage`, `contamination`, `qc_pass` columns are recommended;
the pipeline tolerates them being missing but loses corresponding v2
features (ancestry stratification needs `ancestry`; modifier band
calibration needs `severity` for downstream analysis).

---

## 3. v1 retrospective (carryover)

Brief; full detail in v1's `pipeline_methods.md` §7.

### 3.1 What worked

- Class-conditional tiering with per-class evidence paths.
- AM + ESM1b orthogonal-agreement Tier 1 promotion.
- NF anchor-gene override for the four core genes.
- Per-sample annotation parallelism (Stage 1).
- Provenance log + validation suite.

### 3.2 What v2 fixes

- **ClinVar star parser bug** — Python fix in repo; takes effect on fresh
  Stage 1 run. v2 will have real `clin_stars` populated, can re-enable the
  `≥ 1 star` filter in `tier_variants.R`.
- **VEP/VCF indel encoding mismatch** — `bcftools norm -m-` in Stage 0
  normalizes per-sample VCFs to bi-allelic left-anchored before Stage 1.
  Indel carrier joins work after this lands. Figs 03 / 08 will produce
  sensible output without manual fixes.
- **SYMBOL-keyed gene-level projection** — already fixed in
  `csq_to_wide_tab.py`. Fresh Stage 1 run will populate `loeuf` / `pli`
  correctly; no patch script needed.
- **NF override ordering bug** — already fixed in `tier_variants.R`. Carries
  forward unchanged.

### 3.3 What v2 doesn't fix (deferred)

- SpliceAI integration.
- Conservation as primary tier evidence.
- Regulatory annotations.
- Stratified cohort_common rule (depends on ancestry metadata).

---

## 4. Two-week sprint plan

Calendar assumes one engineer working full-time. Adjust if part-time.

### Week 1: Infrastructure + annotation pipeline

| Day | Tasks | Definition of done |
|-----|-------|---------------------|
| Mon | Implement Stage 0 (`bcftools norm -m-` + chrom rename) as a per-sample script. Test on 3 sample VCFs (different ancestries) | Stage 0 output VCFs are bi-allelic, Ensembl chroms, indexed |
| Mon-Tue | Refactor `cohort_pipeline.sh` into a batched orchestrator (`run_batch.sh <batch_id>` per batch + `run_cohort_stages.sh` for post-batch stages). Add state sentinels and resumability | Resumable across batch failures; --retry-batch flag works |
| Wed | Add LOFTEE post-VEP pass: run `ensembl-vep --plugin LoF` on fastVEP-annotated VCFs, parse the LOF column, join `loftee_lof` / `loftee_filter` columns into the wide tab | Per-sample wide tabs have populated `loftee_lof`; HC / LC distribution matches Ensembl VEP+LOFTEE on a 1000-variant validation set within 1% |
| Wed-Thu | Re-enable LOFTEE filter in `tier_variants.R` pLoF Tier 1 rule (already plumbed; just remove the `# currently always TRUE` comment and verify) | Tier 1 pLoF excludes LOFTEE-LC; validation passes |
| Thu | Re-enable ClinVar star filter in `tier_variants.R` (one-line revert: `T0(stars >= 1) | is.na(stars)`). Add validation check that `clin_stars` is not silently zero | Validation passes; check would have caught v1 bug |
| Thu-Fri | Bake disk-reclamation step into batched orchestrator. Test full cycle on 2 batches × 3 samples each | Local disk after batch_2 completes contains only wide tabs, no source VCFs |
| Fri | Update INSTALL.md with batched manifest format + cleanup of v1 documentation around removed cohort_common rule | INSTALL.md describes the new manifest + batched run |

**End of Week 1 deliverable:** End-to-end batched pipeline runs successfully
on a 2-batch × 3-sample test cohort. ClinVar stars populated, LOFTEE
filtering active, indels properly encoded.

### Week 2: Scaling, ancestry, run, hand-off

| Day | Tasks | Definition of done |
|-----|-------|---------------------|
| Mon | Shard `cohort.genotypes.tsv` by chromosome in Stage 3.5. Update `cohort_figures.R` genotype-streaming to handle sharded files | Per-chrom genotype files produced; figs 03 / 06 / 07 / 08 work against them |
| Mon-Tue | Parallelize Stage 2 cohort merge by chromosome. Parallelize Stage 3.5 cohort summary by chromosome | Stage 2 + Stage 3.5 wall-clock < 15 min each on v1 209-sample data |
| Tue | If ancestry metadata is available: implement per-ancestry popmax lookup + ancestry-stratified `cohort_af` columns. If not: skip and document as v3 carry-over | Either: feature lands and works on test, OR: documented decision to defer with clear reasoning |
| Wed | Re-run v1 209-sample data through v2 pipeline as a regression test. Compare tier distributions vs v1 final state | v1 cohort re-run produces tiering output close to v1 final; differences explained by LOFTEE + correct ClinVar stars |
| Wed-Thu | Begin 500-sample run. Batches launch sequentially, each batch 50 samples | First 3-4 batches complete; throughput projections firm |
| Thu-Fri | Continue 500-sample run. Monitor disk reclamation working as expected | Run completes; final tarball produced |
| Fri | Hand-off: cover email + tarball share + final provenance log | Collaborators have outputs |

**End of Week 2 deliverable:** 500-sample tiered cohort with figures + methods
+ tarball, shipped to collaborators.

---

## 5. Engineering tasks in detail

### 5.1 Stage 0 — VCF normalization (new)

```bash
# scripts/stage0_normalize_vcf.sh <input.vcf.gz> <output.vcf.gz>
# Per-sample step that runs in the batch loop before fastVEP annotation.

bcftools norm -m- -f $FASTA $INPUT |
  bcftools annotate --rename-chrs $CHR_REMAP - |
  bgzip > $OUTPUT
tabix -p vcf $OUTPUT
```

Adds ~2-5 min per sample. Output is bi-allelic split, Ensembl-style chroms,
indexed.

### 5.2 Batched orchestrator structure

```
scripts/
├── stage0_normalize_vcf.sh           # per-sample, called by run_batch.sh
├── run_batch.sh                       # per-batch: download → norm → annotate → cleanup
├── run_cohort_stages.sh               # after all batches: stages 2, 3, 3.5, 4, 5
├── cohort_pipeline.sh                 # top-level orchestrator: calls run_batch.sh
│                                      # in a loop, then run_cohort_stages.sh
└── state/                             # batch sentinel files
```

### 5.3 Post-VEP LOFTEE pass

```bash
# Run AFTER fastVEP annotation, BEFORE csq_to_wide_tab.py
vep \
  --input_file $FASTVEP_VCF \
  --plugin LoF,loftee_path:$LOFTEE_DIR,human_ancestor_fa:$HUMAN_ANCESTOR \
  --vcf --output_file $LOFTEE_VCF \
  --no_stats --offline --cache --dir_cache $VEP_CACHE

# Then csq_to_wide_tab.py picks up the LoF column from the new CSQ string
```

Single-threaded VEP call adds ~5 min per sample. Could be parallelized
across cores within a batch.

### 5.4 Per-chromosome sharding for genotypes

`cohort.genotypes.chr1.tsv.gz`, `chr2`, ..., `chrX`, `chrY`, `chrMT`.

Stage 3.5 already splits by variant_key for awk aggregation; just split
the output emit. Stage 4 R tiering doesn't read genotypes directly so
nothing needs to change there. Figures script needs an awk loop that
streams each shard.

### 5.5 Validation additions

`scripts/validate_cohort.sh` gets new checks:

- **clin_stars distribution sanity check** — fail if every variant has
  `clin_stars == 0` (catches the v1 parser bug regression)
- **LOFTEE distribution sanity check** — fail if every pLoF Tier 1 variant
  has empty `loftee_lof` (catches LOFTEE pipeline regression)
- **Indel encoding sanity check** — fail if any cohort.variants.tsv row
  has `alt == "-"` (catches normalization regression)
- **Per-chromosome genotypes coverage** — fail if any chromosome shard is
  missing samples that should be there

---

## 6. Resource estimates

### 6.1 Engineering

| Phase | Person-days |
|-------|-------------|
| Week 1 infrastructure | 5 |
| Week 2 scaling + run + hand-off | 5 |
| **Total** | **10 (full-time, single engineer)** |

### 6.2 Compute / storage at 500 samples

| Resource | Estimate | Notes |
|----------|----------|-------|
| Per-batch peak disk | ~1.5 TB | 50 VCFs at 30 GB each |
| Concurrent annotated tabs (cumulative) | ~100 GB | After all 10 batches; doesn't grow further |
| Stage 1 wall-clock per batch | 30-90 min | Same as v1 per-sample, batches sequential |
| Stage 0 + Stage 1 cumulative (10 batches) | 5-15 hours | Bulk of the run |
| Stage 2 + 3.5 (parallelized) | < 30 min | Down from v1's 60 min |
| Stage 4 R tier | 2-3 min | |
| Total pipeline | 6-18 hours | Mostly Stage 1 |
| EC2 instance | r6i.4xlarge (16 vCPU, 128 GB) | Same as v1; per-batch peak fits |

---

## 7. Open questions for study leads

Need answers before Week 1 starts:

1. **VCF source / download method** — S3, GCS, FTP? Authenticated? Need
   download credentials before batches can run.
2. **Batch ordering preference** — random shuffle (to avoid batch-effect
   bias if some batches happen to fail), or by date/site (to make failure
   recovery easier to reason about)?
3. **Ancestry metadata** — available? Format? If yes, ancestry stratification
   lands in v2; if no, v3.
4. **Phenotype severity scores** — available? Will be needed for downstream
   modifier analysis but doesn't block v2 pipeline.
5. **Family structure** — any related samples? Affects modifier analysis but
   not v2 pipeline.
6. **Failure tolerance** — if 1-2 samples fail Stage 1 (e.g., corrupted
   VCF), should the run halt or proceed without them? Default: proceed
   and document missing samples in provenance.

---

## 8. Risk register

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| 500 VCFs download bandwidth bottleneck | High | Medium | Download in parallel with previous batch's annotation; bigger EC2 network bandwidth instance if needed |
| LOFTEE post-VEP pass slower than expected | Medium | Low | Parallelize across cores; fall back to no-LOFTEE pLoF Tier 1 (v1 baseline) if unsalvageable |
| Sample-level QC issues surface during run (high contamination, low coverage) | Medium | Low | Add a per-sample QC report after Stage 1; defer those samples for re-processing; document in provenance |
| Disk reclamation race condition | Low | High | State sentinels are written atomically (write to .tmp, rename); only delete VCFs when annotation sentinel exists |
| Cohort merge memory blowup at 500 samples | Low | Medium | Chrom-sharded merge keeps per-shard memory bounded; can fall back to disk-backed sort if needed |
| Ancestry metadata arrives partway through Week 2 | Medium | Low | Re-run Stage 4 only after main run completes; no need to redo annotation |
| Indel normalization changes downstream Tier 1 counts vs v1 | Medium | Low | Expected — v1 indels were mis-counted. Document the delta as a v2 improvement, not a regression |

---

## 9. Definition of v2 success

1. End-to-end 500-sample pipeline completes within budget (target: 24 hours
   total wall-clock from first download to final tarball).
2. Local disk never exceeds ~2 TB peak during the run (batched cleanup
   works).
3. Pipeline is resumable — killing it mid-run and restarting picks up where
   it left off.
4. ClinVar `clin_stars` is populated correctly (`stars >= 1` filter active).
5. pLoF Tier 1 LOFTEE-LC fraction < 5%.
6. Indels are properly encoded — fig 03 / 08 show sensible data without
   workarounds.
7. Modifier candidate count is biologically plausible (1,000–10,000+
   candidates across NF/RAS-MAPK genes; not zero like one v1 intermediate
   state).
8. Final tarball includes: cohort.variants.tsv, cohort.genes.tsv, figures
   (all 8), provenance log, methods documentation.
9. Validation suite passes including the new v2 checks (clin_stars sanity,
   LOFTEE sanity, indel encoding sanity).

---

## 10. v3 carry-over (post-v2)

Items deferred for after the 2-week v2 run:

- SpliceAI restoration for non-coding tiering (Tier A from original v2 plan)
- Conservation scores (PhyloP / GERP) as primary tier evidence
- Regulatory annotations (Ensembl Reg Build + ENCODE)
- fastVEP-side LOFTEE Rust integration (replace the post-VEP fallback)
- Ancestry-stratified popmax (if not landed in v2 Week 2)
- Cohort_common rule re-introduction with ancestry awareness
- ENCODE cell-type-specific regulatory overlay for Schwann cells / neural crest
- Phasing / compound-het resolution (requires trio sequencing data)
- Interactive results viewer / dashboard

---

## 11. Next concrete steps

1. **Today / tomorrow:** Schedule 30-minute kickoff with study leads. Get
   answers to §7 open questions. Lock VCF download credentials.
2. **Monday:** Begin Week 1 work. Stage 0 + batched orchestrator first.
3. **End of Week 1:** Test cohort of 2 batches × 3 samples runs end-to-end.
4. **Monday Week 2:** Begin scaling work (genotype sharding, Stage 2/3.5
   parallelization).
5. **Wed Week 2:** Start 500-sample run.
6. **End of Week 2:** Tarball shipped, hand-off email sent.
