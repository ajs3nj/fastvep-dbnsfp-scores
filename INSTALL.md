# Installation and operation

This document covers everything needed to actually *run* the pipeline. If you're here to
*interpret* output, read `README.md` first — that's where the tier definitions,
`tier_reason` taxonomy, and recipe-style filtering recommendations live.

---

## 1. System dependencies

Tested on Ubuntu 22.04 LTS (AWS EC2 r6i family, 64 GB RAM minimum for a 200-sample cohort
WGS run; 256 GB recommended). All dependencies are open-source and conda-installable.

### Required

| Tool | Version | Notes |
|------|---------|-------|
| fastVEP (Huang-lab fork) | 0.x | Built from source with our dbNSFP patch — see §2 |
| Rust + cargo | 1.70+ | Only needed to build fastVEP |
| Python | 3.10+ | For `csq_to_wide_tab.py` and pipeline glue |
| R | 4.2+ | For `tier_variants.R`, `cohort_figures.R` |
| R packages | data.table, ggplot2, scales | `install.packages(c("data.table","ggplot2","scales"))` |
| bcftools | 1.18+ | Chromosome renaming, VCF normalization |
| samtools | 1.18+ | FASTA `.fai` indexing |
| GNU parallel OR xargs | recent | Per-sample annotation parallelism |
| tabix / bgzip | htslib 1.18+ | VCF indexing |
| awk (gawk) | 4.x+ | Cohort merge and Stage 3.5 aggregation. Mawk works but gawk is preferred for `SUBSEP` semantics |
| synapseclient (Python pkg) | 4.x+ | Only if manifest has `syn*` URLs. `pip install synapseclient`. Requires `SYNAPSE_AUTH_TOKEN` env var or `~/.synapseConfig` (see §4.4) |
| aws-cli | 2.x+ | Only if manifest has `s3://` URLs |
| gcloud (Google Cloud SDK) | recent | Only if manifest has `gs://` URLs |

### Reference data

| Resource | Approx size | Source |
|----------|-------------|--------|
| GRCh38 primary assembly FASTA + `.fai` | 3.0 GB | Ensembl release 115 |
| Ensembl GFF3 (matched release) | 1.4 GB | Ensembl release 115 |
| dbNSFP v4.5+ TSV (bgzipped) | 35 GB | https://sites.google.com/site/jpopgen/dbNSFP |
| gnomAD v4.1 sites VCF (popmax-merged) | 80 GB | gnomAD v4.1 release |
| gnomAD v4.1 gene constraint TSV | 92 MB | https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv |
| ClinVar VCF + index | 80 MB | NCBI ClinVar weekly release |
| SpliceAI scores (optional, Stage 1 currently skips) | 40 GB | Illumina precomputed |

### Disk

A 200-sample WGS cohort run needs ~600 GB peak (per-sample fastVEP-annotated VCFs +
intermediate cohort tables + final outputs). Final outputs fit in ~50 GB.

---

## 2. Building fastVEP with dbNSFP patches

The pipeline depends on a patched fastVEP that emits AlphaMissense, ESM1b, and REVEL
alongside SIFT/PolyPhen in the same VEP pass. The patches and integration steps live in
`INTEGRATION.md` (and the standalone `src/dbnsfp_scores.rs` reference module). Follow
those instructions, then come back here.

The fork repo: `ajs3nj/fastVEP` (the Huang-lab/fastVEP base + our patches).

```bash
# Clone and build
git clone git@github.com:ajs3nj/fastVEP.git
cd fastVEP
cargo build --release

# Verify dbNSFP source emits the new fields
./target/release/fastvep sa-info /path/to/dbnsfp.osa | grep -E 'AlphaMissense|ESM1b|REVEL'
```

---

## 3. Repository layout

```
fastvep-dbnsfp-scores/
├── README.md                       # collaborator-facing: output schema + interpretation
├── INSTALL.md                      # this file: setup + how to run
├── INTEGRATION.md                  # fastVEP fork integration steps (Rust)
├── docs/
│   ├── pipeline_methods.md         # comprehensive methods document
│   ├── tiering.md                  # tiering-scheme design spec
│   ├── noncoding_v2_plan.md        # roadmap for v2 non-coding annotation
│   └── dbnsfp_score_columns.md     # dbNSFP column extraction reference
├── R/
│   ├── tier_variants.R             # main tiering logic (data.table)
│   ├── cohort_figures.R            # 8 clinical/presentation figures
│   └── nf_modifier_genes.txt       # curated NF modifier gene list
├── scripts/
│   ├── cohort_pipeline.sh          # end-to-end cohort orchestrator
│   ├── csq_to_wide_tab.py          # VEP CSQ + FV_* INFO -> wide tab converter
│   ├── validate_cohort.sh          # 9-section validation suite
│   ├── patch_constraint_into_cohort.sh   # legacy recovery: join gnomAD constraint
│   ├── cohort.manifest.example.tsv # manifest format example
│   └── ...
├── src/dbnsfp_scores.rs            # standalone Rust module (header-indexed, tested)
├── tools/extract_dbnsfp_scores.py  # Python validation oracle for the Rust path
└── tests/
    ├── sample_dbnsfp.tsv           # tiny synthetic dbNSFP fixture
    ├── example_annotated_tx.tsv    # synthetic merged input for tier_variants.R
    └── example_genotypes.tsv       # synthetic 4-sample genotype file
```

---

## 4. Cohort manifest format

### 4.1 v2 manifest (recommended, batched)

The v2 pipeline reads a TSV manifest that supports remote VCF URLs and per-sample
metadata. Header required.

```tsv
sample_id	batch_id	vcf_url	ancestry	severity	family_id	coverage	contamination	qc_pass
NG104C2PK8	1	s3://nf-cohort/raw/NG104C2PK8.vcf.gz	AFR	3.2	FAM-001	32.5	0.012	TRUE
NG1JSQFDDW	1	s3://nf-cohort/raw/NG1JSQFDDW.vcf.gz	NFE	4.1	FAM-002	31.8	0.008	TRUE
NG18TID2K1	2	s3://nf-cohort/raw/NG18TID2K1.vcf.gz	EAS	3.8	FAM-005	31.2	0.011	TRUE
...
```

Required columns:
- `sample_id` — unique per row
- `batch_id` — samples in the same batch are downloaded, annotated, and cleaned
  up as a unit. Recommended batch size: 50 samples. Optional: if the column is
  missing entirely, all samples are treated as batch 1.
- Either `vcf_url` (for `s3://`, `gs://`, `http://`, `https://`, `ftp://`) OR
  `vcf_path` (local filesystem, symlinked into place — no download)

Optional metadata columns (recommended):
- `ancestry` — self-reported major group (AFR / AMR / ASJ / EAS / FIN / NFE /
  SAS / OTH). Enables ancestry-stratified popmax + cohort_af if provided.
- `severity` — numeric phenotype-severity score (any scale)
- `family_id` — samples sharing a family_id are treated as related in downstream analyses
- `coverage`, `contamination`, `qc_pass` — sample-level QC metrics; used for
  the sample-level QC report in the provenance log

See `scripts/cohort.manifest.v2.example.tsv` for a full template.

### 4.2 v1 manifest (legacy, single-batch)

The original v1 format used only `sample_id` and `vcf_path`. Still supported by
`scripts/cohort_pipeline.sh` for single-batch runs; the v2 batched orchestrator
tolerates it too (all samples go into implicit batch 1, no download step).

### 4.3 Downloading VCFs from Synapse

If your cohort's VCFs live on [Synapse](https://www.synapse.org) (Sage Bionetworks'
data platform), put the Synapse entity ID directly in the `vcf_url` column.

**Manifest format:**

```tsv
sample_id	batch_id	vcf_url	ancestry	...
NG104C2PK8	1	syn12345678	AFR	...
NG1JSQFDDW	1	syn12345679	NFE	...
```

Both `syn12345678` and `syn://syn12345678` are accepted.

**Authentication.** The Synapse Python client (used under the hood by
`synapse get`) needs credentials before it can download. Two options, in
order of preference for automated pipelines:

**(A) Environment variable (preferred for automated runs)** — generate a
Personal Access Token at https://www.synapse.org/#!PersonalAccessTokens
with `Download` scope, then export it before running the pipeline:

```bash
export SYNAPSE_AUTH_TOKEN='<your token here>'
scripts/cohort_pipeline_batched.sh --manifest ... --out-dir ...
```

Store the token in a secrets manager, not in your shell history; add
`.env` files to `.gitignore` if you use them.

**(B) ~/.synapseConfig (interactive / dev machines)** — run
`synapse login -p <username>` once to write the config file. The download
loop will find it automatically. Not recommended for shared compute
because the token lives on disk in plaintext.

**Access control.** Confirm your Synapse account has been granted download
access to the study's parent project *before* starting the run. The v2
orchestrator only knows to fail when the first sample errors — a permission
issue caught 400 samples in is expensive.

**Rate limiting and retries.** `synapse get` is single-threaded per invocation
and Synapse rate-limits aggressive downloaders. The batched orchestrator
downloads samples serially within a batch by default. If a download fails
mid-batch, the state sentinel for that batch's `downloaded` stage is NOT
written, so a rerun retries the whole batch's download step (already-present
files skip cheaply). Persistent failures usually mean an auth problem or a
Synapse-side outage; check the sample's Synapse page in a browser.

**Index files.** Synapse projects sometimes store `.tbi` index files as
separate entities. The download step will pull the `.tbi` alongside the VCF
if `synapse get` returns both; if only the VCF is retrieved, Stage 0's
`tabix -p vcf -f` recreates the index during normalization.

### 4.4 VCF requirements

VCFs must be:
- bgzipped + tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`) — download step will
  index if the URL points at an unindexed file
- Aligned to GRCh38
- Per-sample (joint-called multi-sample VCFs also work; see notes in the
  pipeline script)

**Chromosome naming and indel encoding** are handled by Stage 0 (`bcftools
norm -m- -f FASTA` + optional `--chr-remap`), which runs inside every batch
before fastVEP annotation. Provide `scripts/chr_remap.example.tsv` (or a
custom map) if your input VCFs use UCSC-style chrom names (`chr1`) and your
downstream tools expect Ensembl-style (`1`).

---

## 5. Running the cohort pipeline

### 5.1 v2 batched invocation (recommended)

For a ~500-sample cohort where disk cannot hold all raw VCFs simultaneously:

```bash
scripts/cohort_pipeline_batched.sh \
    --manifest        manifest.v2.tsv \
    --out-dir         /path/to/results \
    --data-dir        /path/to/data \
    --sa-dir          /path/to/sa \
    --cache           /path/to/grch38_115.fvcache \
    --fasta           /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --chr-remap       scripts/chr_remap.example.tsv \
    --threads         16 \
    --modifier-genes  R/nf_modifier_genes.txt \
    --af-max          1e-4 \
    --filter-af       0.01 \
    --loftee --loftee-dir /path/to/loftee \
    --vep-cache /path/to/vep_cache
```

Batching works by looping over distinct `batch_id` values from the manifest.
For each batch:

1. **Download** VCFs into `$OUT_DIR/inputs/batch_<N>/`
2. **Stage 0 normalize** — `bcftools norm -m- -f FASTA` + chrom rename → `$OUT_DIR/normalized/batch_<N>/`
3. **Stage 1a annotate** — fastvep + csq_to_wide_tab.py → `$OUT_DIR/per_sample/<sample>.annotated.tab.gz`
4. **Stage 1b genotypes** — `bcftools query` → `$OUT_DIR/per_sample/<sample>.genotypes.tsv`
5. **(Optional) Stage 1c LOFTEE** — VEP+LOFTEE plugin pass; re-writes the wide tab with `loftee_lof` populated
6. **Stage 0z cleanup** — deletes `inputs/`, `normalized/`, `annotated/` for
   this batch after verifying every sample produced both a non-empty
   `annotated.tab.gz` and `genotypes.tsv`

After all batches finish, cohort stages 2–5 run once against the cumulative
`per_sample/` directory (Stage 2 merge, Stage 3 pre-filter, Stage 3.5 cohort
summary + genotypes, Stage 4 R tier, Stage 5 per-gene burden).

**Resumability:** killing the pipeline mid-run and restarting picks up where
it left off. State sentinels under `$OUT_DIR/state/` mark completed
stage-per-batch pairs. Force-retry a single batch with `--retry-batch <N>`;
process just one batch and stop with `--only-batch <N>`.

**Peak disk per batch** (50 samples, ~30 GB VCFs): ~1.5 TB. Cumulative
annotated tabs across all batches: ~100 GB for 500 samples.

### 5.2 v1 non-batched invocation (legacy)

For small cohorts (<100 samples) or when all VCFs already live on local disk:

```bash
scripts/cohort_pipeline.sh \
    --manifest        my.manifest.tsv \
    --out-dir         /path/to/results \
    --data-dir        /path/to/data \      # GFF3 + FASTA + transcript cache live here
    --sa-dir          /path/to/sa \        # built .osa/.osi/.oga files
    --cache           /path/to/grch38_115.fvcache \
    --fasta           /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --threads         16 \
    --modifier-genes  R/nf_modifier_genes.txt \
    --af-max          1e-4 \
    --filter-af       0.01
```

The pipeline runs in stages and is idempotent — each stage skips if its output already
exists. To re-run from a specific stage, use the `--skip-*` flags or delete the
corresponding output file.

| Stage | Output | Skip flag | Approx time (209-sample cohort) |
|-------|--------|-----------|---------------------------------|
| 1. Per-sample fastVEP annotate | `per_sample_vcf/*.vep.vcf.gz` + `per_sample_tab/*.tab.gz` | `--skip-annotate` | 30–90 min with 16 threads |
| 2. Cohort merge | `cohort.annotated.tab` | `--skip-merge` | 30 min |
| 3. Pre-filter | `cohort.filtered.tab` | (always runs) | 10 min |
| 3.5 Cohort summary | `cohort.filtered.with_cohort.tab`, `cohort.genotypes.tsv` | (always runs) | 30 min |
| 4. Tiering | `cohort.variants.tsv`, `cohort.genes.tsv` | (always runs) | 1–2 min |
| 5. Per-gene burden | per-gene burden columns merged into `cohort.genes.tsv` | (always runs) | 10 min |
| Validation | (stdout) | run separately | 5 min |
| Figures | `figures/*.png` | run separately | 10 min |

Each run writes a `_provenance.log` to the output directory recording the exact command,
git SHA, and reference file paths used.

---

## 6. Common operational variants

### Re-tier only (no fastVEP run)

```bash
Rscript R/tier_variants.R \
    --input          /path/to/cohort.filtered.with_cohort.tab \
    --modifier-genes R/nf_modifier_genes.txt \
    --af-max         1e-4 \
    --out-prefix     /path/to/cohort
```

Writes `<prefix>.variants.tsv` and `<prefix>.genes.tsv`. Useful when iterating on tier
rules or modifier-gene list without re-running annotation.

### Re-render figures only

```bash
Rscript R/cohort_figures.R \
    --variants    /path/to/cohort.variants.tsv \
    --genes       /path/to/cohort.genes.tsv \
    --genotypes   /path/to/cohort.genotypes.tsv \
    --out-dir     /path/to/figures \
    --cohort-name "NF1 batch1 (n=209)"
```

### Validate outputs

```bash
scripts/validate_cohort.sh --out-dir /path/to/results
```

Runs 9 sections (file presence, required columns, tier-distribution shape, class-conditional
rules, ClinVar overrides, modifier-candidate flag, anchor-gene signal, carrier-count
consistency, non-coding gap visibility). Exits non-zero on any FAIL.

### Patch legacy cohort tables for constraint columns

If you have a cohort table built before the `csq_to_wide_tab.py` SYMBOL-keyed fix (which is
the case for the v1 NF1 cohort), the `loeuf` and `pli` columns are silently empty and
`tier_variants.R` will not promote pLoF to Tier 1. Recover without re-running Stage 1:

```bash
scripts/patch_constraint_into_cohort.sh \
    /path/to/gnomad.v4.1.constraint_metrics.tsv \
    /path/to/cohort.filtered.with_cohort.tab
```

The script joins gnomAD v4.1 constraint metrics by MANE-select gene symbol and writes the
result back in place (with a `.preconstraintpatch.bak` of the original).

---

## 7. Getting the repo onto EC2

Plain git workflow. Push from your dev box, pull on EC2:

```bash
# On dev box
git remote add origin git@github.com:<you>/fastvep-dbnsfp-scores.git
git push -u origin main

# On EC2
git clone git@github.com:<you>/fastvep-dbnsfp-scores.git
```

---

## 8. Troubleshooting reference

Common failure modes the pipeline catches via upfront checks or sentinel messages:

| Symptom | Cause | Fix |
|---------|-------|-----|
| "Chromosome naming mismatch" upfront fail | VCFs are UCSC-style (`chr1`) but FASTA is Ensembl-style (`1`) | `bcftools annotate --rename-chrs` |
| Tier 1 is empty across whole cohort | `loeuf` / `pli` columns silently empty (pre-fix `csq_to_wide_tab.py`) | Run `scripts/patch_constraint_into_cohort.sh` |
| `clin_stars = 0` for every variant | Pre-fix `csq_to_wide_tab.py` REVIEW_STATUS parser | Already fixed in repo; for legacy data the R rule auto-relaxes (methods §7.10) |
| Stage 2c silently produces only header | `pipefail` + SIGPIPE from `... | head -1` | Pipeline now wraps the relevant calls in `( set +o pipefail; ... )` subshells |
| Stage 4 R script OOM-killed | Trying to fread the 32 GB genotype table | Stage 3.5 pre-aggregates carrier counts; R never reads genotypes directly |
| xargs/parallel fills `/tmp` and dies | Container `/tmp` on overlay is constrained | Pipeline exports `TMPDIR=$OUT_DIR/_parallel_tmp` |
| sample-gene hotspot figure dominated by ZNF/KRT/MUC | Short-read mis-mapping in paralog clusters | `R/cohort_figures.R` filters these families via `QC_ARTIFACT_GENE_PATTERNS` |

For anything not on this list, `docs/pipeline_methods.md` §7 (Known gotchas and limitations)
is the comprehensive reference.

---

## 9. Pointers for further work

- New scoring sources (SpliceAI restore, conservation, regulatory, LOFTEE): see
  `docs/noncoding_v2_plan.md` for a tiered implementation plan (Tier A–D priorities).
- Tier-rule modifications: edit `R/tier_variants.R` `CFG` block at the top, then re-tier
  per §6. The validate script's class-conditional checks will tell you if a rule change
  breaks the invariants.
- Adding a new dbNSFP score: see `INTEGRATION.md` for the fastVEP-side wiring plus
  `docs/dbnsfp_score_columns.md` for the column-extraction reference. Add the new column to
  `OUT_COLS` and the `FV_DBNSFP` extraction block in `csq_to_wide_tab.py`, then surface it
  in `tier_variants.R`.
