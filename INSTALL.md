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

The pipeline takes a TSV manifest of `(sample_id, vcf_path)` rows. Header required.

```tsv
sample_id	vcf_path
NG1F9ZXJ42	/data/nf1/inputs/batch1/NG1F9ZXJ42.vcf.gz
NG1CEFH8R3	/data/nf1/inputs/batch1/NG1CEFH8R3.vcf.gz
...
```

VCFs must be:
- bgzipped + tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`)
- Aligned to GRCh38 (Ensembl-style chromosome names: `1`, `2`, ..., `X`, `Y`, `MT`)
- Per-sample (joint-called multi-sample VCFs also work; see notes in `cohort_pipeline.sh`)

**Chromosome naming.** If your VCFs use UCSC-style names (`chr1`, `chrX`), they must
be renamed to Ensembl-style before annotation. The pipeline includes an upfront check
that fails fast on mismatch. Use:

```bash
bcftools annotate --rename-chrs chr_remap.tsv input.vcf.gz -Oz -o input.renamed.vcf.gz
```

`chr_remap.tsv` is a 2-column TSV: `chr1\t1`, `chr2\t2`, etc. (See
`scripts/chr_remap.example.tsv`.)

---

## 5. Running the cohort pipeline

End-to-end:

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
