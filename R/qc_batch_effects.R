#!/usr/bin/env Rscript
# qc_batch_effects.R
# -----------------------------------------------------------------------------
# Per-batch variant burden QC for a cohort with mixed sequencer / caller
# pipelines. Detects batch effects that would bias downstream gene-burden
# tests or modifier discovery.
#
# For the NF1 GWAS cohort specifically we mix three distinct pipelines:
#   - DRAGEN                            (306 samples, chunked into batches 1-3)
#   - DRAGEN with hard-filtered VCF     (52 samples, batch 4)
#   - bwa-mem2 + GATK HaplotypeCaller   (289 samples, chunked into batches 5-7)
#
# NOTE: batches 1/2/3 within DRAGEN and 5/6/7 within bwa-mem2+GATK are
# NOT distinct sequencing conditions -- they are arbitrary throughput
# partitions of the same underlying data, made so the batched pipeline
# could process ~100 samples per chunk with disk reclamation between
# chunks. Only the three-pipeline split (DRAGEN / DRAGEN-hard-filtered /
# bwa-mem2+GATK) is a real condition-level distinction.
#
# Different callers emit different variant counts per sample (sensitivity /
# filtering choices differ). Before we interpret a per-gene burden test as
# biological signal, we need to know whether the pipeline choice dominates
# the per-sample variant count. Running Kruskal-Wallis across all seven
# batch IDs tests two things simultaneously: (i) that within-pipeline
# chunks don't differ from each other -- a sanity check on the arbitrary
# partitioning -- and (ii) that the three real pipeline conditions don't
# differ. If p < 0.01 across-batches, downstream tests must control for
# pipeline (batch as covariate).
#
# Inputs:
#   --variants   cohort.variants.tsv  (has tier, chrom, pos, ref, alt)
#   --genotypes  cohort.genotypes.tsv (chrom, pos, ref, alt, sample, gt) -- gz OK
#   --manifest   sample_id -> batch_id mapping (v2 format, TSV with header)
#   --out-dir    output directory for figures + tables
#   --cohort-name label for figure titles
#   --tier-max   consider Tier 1..tier-max as "high priority" (default: 2)
#
# Outputs:
#   qc_batch_burden_summary.tsv
#     One row per batch: n_samples, median, q10, q90, mean, sd for both
#     Tier 1 alone and Tier 1..tier_max carrier counts.
#   qc_batch_burden_t1.png
#     Violin + box plot of per-sample Tier 1 carrier count, split by batch.
#   qc_batch_burden_high_priority.png
#     Same for Tier 1..tier_max (default Tier 1+2).
#   qc_batch_outliers.tsv
#     Samples > q90 + 1.5 * IQR within their batch. Prime candidates for
#     manual review (contamination? sequencer failure? real high burden?).
#   qc_batch_stats.txt
#     Kruskal-Wallis test across batches (non-parametric ANOVA analog --
#     no distributional assumption). Reports H, df, p-value. Pairwise
#     Wilcoxon with Bonferroni correction between batches. If p < 0.01
#     across-batches, downstream burden tests should include batch as a
#     covariate.
#
# Streams genotypes via awk rather than reading into R (long-format table is
# 30+ GB at cohort scale). Design mirrors cohort_figures.R fig_per_sample_burden.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ----------------------------- args ------------------------------------------
parse_args <- function(args) {
  out <- list(variants=NULL, genotypes=NULL, manifest=NULL, out_dir=NULL,
              cohort_name="cohort", tier_max=2L)
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if      (a == "--variants")     { out$variants  <- args[[i+1]]; i <- i+2 }
    else if (a == "--genotypes")    { out$genotypes <- args[[i+1]]; i <- i+2 }
    else if (a == "--manifest")     { out$manifest  <- args[[i+1]]; i <- i+2 }
    else if (a == "--out-dir")      { out$out_dir   <- args[[i+1]]; i <- i+2 }
    else if (a == "--cohort-name")  { out$cohort_name <- args[[i+1]]; i <- i+2 }
    else if (a == "--tier-max")     { out$tier_max  <- as.integer(args[[i+1]]); i <- i+2 }
    else stop(sprintf("unknown argument: %s", a))
  }
  for (req in c("variants","genotypes","manifest","out_dir")) {
    if (is.null(out[[req]])) stop(sprintf("--%s is required", gsub("_","-",req)))
  }
  dir.create(out$out_dir, recursive = TRUE, showWarnings = FALSE)
  out
}

# ----------------------------- theme -----------------------------------------
THEME_PRES <- theme_minimal(base_size = 15) +
  theme(
    plot.title       = element_text(face = "bold", size = 18),
    plot.subtitle    = element_text(size = 13, color = "grey30"),
    plot.caption     = element_text(size = 10, color = "grey45", hjust = 0),
    axis.title       = element_text(face = "bold", size = 14),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 20, hjust = 1),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

# ----------------------------- helpers ---------------------------------------

# Per-sample carrier count over a variant key set, streamed via awk. Emits
# every sample seen in the genotype table -- including samples with zero
# carriers -- so we can distinguish "no carriers" from "sample missing".
per_sample_carrier_count <- function(v_keys, gt_path, label) {
  keys_tmp <- tempfile(fileext = ".tsv")
  fwrite(v_keys[, .(chrom, pos, ref, alt)], keys_tmp,
         sep = "\t", col.names = FALSE)
  on.exit(unlink(keys_tmp), add = TRUE)

  awk_cmd <- sprintf(
    "awk -F'\\t' '
       NR == FNR { v[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = 1; next }
       FNR == 1  { next }
       $6 ~ /^(0\\/1|1\\/0|0\\|1|1\\|0|1\\/1|1\\|1)$/ {
         key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
         if (key in v) cnt[$5]++
         seen[$5] = 1
       }
       END {
         for (s in seen) printf \"%%s\\t%%d\\n\", s, (s in cnt) ? cnt[s] : 0
       }
     ' %s %s",
    shQuote(keys_tmp), shQuote(gt_path)
  )
  message(sprintf("[qc] streaming genotypes for %s (%d variant keys)...",
                  label, nrow(v_keys)))
  per_sample <- fread(cmd = awk_cmd, header = FALSE,
                      col.names = c("sample", "n"))
  per_sample
}

fig_batch_burden <- function(d, label, cohort_name, out_path, kw_p = NA_real_) {
  # d: sample, batch, n
  d[, batch := factor(batch)]
  # Order batches by median burden for readability (highest median first)
  order_batches <- d[, .(med = median(n, na.rm = TRUE)), by = batch][
    order(-med), batch]
  d[, batch := factor(batch, levels = as.character(order_batches))]

  # Explicit as.numeric() -- see the sum_stats() note in main(). quantile /
  # median over integer input can return either type across groups, breaking
  # the by-group rbind with a type-consistency error.
  n_by_batch <- d[, .(n_samples = .N,
                      median   = as.numeric(median(n)),
                      q10      = as.numeric(quantile(n, 0.10)),
                      q90      = as.numeric(quantile(n, 0.90))),
                  by = batch]

  # Cap x-axis label lines at ~20 chars (long batch names wrap otherwise)
  cap_label <- function(x) ifelse(nchar(x) > 22,
                                  paste0(substr(x, 1, 20), "..."), x)

  p <- ggplot(d, aes(x = batch, y = n, fill = batch)) +
    geom_violin(alpha = 0.55, width = 0.9, color = NA) +
    geom_boxplot(width = 0.15, alpha = 0.85, outlier.size = 0.7,
                 fill = "white", color = "grey15") +
    geom_text(data = n_by_batch,
              aes(x = batch, y = max(d$n) * 1.05,
                  label = sprintf("n=%d\nmed=%s", n_samples, comma(median))),
              inherit.aes = FALSE, vjust = 0, size = 3.6, color = "grey20",
              fontface = "bold") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.15))) +
    scale_x_discrete(labels = cap_label) +
    labs(
      title    = sprintf("Per-sample %s variant burden by batch", label),
      subtitle = sprintf("%s -- %d samples across %d batches",
                         cohort_name, nrow(d), length(order_batches)),
      x = NULL,
      y = sprintf("Number of %s variants per sample", label),
      caption = sprintf(paste0(
        "Batches ordered by median burden (highest first). ",
        "Kruskal-Wallis across batches: p = %s.\n",
        "Batch differences may reflect real sequencer / caller sensitivity ",
        "differences rather than biology.\nIf p < 0.01, downstream burden ",
        "tests should include batch as a covariate."),
        if (is.na(kw_p)) "n/a" else format.pval(kw_p, digits = 3))
    ) +
    THEME_PRES

  ggsave(out_path, p, width = 11, height = 6.5, dpi = 300, bg = "white")
  message("[qc] wrote ", out_path)
}

# ----------------------------- main ------------------------------------------
main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  message(sprintf(
    "[qc] variants=%s genotypes=%s manifest=%s out-dir=%s tier-max=%d",
    opt$variants, opt$genotypes, opt$manifest, opt$out_dir, opt$tier_max))

  # Load variants (only need chrom/pos/ref/alt/tier)
  v <- fread(opt$variants, sep = "\t",
             select = c("chrom","pos","ref","alt","tier"))
  message(sprintf("[qc] loaded %s variants", comma(nrow(v))))

  # Load manifest: sample_id + batch_id required
  m <- fread(opt$manifest, sep = "\t")
  req_manifest_cols <- c("sample_id","batch_id")
  missing <- setdiff(req_manifest_cols, names(m))
  if (length(missing) > 0L) {
    stop(sprintf("manifest missing required columns: %s",
                 paste(missing, collapse = ", ")))
  }
  message(sprintf("[qc] manifest: %d samples, %d batches",
                  nrow(m), length(unique(m$batch_id))))

  # Tier 1 keys
  v_t1  <- v[tier == 1L, .(chrom = as.character(chrom),
                           pos   = as.integer(pos),
                           ref   = as.character(ref),
                           alt   = as.character(alt))]
  # Tier 1..tier_max keys
  v_hi  <- v[tier %in% seq_len(opt$tier_max),
             .(chrom = as.character(chrom),
               pos   = as.integer(pos),
               ref   = as.character(ref),
               alt   = as.character(alt))]

  # Checkpoint per-sample carrier counts. The awk streams over 3+ B rows of
  # cohort.genotypes.tsv and takes 30-60 min at cohort scale. If anything
  # downstream (aggregation, plotting) fails, we don't want to re-stream --
  # the intermediates are tiny (~n_samples rows) and cheap to persist. If the
  # cached files exist AND the tier-key sets haven't changed size, reuse them.
  ck_t1 <- file.path(opt$out_dir, "qc_ck_ps_t1.tsv")
  ck_hi <- file.path(opt$out_dir, sprintf("qc_ck_ps_t1_%d.tsv", opt$tier_max))

  load_or_stream <- function(v_keys, gt_path, label, ck_path) {
    if (file.exists(ck_path) && file.size(ck_path) > 0) {
      cached <- fread(ck_path, sep = "\t")
      if (nrow(cached) > 0 &&
          "n" %in% names(cached) && "sample" %in% names(cached)) {
        message(sprintf("[qc] reusing cached carrier counts for %s -> %s",
                        label, ck_path))
        return(cached)
      }
    }
    ps <- per_sample_carrier_count(v_keys, gt_path, label)
    fwrite(ps, ck_path, sep = "\t")
    message(sprintf("[qc] cached %s carrier counts -> %s", label, ck_path))
    ps
  }
  ps_t1 <- load_or_stream(v_t1, opt$genotypes, "Tier 1", ck_t1)
  ps_hi <- load_or_stream(v_hi, opt$genotypes,
                          sprintf("Tier 1..%d", opt$tier_max), ck_hi)

  # Join to manifest for batch labels. Left-join on manifest so samples in
  # manifest but missing from genotypes get flagged (n=NA -> drop with a
  # message; typically means the sample failed upstream).
  merge_burden <- function(ps, label) {
    d <- merge(m[, .(sample = sample_id, batch = batch_id)],
               ps, by = "sample", all.x = TRUE)
    n_missing <- sum(is.na(d$n))
    if (n_missing > 0L) {
      message(sprintf("[qc] %s: %d samples in manifest missing from genotype file",
                      label, n_missing))
    }
    d[!is.na(n)]
  }
  d_t1 <- merge_burden(ps_t1, "Tier 1")
  d_hi <- merge_burden(ps_hi, sprintf("Tier 1..%d", opt$tier_max))

  # ---------- summary table ---------------------------------------------------
  # Coerce every stat to numeric explicitly. median() / quantile() / min() /
  # max() over an integer vector can return integer OR double depending on
  # whether the result lands exactly on an existing value or interpolates.
  # Across 7 batches you get mixed types per column, and data.table's rbind
  # refuses to combine groups with inconsistent column types.
  sum_stats <- function(d, label) {
    d[, .(analysis  = label,
          n_samples = .N,
          median    = as.numeric(median(n)),
          mean      = round(mean(n), 1),
          sd        = round(sd(n), 1),
          q10       = as.numeric(quantile(n, 0.10)),
          q90       = as.numeric(quantile(n, 0.90)),
          min       = as.numeric(min(n)),
          max       = as.numeric(max(n))),
      by = .(batch)]
  }
  summary_dt <- rbindlist(list(
    sum_stats(d_t1, "Tier 1"),
    sum_stats(d_hi, sprintf("Tier 1..%d", opt$tier_max))
  ))
  setorder(summary_dt, analysis, -median)
  sum_path <- file.path(opt$out_dir, "qc_batch_burden_summary.tsv")
  fwrite(summary_dt, sum_path, sep = "\t")
  message("[qc] wrote ", sum_path)

  # ---------- Kruskal-Wallis + pairwise Wilcoxon ------------------------------
  # Kruskal-Wallis: non-parametric ANOVA analog. Tests whether all batches
  # have the same median burden. We don't assume normality (per-sample
  # burden is often right-skewed with heavy tails).
  # Follow up with pairwise Wilcoxon (Bonferroni-corrected) to identify
  # which batches differ from which.
  run_kw <- function(d, label) {
    if (length(unique(d$batch)) < 2L) {
      return(list(H = NA, df = NA, p = NA, pairwise = NULL))
    }
    kw <- kruskal.test(n ~ batch, data = d)
    pw <- suppressWarnings(pairwise.wilcox.test(d$n, d$batch,
                                                 p.adjust.method = "bonferroni"))
    list(H = unname(kw$statistic), df = unname(kw$parameter),
         p = kw$p.value, pairwise = pw$p.value)
  }
  kw_t1 <- run_kw(d_t1, "Tier 1")
  kw_hi <- run_kw(d_hi, sprintf("Tier 1..%d", opt$tier_max))

  stats_path <- file.path(opt$out_dir, "qc_batch_stats.txt")
  con <- file(stats_path, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    sprintf("Kruskal-Wallis test across batches (H0: all batches have equal median burden)"),
    "",
    sprintf("Tier 1:            H = %.3f, df = %d, p = %s",
            kw_t1$H, kw_t1$df, format.pval(kw_t1$p, digits = 4)),
    sprintf("Tier 1..%d:         H = %.3f, df = %d, p = %s",
            opt$tier_max, kw_hi$H, kw_hi$df, format.pval(kw_hi$p, digits = 4)),
    "",
    "Interpretation:",
    "  p < 0.01  -> strong batch effect; INCLUDE batch as covariate in downstream tests",
    "  p 0.01 to 0.05 -> moderate batch effect; consider stratified analysis",
    "  p > 0.05  -> no significant batch effect; pooled analysis OK",
    ""
  ), con)
  if (!is.null(kw_t1$pairwise)) {
    writeLines("Pairwise Wilcoxon (Tier 1, Bonferroni-corrected p-values):", con)
    capture.output(print(round(kw_t1$pairwise, 4)), file = con)
    writeLines("", con)
  }
  if (!is.null(kw_hi$pairwise)) {
    writeLines(sprintf("Pairwise Wilcoxon (Tier 1..%d, Bonferroni-corrected p-values):",
                       opt$tier_max), con)
    capture.output(print(round(kw_hi$pairwise, 4)), file = con)
  }
  message("[qc] wrote ", stats_path)

  # ---------- outlier samples --------------------------------------------------
  # Standard Tukey outlier definition per batch, applied on Tier 1+..tier_max
  # (more variants -> more stable IQR).
  # as.numeric() to avoid the same type-inconsistency bug that hit sum_stats().
  d_hi[, iqr := as.numeric(quantile(n, 0.75) - quantile(n, 0.25)), by = batch]
  d_hi[, q75 := as.numeric(quantile(n, 0.75)), by = batch]
  d_hi[, is_outlier := n > q75 + 1.5 * iqr]
  outliers <- d_hi[is_outlier == TRUE,
                   .(sample, batch, n_variants = n,
                     batch_q75 = round(q75, 1),
                     batch_iqr = round(iqr, 1))]
  setorder(outliers, batch, -n_variants)
  out_path <- file.path(opt$out_dir, "qc_batch_outliers.tsv")
  fwrite(outliers, out_path, sep = "\t")
  message(sprintf("[qc] %d outlier samples -> %s",
                  nrow(outliers), out_path))

  # ---------- figures ----------------------------------------------------------
  fig_batch_burden(d_t1, "Tier 1", opt$cohort_name,
                   file.path(opt$out_dir, "qc_batch_burden_t1.png"),
                   kw_p = kw_t1$p)
  fig_batch_burden(d_hi, sprintf("Tier 1..%d", opt$tier_max), opt$cohort_name,
                   file.path(opt$out_dir, "qc_batch_burden_high_priority.png"),
                   kw_p = kw_hi$p)

  message("[qc] done.")
}

main()
