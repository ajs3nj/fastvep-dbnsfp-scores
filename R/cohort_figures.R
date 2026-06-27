#!/usr/bin/env Rscript
# cohort_figures.R
# -----------------------------------------------------------------------------
# Clinical / presentation-quality summary figures from a cohort tiering run.
#
# Inputs:
#   --variants   cohort.variants.tsv  (one row per unique tiered variant)
#   --genes      cohort.genes.tsv     (one row per gene with rollup + burden)
#   --genotypes  cohort.genotypes.tsv (optional; needed for per-sample burden)
#   --out-dir    output directory for PNGs (default: <variants dir>/figures)
#   --cohort-name label for figure titles (default: "cohort")
#
# Outputs (PNG, 300 dpi, 10x6 in unless noted):
#   01_tier_distribution.png        Cohort-wide variants by tier (1-5)
#   02_variant_class_by_tier.png    Stacked bar: variant_class composition per tier
#   03_per_sample_burden.png        Histogram of Tier 1+2 variants per sample
#   04_modifier_landscape.png       Top 20 genes by n_modifier_candidates
#   05_nf_anchor_landscape.png      NF1/NF2/SMARCB1/LZTR1 variants by tier (positive control)
#   06_variant_recurrence.png       n_carriers distribution -- private vs recurrent
#   07_top_recurrent_t12.png        Top 30 most-recurrent Tier 1+2 variants in cohort
#   08_sample_gene_hotspots.png     (sample, gene) pairs with >=2 rare Tier 1-3 variants
#                                   (compound-het / multi-hit candidates)
#
# Dependencies: data.table, ggplot2, scales. (All available via conda-forge.)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ----------------------------- args ------------------------------------------
parse_args <- function(args) {
  out <- list(variants=NULL, genes=NULL, genotypes=NULL, out_dir=NULL,
              cohort_name="cohort")
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if      (a == "--variants")    { out$variants  <- args[[i+1]]; i <- i+2 }
    else if (a == "--genes")       { out$genes     <- args[[i+1]]; i <- i+2 }
    else if (a == "--genotypes")   { out$genotypes <- args[[i+1]]; i <- i+2 }
    else if (a == "--out-dir")     { out$out_dir   <- args[[i+1]]; i <- i+2 }
    else if (a == "--cohort-name") { out$cohort_name <- args[[i+1]]; i <- i+2 }
    else stop(sprintf("unknown argument: %s", a))
  }
  if (is.null(out$variants)) stop("--variants is required")
  if (is.null(out$genes))    stop("--genes is required")
  if (is.null(out$out_dir))  out$out_dir <- file.path(dirname(out$variants), "figures")
  dir.create(out$out_dir, recursive = TRUE, showWarnings = FALSE)
  out
}

# ---------------------------- theme ------------------------------------------
# Slides-ready: clean, light grid, large labels, presentation typography.
THEME_PRES <- theme_minimal(base_size = 16) +
  theme(
    plot.title         = element_text(face = "bold", size = 20, margin = margin(b = 6)),
    plot.subtitle      = element_text(size = 14, color = "grey30", margin = margin(b = 14)),
    plot.caption       = element_text(size = 11, color = "grey45", hjust = 0, margin = margin(t = 12)),
    axis.title         = element_text(face = "bold", size = 15),
    axis.text          = element_text(size = 13),
    legend.title       = element_text(face = "bold", size = 13),
    legend.text        = element_text(size = 12),
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold", size = 14)
  )

# Tier color palette: red (most pathogenic) -> green (benign), with a neutral middle.
TIER_COLORS <- c("1" = "#c0392b",   # red
                 "2" = "#e67e22",   # orange
                 "3" = "#f1c40f",   # yellow
                 "4" = "#95a5a6",   # grey (uncertain)
                 "5" = "#27ae60")   # green
TIER_LABELS <- c("1" = "Tier 1\n(most likely\npathogenic)",
                 "2" = "Tier 2\n(likely\npathogenic)",
                 "3" = "Tier 3\n(uncertain,\nleaning damaging)",
                 "4" = "Tier 4\n(uncertain,\nno info)",
                 "5" = "Tier 5\n(very likely\nbenign)")

# Variant-class palette: ordered by severity, distinct categorical hues.
CLASS_ORDER  <- c("pLoF", "splice_noncanonical", "missense", "inframe",
                  "noncoding", "regulatory", "other")
CLASS_LABELS <- c(pLoF                  = "pLoF",
                  splice_noncanonical   = "Splice (non-canonical)",
                  missense              = "Missense",
                  inframe               = "In-frame indel",
                  noncoding             = "Non-coding (intronic/UTR)",
                  regulatory            = "Regulatory",
                  other                 = "Other")
CLASS_COLORS <- c(pLoF                  = "#8e44ad",   # purple
                  splice_noncanonical   = "#2980b9",   # blue
                  missense              = "#16a085",   # teal
                  inframe               = "#f39c12",   # amber
                  noncoding             = "#7f8c8d",   # grey
                  regulatory            = "#d35400",   # rust
                  other                 = "#bdc3c7")   # light grey

# ----------------------------- figures ---------------------------------------

fig_tier_distribution <- function(v, cohort_name, out_path) {
  d <- v[!is.na(tier), .N, by = .(tier)]
  d[, tier := factor(tier, levels = c(1L,2L,3L,4L,5L))]
  setorder(d, tier)
  d[, pct := N / sum(N) * 100]

  p <- ggplot(d, aes(x = tier, y = N, fill = as.character(tier))) +
    geom_col(width = 0.75, color = "white") +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", comma(N), pct)),
              vjust = -0.4, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = TIER_COLORS, labels = TIER_LABELS, name = "Tier") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    labs(
      title    = "Variant tier distribution",
      subtitle = sprintf("%s — %s tiered variants across %s samples",
                         cohort_name,
                         comma(nrow(v)),
                         if ("n_carriers" %in% names(v))
                           comma(max(v$n_carriers, na.rm = TRUE)) else "the cohort"),
      x = NULL,
      y = "Number of variants",
      caption = "Tier reflects pathogenicity likelihood. Tier 1 = most likely pathogenic; Tier 5 = very likely benign.\nClass-conditional logic: per-class evidence rules (see docs/tiering.md)."
    ) +
    THEME_PRES + theme(legend.position = "none")

  ggsave(out_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

fig_variant_class_by_tier <- function(v, cohort_name, out_path) {
  cls_col <- if ("variant_class" %in% names(v)) "variant_class" else NA
  if (is.na(cls_col)) {
    message("[fig] skipping variant_class_by_tier (no variant_class column in input)")
    return(invisible())
  }
  d <- v[!is.na(tier) & !is.na(get(cls_col)),
         .N, by = .(tier, variant_class = get(cls_col))]
  d[, tier := factor(tier, levels = c(1L,2L,3L,4L,5L))]
  d[, variant_class := factor(variant_class,
                              levels = intersect(CLASS_ORDER, unique(variant_class)))]
  setorder(d, tier)

  p <- ggplot(d, aes(x = tier, y = N, fill = variant_class)) +
    geom_col(position = "fill", color = "white", width = 0.75) +
    scale_fill_manual(values = CLASS_COLORS, labels = CLASS_LABELS,
                      name = "Variant class", drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(labels = sub("\n.*", "", TIER_LABELS)) +
    labs(
      title    = "Variant class composition by tier",
      subtitle = sprintf("%s — what kinds of variants land in each tier", cohort_name),
      x = "Tier",
      y = "Proportion of variants",
      caption = "pLoF and missense dominate higher tiers; non-coding accumulates in Tier 4 where SpliceAI / conservation signal is absent."
    ) +
    THEME_PRES

  ggsave(out_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

fig_per_sample_burden <- function(v, gt_path, cohort_name, out_path) {
  if (is.null(gt_path)) {
    message("[fig] skipping per_sample_burden (no --genotypes path provided)")
    return(invisible())
  }

  # Stream the genotype file via awk rather than fread()ing it -- at cohort scale
  # the long-format genotype TSV is 30+ GB and OOM-kills R. Awk holds only:
  #   - the Tier 1+2 variant key set (~thousands of keys, ~MB)
  #   - per-sample carrier counts (a few hundred samples, tiny)
  v_t12 <- v[tier %in% c(1L, 2L),
             .(chrom = as.character(get("chrom")),
               pos   = as.integer(get("pos")),
               ref   = as.character(get("ref")),
               alt   = as.character(get("alt")))]
  if (nrow(v_t12) == 0L) {
    message("[fig] skipping per_sample_burden (no Tier 1+2 variants)")
    return(invisible())
  }
  keys_tmp <- tempfile(fileext = ".tsv")
  fwrite(v_t12, keys_tmp, sep = "\t", col.names = FALSE)
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
  per_sample <- fread(cmd = awk_cmd, header = FALSE,
                      col.names = c("sample", "n_t12"))
  if (nrow(per_sample) == 0L) {
    message("[fig] skipping per_sample_burden (no carriers found in genotype file)")
    return(invisible())
  }

  med <- median(per_sample$n_t12)
  q90 <- quantile(per_sample$n_t12, 0.9)

  p <- ggplot(per_sample, aes(x = n_t12)) +
    geom_histogram(binwidth = 1, fill = TIER_COLORS["2"], color = "white", boundary = 0) +
    geom_vline(xintercept = med, linetype = "dashed", color = "grey30", linewidth = 0.6) +
    annotate("text", x = med, y = Inf,
             label = sprintf("median = %g", med),
             vjust = 1.6, hjust = -0.1, size = 4.5, color = "grey20") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(breaks = pretty_breaks(n = 8)) +
    labs(
      title    = "Per-sample burden: Tier 1 + Tier 2 variants",
      subtitle = sprintf("%s — %d samples, distribution of high-priority variants per sample",
                         cohort_name, nrow(per_sample)),
      x = "Number of Tier 1 or Tier 2 variants",
      y = "Number of samples",
      caption = sprintf("Median %g  |  90th percentile %g.\nOutlier samples may reflect contamination, ancestry effects, or genuine high burden — worth manual review.",
                        med, q90)
    ) +
    THEME_PRES

  ggsave(out_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

fig_modifier_landscape <- function(g, cohort_name, out_path, top_n = 20,
                                   modifier_genes = NULL) {
  if (!"n_modifier_candidates" %in% names(g)) {
    message("[fig] skipping modifier_landscape (no n_modifier_candidates column)")
    return(invisible())
  }
  d <- g[n_modifier_candidates > 0,
         .(gene, n_modifier_candidates,
           best_tier = if ("best_tier" %in% names(g)) best_tier else NA_integer_)]
  setorder(d, -n_modifier_candidates)
  d <- head(d, top_n)
  d[, gene := factor(gene, levels = rev(d$gene))]

  # Highlight known modifier genes if a list is provided
  if (!is.null(modifier_genes)) {
    mg <- toupper(modifier_genes)
    d[, is_anchor := toupper(as.character(gene)) %in% mg]
  } else {
    # Default highlight: NF / RAS-MAPK anchor set
    d[, is_anchor := toupper(as.character(gene)) %in% toupper(c(
      "NF1","NF2","SMARCB1","LZTR1",
      "KRAS","NRAS","HRAS","BRAF","RAF1","ARAF",
      "MAP2K1","MAP2K2","MAPK1","MAPK3","PTPN11","SPRED1","SPRED2","RASA1",
      "CDKN2A","CDKN2B","ATM","TP53","EED","SUZ12"))]
  }

  p <- ggplot(d, aes(x = n_modifier_candidates, y = gene, fill = is_anchor)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = comma(n_modifier_candidates)),
              hjust = -0.15, size = 4.3, color = "grey20") +
    scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#7f8c8d"),
                      labels = c(`TRUE` = "NF / RAS-MAPK anchor or modifier gene",
                                 `FALSE` = "Other"),
                      name = NULL) +
    scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    labs(
      title    = "Top genes by modifier candidate count",
      subtitle = sprintf("%s — genes with the most variants flagged as modifier candidates", cohort_name),
      x = "Number of modifier-candidate variants",
      y = NULL,
      caption = "Modifier candidates: rare (1e-4 < popmax <= 5e-2) or sub-Mendelian with class-appropriate moderate signal.\nHighlighted genes are NF / RAS-MAPK pathway anchors (see R/nf_modifier_genes.txt)."
    ) +
    THEME_PRES + theme(legend.position = "top")

  ggsave(out_path, p, width = 10, height = 8, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

# ----------------------------------------------------------------------------
# Figure 05 -- NF-spectrum tumor suppressor variant landscape (positive control)
#
# A 209-sample NF1 cohort should have substantial Tier 1-2 signal in the four
# anchor genes (NF1, NF2, SMARCB1, LZTR1). If this figure is empty or
# Tier-1-poor for any anchor, something is wrong upstream (constraint columns
# missing, ClinVar parsing broken, NF anchor override not firing, etc.).
# Treats the figure as a confidence check for the whole pipeline.
# ----------------------------------------------------------------------------
fig_nf_anchor_landscape <- function(v, cohort_name, out_path) {
  NF_ANCHORS <- c("NF1", "NF2", "SMARCB1", "LZTR1")
  d <- v[gene %in% NF_ANCHORS & !is.na(tier),
         .N, by = .(gene, tier)]
  if (nrow(d) == 0L) {
    message("[fig] skipping nf_anchor_landscape (no anchor-gene variants)")
    return(invisible())
  }
  d[, gene := factor(gene, levels = NF_ANCHORS)]
  d[, tier := factor(tier, levels = c(1L, 2L, 3L, 4L, 5L))]

  totals <- v[gene %in% NF_ANCHORS, .(total = .N), by = .(gene)]
  totals[, gene := factor(gene, levels = NF_ANCHORS)]
  t12 <- v[gene %in% NF_ANCHORS & tier %in% c(1L, 2L),
           .(t12_count = .N), by = .(gene)]
  t12[, gene := factor(gene, levels = NF_ANCHORS)]

  p <- ggplot(d, aes(y = gene, x = N, fill = tier)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    geom_text(data = totals, aes(y = gene, x = total, label = paste0("n=", total)),
              inherit.aes = FALSE, hjust = -0.15, size = 4.2, fontface = "bold",
              color = "grey20") +
    scale_fill_manual(values = TIER_COLORS, labels = TIER_LABELS,
                      drop = FALSE, name = "Tier") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)), labels = comma) +
    scale_y_discrete(limits = rev(NF_ANCHORS)) +
    labs(
      title    = "Positive control: NF-spectrum tumor suppressor variant landscape",
      subtitle = sprintf("%s -- variants by tier in the four NF anchor genes (NF1, NF2, SMARCB1, LZTR1)",
                         cohort_name),
      x = "Number of distinct variants in cohort",
      y = NULL,
      caption = paste(
        sprintf("Tier 1+2 totals -- %s",
                paste(sprintf("%s: %d", t12$gene, t12$t12_count), collapse = "  |  ")),
        "Tier 1 hits in these genes are immediate candidates for re-contact / phenotype review.",
        "If any anchor is Tier-1-empty, double-check that constraint columns and the NF",
        "anchor override (see tier_variants.R) are populated for this run.",
        sep = "\n"
      )
    ) +
    THEME_PRES

  ggsave(out_path, p, width = 11, height = 6.5, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

# ----------------------------------------------------------------------------
# Figure 06 -- Variant recurrence across samples
#
# Most rare variants in a small cohort are private (carriers = 1). Recurrent
# variants are either (a) real founder/recurrent mutations, (b) sequencing
# artifacts in low-complexity regions, or (c) common-but-mis-annotated
# variants slipping through the AF filter. This figure shows the recurrence
# distribution split by tier so collaborators can see what fraction of
# Tier 1-2 variants are recurrent (worth manual review) vs private.
# ----------------------------------------------------------------------------
fig_variant_recurrence <- function(v, cohort_name, out_path) {
  if (!"n_carriers" %in% names(v)) {
    message("[fig] skipping variant_recurrence (no n_carriers column)")
    return(invisible())
  }
  d <- v[!is.na(n_carriers) & n_carriers > 0]
  d[, tier_cat := fifelse(tier %in% c(1L, 2L), "Tier 1-2 (high priority)",
                  fifelse(tier %in% c(3L, 4L), "Tier 3-4 (uncertain)",
                                                "Tier 5 (likely benign)"))]
  d[, tier_cat := factor(tier_cat, levels = c("Tier 1-2 (high priority)",
                                              "Tier 3-4 (uncertain)",
                                              "Tier 5 (likely benign)"))]
  agg <- d[, .(n_variants = .N), by = .(n_carriers, tier_cat)]

  # Stats for caption
  t12_private <- d[tier_cat == "Tier 1-2 (high priority)" & n_carriers == 1, .N]
  t12_total   <- d[tier_cat == "Tier 1-2 (high priority)", .N]
  t12_recur   <- t12_total - t12_private
  max_carrier <- max(d$n_carriers)

  cat_colors <- c("Tier 1-2 (high priority)" = "#c0392b",
                  "Tier 3-4 (uncertain)"     = "#f1c40f",
                  "Tier 5 (likely benign)"   = "#27ae60")

  p <- ggplot(agg, aes(x = n_carriers, y = n_variants, color = tier_cat)) +
    geom_step(linewidth = 1.1, direction = "mid") +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_x_log10(labels = comma, breaks = c(1, 2, 5, 10, 20, 50, 100, 200)) +
    scale_y_log10(labels = comma) +
    scale_color_manual(values = cat_colors, name = "Tier category") +
    labs(
      title    = "Variant recurrence across the cohort",
      subtitle = sprintf("%s -- distribution of variants by number of carrier samples",
                         cohort_name),
      x = "Number of carrier samples (log)",
      y = "Number of variants (log)",
      caption = sprintf(paste(
        "Tier 1-2: %s of %s variants are private (carriers = 1); %s are recurrent.",
        "Recurrent high-confidence variants warrant inspection for founder effects",
        "vs sequencing artifacts (low-complexity regions, mappability gaps).",
        "Max carrier count in cohort: %d samples.",
        sep = "\n"),
        comma(t12_private), comma(t12_total), comma(t12_recur), max_carrier
      )
    ) +
    THEME_PRES +
    theme(legend.position = "top")

  ggsave(out_path, p, width = 11, height = 6.5, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

# ----------------------------------------------------------------------------
# Figure 07 -- Top recurrent Tier 1+2 variants
#
# Companion to Fig 06: the actual list of which recurrent Tier 1-2 variants
# are most carriers-rich. Use this to spot:
#   - Founder mutations (one variant, many carriers, biologically plausible gene)
#   - QC artifacts (variants in well-known mappability-poor genes; common
#     ClinVar entries that escaped the AF filter)
#   - Population-stratified variants (ancestry-specific common alleles
#     mis-annotated as rare due to missing-popmax)
# ----------------------------------------------------------------------------
fig_top_recurrent_variants <- function(v, cohort_name, out_path, top_n = 30) {
  if (!"n_carriers" %in% names(v)) {
    message("[fig] skipping top_recurrent_variants (no n_carriers column)")
    return(invisible())
  }
  d <- v[tier %in% c(1L, 2L) & !is.na(n_carriers) & n_carriers >= 2]
  if (nrow(d) == 0L) {
    message("[fig] skipping top_recurrent_variants (no recurrent Tier 1-2 variants)")
    return(invisible())
  }
  setorder(d, -n_carriers, gene)
  d <- head(d, top_n)
  d[, label := sprintf("%-8s  %s:%s  %s>%s",
                       ifelse(is.na(gene) | gene == "", ".", gene),
                       chrom, pos, ref, alt)]
  d[, label := factor(label, levels = rev(unique(label)))]
  d[, tier := factor(tier, levels = c(1L, 2L, 3L, 4L, 5L))]

  p <- ggplot(d, aes(x = n_carriers, y = label, fill = tier)) +
    geom_col() +
    geom_text(aes(label = paste0("n=", n_carriers)),
              hjust = -0.15, size = 3.6, color = "grey20") +
    scale_fill_manual(values = TIER_COLORS, labels = TIER_LABELS,
                      drop = FALSE, name = "Tier") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
    labs(
      title    = sprintf("Top %d recurrent Tier 1-2 variants in cohort", nrow(d)),
      subtitle = sprintf("%s -- candidates for founder effects, recurrent hits, or QC review",
                         cohort_name),
      x = "Number of carrier samples",
      y = NULL,
      caption = paste(
        "Variants present in two or more samples at Tier 1 or Tier 2 priority.",
        "Inspect top hits for: real founder mutations (concentrated by ancestry?),",
        "QC artifacts (low-complexity / repetitive regions, mappability gaps),",
        "or ancestry-private common alleles that slipped the gnomAD popmax filter.",
        sep = "\n"
      )
    ) +
    THEME_PRES +
    theme(axis.text.y = element_text(size = 10, family = "mono"))

  ggsave(out_path, p, width = 12, height = 0.32 * nrow(d) + 2.5, dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

# ----------------------------------------------------------------------------
# Figure 08 -- Per-sample, per-gene hotspots (compound-het / multi-hit)
#
# For each (sample, gene) pair, count how many distinct rare Tier 1-3 variants
# the sample carries in that gene. Pairs with >=2 distinct variants are
# candidates for:
#   - Compound heterozygous recessive effects (if both variants het)
#   - Multi-hit somatic events (e.g. two-hit tumor suppressor inactivation)
#   - Phase artifacts (same haplotype, double-counted)
# Phasing requires trio data or long-read sequencing for confirmation, but
# these pairs are the right manual-review starting point.
#
# Streams the 30+ GB genotype table via awk; never reads it into R.
# ----------------------------------------------------------------------------
fig_sample_gene_hotspots <- function(v, gt_path, cohort_name, out_path, top_n = 30) {
  if (is.null(gt_path)) {
    message("[fig] skipping sample_gene_hotspots (no --genotypes path provided)")
    return(invisible())
  }

  # Lookup table: variant_key -> gene (Tier 1-3 only, gene-annotated)
  vt <- v[!is.na(gene) & gene != "" & tier %in% c(1L, 2L, 3L),
          .(chrom = as.character(chrom),
            pos   = as.integer(pos),
            ref   = as.character(ref),
            alt   = as.character(alt),
            gene  = as.character(gene))]
  if (nrow(vt) == 0L) {
    message("[fig] skipping sample_gene_hotspots (no gene-annotated Tier 1-3 variants)")
    return(invisible())
  }
  keys_tmp <- tempfile(fileext = ".tsv")
  fwrite(vt, keys_tmp, sep = "\t", col.names = FALSE)
  on.exit(unlink(keys_tmp), add = TRUE)

  # Stream genotypes: $1=chrom, $2=pos, $3=ref, $4=alt, $5=sample, $6=gt
  # For each carrier (het or hom_alt) line where variant is in our Tier 1-3 set,
  # increment count for (sample, gene). Emit (sample, gene, count).
  awk_cmd <- sprintf(
    "awk -F'\\t' '
       NR == FNR { vg[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = $5; next }
       FNR == 1  { next }
       $6 ~ /^(0\\/1|1\\/0|0\\|1|1\\|0|1\\/1|1\\|1)$/ {
         k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
         if (k in vg) cnt[$5 SUBSEP vg[k]]++
       }
       END {
         for (k in cnt) {
           split(k, a, SUBSEP)
           if (cnt[k] >= 2) printf \"%%s\\t%%s\\t%%d\\n\", a[1], a[2], cnt[k]
         }
       }
     ' %s %s",
    shQuote(keys_tmp), shQuote(gt_path)
  )

  message("[fig] streaming genotypes for sample-gene hotspots...")
  pairs <- tryCatch(
    fread(cmd = awk_cmd, header = FALSE,
          col.names = c("sample", "gene", "n_variants")),
    error = function(e) { message("[fig] hotspot awk failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(pairs) || nrow(pairs) == 0L) {
    message("[fig] skipping sample_gene_hotspots (no (sample, gene) pair has >=2 rare Tier 1-3 variants)")
    return(invisible())
  }
  setorder(pairs, -n_variants, gene, sample)
  full_n   <- nrow(pairs)                       # capture BEFORE cropping for caption
  full_max <- max(pairs$n_variants)
  pairs <- head(pairs, top_n)
  pairs[, label := sprintf("%-12s  carrier:  %s", gene, sample)]
  pairs[, label := factor(label, levels = rev(unique(label)))]

  # Color by anchor-gene membership to make NF anchors pop visually
  NF_ANCHORS <- c("NF1", "NF2", "SMARCB1", "LZTR1")
  pairs[, anchor := ifelse(gene %in% NF_ANCHORS, "NF anchor gene", "Other")]
  pairs[, anchor := factor(anchor, levels = c("NF anchor gene", "Other"))]
  anchor_colors <- c("NF anchor gene" = "#c0392b", "Other" = "#8e44ad")

  p <- ggplot(pairs, aes(x = n_variants, y = label, fill = anchor)) +
    geom_col() +
    geom_text(aes(label = paste0(n_variants, " variants")),
              hjust = -0.12, size = 3.6, color = "grey20") +
    scale_fill_manual(values = anchor_colors, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22)),
                       breaks = pretty_breaks(n = 6)) +
    labs(
      title    = sprintf("Top %d of %s sample-gene hotspots (>=2 rare Tier 1-3 variants)",
                         nrow(pairs), comma(full_n)),
      subtitle = sprintf("%s -- compound-heterozygous / multi-hit candidates  (max in cohort: %d variants)",
                         cohort_name, full_max),
      x = "Distinct rare Tier 1-3 variants in same gene, same sample",
      y = NULL,
      caption = paste(
        "Each row: one sample carrying multiple distinct rare predicted-damaging",
        "variants in the same gene. Candidates for compound-heterozygous recessive",
        "effects (if both het) or multi-hit somatic two-hit events. Phase by trio",
        "genotypes or long-read sequencing before clinical interpretation.",
        "NF anchor-gene hotspots highlighted in red.",
        sep = "\n"
      )
    ) +
    THEME_PRES +
    theme(axis.text.y = element_text(size = 10, family = "mono"),
          legend.position = "top")

  ggsave(out_path, p, width = 13, height = 0.34 * nrow(pairs) + 2.8,
         dpi = 300, bg = "white")
  message("[fig] wrote ", out_path)
}

# ------------------------------- main ----------------------------------------
main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  message(sprintf("[cohort_figures] variants=%s genes=%s out_dir=%s",
                  opt$variants, opt$genes, opt$out_dir))

  v <- fread(opt$variants, sep = "\t", na.strings = c("", ".", "NA"))
  g <- fread(opt$genes,    sep = "\t", na.strings = c("", ".", "NA"))
  # Don't fread the genotype file -- it's 30+ GB at cohort scale and OOMs R.
  # fig_per_sample_burden() streams it through awk via system() instead.
  gt_path <- if (!is.null(opt$genotypes) && file.exists(opt$genotypes)) opt$genotypes else NULL

  fig_tier_distribution     (v, opt$cohort_name, file.path(opt$out_dir, "01_tier_distribution.png"))
  fig_variant_class_by_tier (v, opt$cohort_name, file.path(opt$out_dir, "02_variant_class_by_tier.png"))
  fig_per_sample_burden     (v, gt_path, opt$cohort_name, file.path(opt$out_dir, "03_per_sample_burden.png"))
  fig_modifier_landscape    (g, opt$cohort_name, file.path(opt$out_dir, "04_modifier_landscape.png"))
  fig_nf_anchor_landscape   (v, opt$cohort_name, file.path(opt$out_dir, "05_nf_anchor_landscape.png"))
  fig_variant_recurrence    (v, opt$cohort_name, file.path(opt$out_dir, "06_variant_recurrence.png"))
  fig_top_recurrent_variants(v, opt$cohort_name, file.path(opt$out_dir, "07_top_recurrent_t12.png"))
  fig_sample_gene_hotspots  (v, gt_path, opt$cohort_name, file.path(opt$out_dir, "08_sample_gene_hotspots.png"))

  message("[cohort_figures] done. wrote figures to: ", opt$out_dir)
}

if (sys.nframe() == 0L) main()
