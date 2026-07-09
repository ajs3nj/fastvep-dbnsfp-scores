#!/usr/bin/env Rscript
# plot_af_source_stratification.R
# -----------------------------------------------------------------------------
# Stratified histogram of af_used colored by af_source. Documents which AF
# source carried the rarity information for the run and where the rarity /
# modifier thresholds sit relative to the actual distribution.
#
# Why this matters for the NF1 GWAS run specifically: gnomAD per-variant
# sources weren't loaded in the fastVEP SA dir for this cohort, so af_used
# is populated entirely by the cohort_af fallback. A reviewer looking at
# the tier table needs to know that up front. This figure makes it visible.
#
# Also visualizes the adaptive rarity threshold (3/(2N) for cohort_af source
# vs 1e-4 for gnomAD sources), so the "why is Tier 1 bigger than usual"
# question has an immediate answer.
#
# Inputs:
#   --variants   cohort.variants.tsv  (needs af_used, af_source columns)
#   --out-dir    output directory
#   --cohort-name label for figure title
#   --af-max     rarity threshold used (default: 1e-4)
#   --af-mod-upper  modifier band upper (default: 5e-2)
#
# Outputs:
#   af_source_distribution.png      Stacked histogram of log10(af_used) by
#                                    af_source, with threshold lines.
#   af_source_summary.tsv            One row per af_source: n variants,
#                                    min/median/max af_used, contribution %.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ----------------------------- args ------------------------------------------
parse_args <- function(args) {
  out <- list(variants=NULL, out_dir=NULL, cohort_name="cohort",
              af_max=1e-4, af_mod_upper=5e-2)
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if      (a == "--variants")      { out$variants  <- args[[i+1]]; i <- i+2 }
    else if (a == "--out-dir")       { out$out_dir   <- args[[i+1]]; i <- i+2 }
    else if (a == "--cohort-name")   { out$cohort_name <- args[[i+1]]; i <- i+2 }
    else if (a == "--af-max")        { out$af_max <- as.numeric(args[[i+1]]); i <- i+2 }
    else if (a == "--af-mod-upper")  { out$af_mod_upper <- as.numeric(args[[i+1]]); i <- i+2 }
    else stop(sprintf("unknown argument: %s", a))
  }
  if (is.null(out$variants)) stop("--variants is required")
  if (is.null(out$out_dir))  out$out_dir <- file.path(dirname(out$variants), "figures")
  dir.create(out$out_dir, recursive = TRUE, showWarnings = FALSE)
  out
}

# ----------------------------- theme + colors --------------------------------
# Source colors: gnomAD sources in a blue-teal cascade (preferred to less
# preferred), cohort_af in orange (fallback flag), none in grey (missing).
SOURCE_ORDER <- c("gnomad_faf",
                  "gnomad_popmax_af",
                  "gnomad_af",
                  "cohort_af (gnomAD sources missing)",
                  "none")
SOURCE_LABELS <- c(
  "gnomad_faf"                                = "gnomAD FAF (preferred)",
  "gnomad_popmax_af"                          = "gnomAD popmax AF",
  "gnomad_af"                                 = "gnomAD overall AF",
  "cohort_af (gnomAD sources missing)"        = "Cohort AF (fallback)",
  "none"                                      = "None (NA)"
)
SOURCE_COLORS <- c(
  "gnomad_faf"                                = "#1b6ca8",   # deep blue
  "gnomad_popmax_af"                          = "#2c8fbc",   # medium blue
  "gnomad_af"                                 = "#5db4c9",   # light teal
  "cohort_af (gnomAD sources missing)"        = "#e67e22",   # orange (fallback)
  "none"                                      = "#95a5a6"    # grey (missing)
)

THEME_PRES <- theme_minimal(base_size = 15) +
  theme(
    plot.title       = element_text(face = "bold", size = 18),
    plot.subtitle    = element_text(size = 13, color = "grey30"),
    plot.caption     = element_text(size = 10, color = "grey45", hjust = 0),
    axis.title       = element_text(face = "bold", size = 14),
    axis.text        = element_text(size = 12),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

# ----------------------------- main ------------------------------------------
main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  message(sprintf("[af-fig] variants=%s out-dir=%s", opt$variants, opt$out_dir))

  # Load only the columns we need. af_max_effective is optional (added by the
  # adaptive-threshold patch); handle both old and new tier tables.
  hdr <- names(fread(opt$variants, nrows = 0, sep = "\t"))
  cols <- intersect(c("af_used","af_source","af_max_effective","variant_class","tier"), hdr)
  v <- fread(opt$variants, sep = "\t", select = cols)
  message(sprintf("[af-fig] loaded %s variants", comma(nrow(v))))

  if (!"af_source" %in% names(v)) stop("cohort.variants.tsv missing af_source column")
  if (!"af_used"   %in% names(v)) stop("cohort.variants.tsv missing af_used column")

  # ----- summary table --------------------------------------------------------
  # Normalize af_source values: any variant of "cohort_af..." maps to the
  # canonical label so we don't get 2+ rows if the label ever varies.
  v[, af_source_canonical := ifelse(grepl("^cohort_af", af_source),
                                     "cohort_af (gnomAD sources missing)",
                                     af_source)]
  v[, af_source_canonical := ifelse(af_source_canonical %in% SOURCE_ORDER,
                                     af_source_canonical, "none")]
  v[, af_source_canonical := factor(af_source_canonical,
                                     levels = SOURCE_ORDER)]

  af_num <- suppressWarnings(as.numeric(v$af_used))
  summary_dt <- v[, .(
      n_variants   = .N,
      pct_of_total = .N / nrow(v) * 100,
      n_rare_1e4   = sum(!is.na(af_num) & af_num <= 1e-4, na.rm = TRUE) +
                     sum(is.na(af_num), na.rm = TRUE),
      n_na         = sum(is.na(af_num)),
      min_af_used  = suppressWarnings(min(af_num, na.rm = TRUE)),
      median_af    = suppressWarnings(median(af_num, na.rm = TRUE)),
      max_af_used  = suppressWarnings(max(af_num, na.rm = TRUE))
    ),
    by = .(af_source = af_source_canonical)]
  setorder(summary_dt, -n_variants)

  # Clean up infinites in min/max where entire source is NA
  for (col in c("min_af_used","median_af","max_af_used")) {
    summary_dt[!is.finite(get(col)), (col) := NA_real_]
  }

  sum_path <- file.path(opt$out_dir, "af_source_summary.tsv")
  fwrite(summary_dt, sum_path, sep = "\t")
  message("[af-fig] wrote ", sum_path)
  print(summary_dt)

  # ----- histogram data -------------------------------------------------------
  # For plotting: variants with af_used = NA (af_source == "none") get placed
  # in a dedicated leftmost bin labeled "NA". Everything else uses log10.
  # Use a small pseudo-value below the display range to keep them in the
  # geom_histogram bin logic, then relabel the x-axis.
  d <- copy(v)
  d[, af_num := suppressWarnings(as.numeric(af_used))]
  # A left-of-range marker for NA
  na_marker <- -6.5   # placed below x-min of -5
  d[, log_af := ifelse(is.na(af_num), na_marker,
                       ifelse(af_num < 1e-6, log10(1e-6), log10(af_num)))]

  # ----- figure ---------------------------------------------------------------
  # Adaptive threshold for cohort_af: 3/(2N). Detect N from data if af_max_effective
  # is present (uses the max), otherwise skip that annotation.
  adaptive_thresh <- if ("af_max_effective" %in% names(v)) {
    suppressWarnings(max(v[af_source_canonical == "cohort_af (gnomAD sources missing)",
                             as.numeric(af_max_effective)], na.rm = TRUE))
  } else NA_real_

  # Total counts per source for legend annotations
  legend_lbl <- summary_dt[, .(af_source,
                               lbl = sprintf("%s (n=%s, %.1f%%)",
                                             SOURCE_LABELS[as.character(af_source)],
                                             comma(n_variants),
                                             pct_of_total))]
  # Build a named vector for scale_fill_manual labels
  fill_labels <- setNames(legend_lbl$lbl, as.character(legend_lbl$af_source))

  # X-axis: log10 range from -6 to 0, plus NA bin at -6.5
  x_breaks <- c(na_marker, -6, -5, -4, -3, -2, -1, 0)
  x_labels <- c("NA",
                expression(10^-6),
                expression(10^-5),
                expression(10^-4),
                expression(10^-3),
                expression(10^-2),
                expression(10^-1),
                expression(10^0))

  p <- ggplot(d, aes(x = log_af, fill = af_source_canonical)) +
    geom_histogram(binwidth = 0.15, color = "white", linewidth = 0.15) +
    # Threshold lines
    geom_vline(xintercept = log10(opt$af_max), linetype = "dashed",
               color = "grey20", linewidth = 0.6) +
    annotate("text", x = log10(opt$af_max), y = Inf, hjust = -0.05, vjust = 1.5,
             label = sprintf("gnomAD rare gate\n(%.0e)", opt$af_max),
             size = 3.4, color = "grey20", fontface = "italic") +
    geom_vline(xintercept = log10(opt$af_mod_upper), linetype = "dashed",
               color = "grey20", linewidth = 0.6) +
    annotate("text", x = log10(opt$af_mod_upper), y = Inf, hjust = -0.05, vjust = 1.5,
             label = sprintf("Modifier band upper\n(%.0e)", opt$af_mod_upper),
             size = 3.4, color = "grey20", fontface = "italic") +
    scale_fill_manual(values = SOURCE_COLORS, labels = fill_labels,
                      name = "af_source", drop = FALSE,
                      breaks = SOURCE_ORDER) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels,
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05)))

  # Adaptive threshold line (only if present)
  if (is.finite(adaptive_thresh) && adaptive_thresh > opt$af_max) {
    p <- p +
      geom_vline(xintercept = log10(adaptive_thresh),
                 linetype = "dotted", color = "#c0392b", linewidth = 0.7) +
      annotate("text", x = log10(adaptive_thresh), y = Inf,
               hjust = 1.05, vjust = 1.5,
               label = sprintf("Adaptive rare\n(cohort_af: %.1e)",
                               adaptive_thresh),
               size = 3.4, color = "#c0392b", fontface = "italic")
  }

  p <- p + labs(
    title    = "Effective allele frequency (af_used) distribution by source",
    subtitle = sprintf("%s — which AF source carried the rarity information", opt$cohort_name),
    x = expression(paste("af_used  (log"[10], " scale;  NA = af_source \"none\")")),
    y = "Number of variants",
    caption = paste(
      "af_used cascade: gnomad_faf > gnomad_popmax_af > gnomad_af > cohort_af.",
      "When gnomAD per-variant sources aren't loaded, cohort_af is used as fallback.",
      "The adaptive rare gate (3/(2N) for cohort_af) recovers cohort singletons that",
      "the strict 1e-4 gnomAD-scale threshold would demote.",
      sep = "\n"
    )
  ) + THEME_PRES

  fig_path <- file.path(opt$out_dir, "af_source_distribution.png")
  ggsave(fig_path, p, width = 12, height = 6.5, dpi = 300, bg = "white")
  message("[af-fig] wrote ", fig_path)
}

main()
