#!/usr/bin/env Rscript
# tier_variants.R
# -----------------------------------------------------------------------------
# Two parallel classifications on a per-(variant, transcript) annotated table:
#   (1) Causation tier 1..5 -- class-conditional pathogenicity gradient,
#       no cross-tool consensus stacking. See docs/tiering.md for the design.
#   (2) Modifier candidate flag -- independent, genome-wide by default, optional
#       --modifier-genes restriction. Same variant can be tier=4 AND modifier.
#
# Inputs:
#   --input       per-(variant, transcript) TSV (required)
#   --out-prefix  output file prefix (required)
#   --mode        "research" (default) | "acmg"
#   --af-max      rarity gate, default 1e-4
#   --genotypes   optional long TSV: chrom, pos, ref, alt, sample, gt
#   --modifier-genes  optional one-per-line file (HGNC symbols or Ensembl IDs)
#
# Outputs:
#   <out-prefix>.variants.tsv   one row per variant, tier + modifier + cohort
#   <out-prefix>.genes.tsv      per-gene summary with description + burden
#
# Authored without an R toolchain in the sandbox -- tested by mirroring the
# logic in Python against tests/example_annotated_tx.tsv (see docs/tiering.md).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

# ============================ configuration ==================================
CFG <- list(
  mode               = "research",
  af_max             = 1e-4,
  af_modifier_upper  = 0.05,
  loeuf_constrained  = 0.35,
  pli_hi             = 0.90,
  am_path            = 0.564,
  am_strong          = 0.85,
  am_benign          = 0.34,
  revel_pp3_supp     = 0.644,
  revel_pp3_mod      = 0.773,
  revel_pp3_strong   = 0.932,
  revel_bp4          = 0.290,
  esm1b_damaging     = -7.5,
  spliceai_flag      = 0.20,
  spliceai_likely    = 0.50,
  spliceai_high      = 0.80,
  acmg_primary       = "alphamissense"
)

# ----------------------------- ontologies ------------------------------------
SO_RANK <- c(
  "transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained",
  "frameshift_variant","stop_lost","start_lost","transcript_amplification",
  "inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant",
  "splice_region_variant","splice_donor_5th_base_variant","splice_donor_region_variant",
  "splice_polypyrimidine_tract_variant","incomplete_terminal_codon_variant","start_retained_variant",
  "stop_retained_variant","synonymous_variant","coding_sequence_variant","mature_miRNA_variant",
  "5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant",
  "NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant",
  "downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant",
  "regulatory_region_ablation","regulatory_region_amplification","feature_elongation",
  "regulatory_region_variant","feature_truncation","intergenic_variant"
)

LOF_TERMS   <- c("transcript_ablation","splice_acceptor_variant","splice_donor_variant",
                 "stop_gained","frameshift_variant","start_lost")
SPLICE_NCAN <- c("splice_region_variant","splice_donor_5th_base_variant",
                 "splice_donor_region_variant","splice_polypyrimidine_tract_variant")
MISSENSE    <- c("missense_variant","protein_altering_variant")
INFRAME     <- c("inframe_insertion","inframe_deletion")
NONCODING   <- c("synonymous_variant","start_retained_variant","stop_retained_variant",
                 "incomplete_terminal_codon_variant","5_prime_UTR_variant","3_prime_UTR_variant",
                 "intron_variant","non_coding_transcript_exon_variant",
                 "non_coding_transcript_variant","mature_miRNA_variant","coding_sequence_variant")
REGULATORY  <- c("regulatory_region_variant","TFBS_ablation","TFBS_amplification",
                 "TF_binding_site_variant","regulatory_region_ablation",
                 "regulatory_region_amplification")

# ------------------------ input column mapping --------------------------------
# Edit the RHS if your annotated table uses different headers.
COLS <- list(
  chrom="chrom", pos="pos", ref="ref", alt="alt",
  gene="gene", gene_id="gene_id", transcript="transcript",
  mane_select="mane_select", canonical="canonical",
  biotype="biotype",
  consequence="consequence", impact="impact",
  exon="exon", intron="intron",
  hgvsg="hgvsg", hgvsc="hgvsc", hgvsp="hgvsp",
  existing_variation="existing_variation",
  domains="domains", variant_class="variant_class",
  pheno="pheno", gene_pheno="gene_pheno",
  motif_name="motif_name",
  loftee_lof="loftee_lof",
  gnomad_af="gnomad_af", gnomad_faf="gnomad_faf", gnomad_popmax_af="gnomad_popmax_af",
  clin_sig="clin_sig", clin_stars="clin_stars",
  loeuf="loeuf", pli="pli", mis_z="mis_z", syn_z="syn_z",
  omim_phenotype="omim_phenotype",
  spliceai_ds_max="spliceai_ds_max",
  alphamissense="alphamissense", am_class="am_class",
  esm1b="esm1b", revel="revel",
  sift="sift", polyphen="polyphen",
  nmd_escape="nmd_escape"
)

# ============================= helpers ========================================
num_col  <- function(dt, col) suppressWarnings(as.numeric(dt[[col]]))
safe_min <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x) else NA_real_ }
safe_max <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }
T0       <- function(x) x %in% TRUE   # NA-safe TRUE
truthy   <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  tolower(as.character(x)) %in% c("true","t","yes","y","1")
}
`%||%` <- function(a, b) if (is.null(a) || all(is.na(a)) || all(a == "")) b else a

parse_args <- function(args) {
  out <- list(input=NULL, out_prefix="tiered", mode=CFG$mode, af_max=CFG$af_max,
              genotypes=NULL, modifier_genes=NULL)
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if      (a == "--input")           { out$input <- args[[i+1]]; i <- i+2 }
    else if (a == "--out-prefix")      { out$out_prefix <- args[[i+1]]; i <- i+2 }
    else if (a == "--mode")            { out$mode <- args[[i+1]]; i <- i+2 }
    else if (a == "--af-max")          { out$af_max <- as.numeric(args[[i+1]]); i <- i+2 }
    else if (a == "--genotypes")       { out$genotypes <- args[[i+1]]; i <- i+2 }
    else if (a == "--modifier-genes")  { out$modifier_genes <- args[[i+1]]; i <- i+2 }
    else stop(sprintf("unknown argument: %s", a))
  }
  if (is.null(out$input)) stop("--input is required")
  out
}

so_rank_vec <- function(conseq) {
  vapply(strsplit(as.character(conseq), "[&,]"), function(terms) {
    r <- match(trimws(terms), SO_RANK)
    if (all(is.na(r))) length(SO_RANK) + 1L else as.integer(min(r, na.rm = TRUE))
  }, integer(1))
}
most_severe_term <- function(conseq) {
  vapply(strsplit(as.character(conseq), "[&,]"), function(terms) {
    terms <- trimws(terms); r <- match(terms, SO_RANK)
    if (all(is.na(r))) terms[1] else terms[which.min(r)]
  }, character(1))
}
classify_term <- function(term) {
  fcase(
    term %in% LOF_TERMS,   "pLoF",
    term %in% SPLICE_NCAN, "splice_noncanonical",
    term %in% MISSENSE,    "missense",
    term %in% INFRAME,     "inframe",
    term %in% NONCODING,   "noncoding",
    term %in% REGULATORY,  "regulatory",
    default = "other"
  )
}

# Best-effort NMD-escape heuristic from VEP EXON "i/N" cell. Returns TRUE only if
# last-exon (i == N). Better-grade signal: pass a precomputed nmd_escape column.
parse_last_exon <- function(exon_cell) {
  vapply(as.character(exon_cell), function(s) {
    if (is.na(s) || s == "" || s == ".") return(NA)
    p <- strsplit(s, "/", fixed = TRUE)[[1]]
    if (length(p) != 2) return(NA)
    idx <- suppressWarnings(as.integer(p[1])); tot <- suppressWarnings(as.integer(p[2]))
    if (is.na(idx) || is.na(tot)) return(NA)
    isTRUE(idx == tot)
  }, logical(1), USE.NAMES = FALSE)
}

ensure_cols <- function(dt) {
  for (nm in unlist(COLS)) if (!nm %in% names(dt)) dt[, (nm) := NA]
  dt
}

# =========================== load + collapse ==================================
collapse_to_variant <- function(dt) {
  C <- COLS
  dt[, sev_rank := so_rank_vec(get(C$consequence))]
  dt[, is_mane  := truthy(get(C$mane_select))]
  dt[, is_canon := truthy(get(C$canonical))]
  key <- c(C$chrom, C$pos, C$ref, C$alt)

  # Worst-across-any-transcript score (kept; not used for primary tier).
  dt[, am_all_helper  := num_col(dt, C$alphamissense)]
  dt[, esm_all_helper := num_col(dt, C$esm1b)]
  worst <- dt[, .(max_alphamissense_any_tx = safe_max(am_all_helper),
                  min_esm1b_any_tx         = safe_min(esm_all_helper)), by = key]

  setorderv(dt, c(key, "is_mane", "is_canon", "sev_rank"),
                c(rep(1L, length(key)), -1L, -1L, 1L))
  rep_tx <- dt[, .SD[1L], by = key]
  rep_tx[, c("sev_rank","is_mane","is_canon","am_all_helper","esm_all_helper") := NULL]
  rep_tx <- merge(rep_tx, worst, by = key, all.x = TRUE)

  rep_tx[, ms_term := most_severe_term(get(C$consequence))]
  rep_tx[, variant_class := classify_term(ms_term)]
  rep_tx[]
}

# ============================ evidence flags ==================================
add_evidence <- function(v) {
  C <- COLS

  # Rarity: prefer FAF, fall back to popmax, fall back to overall gnomAD.
  faf  <- num_col(v, C$gnomad_faf)
  pmax <- num_col(v, C$gnomad_popmax_af)
  af   <- num_col(v, C$gnomad_af)
  af_used <- ifelse(!is.na(faf), faf, ifelse(!is.na(pmax), pmax, af))
  v[, af_used := af_used]
  v[, rare := is.na(af_used) | af_used <= CFG$af_max]

  v[, is_lof := variant_class == "pLoF"]
  v[, constrained := T0(num_col(v, C$loeuf) <  CFG$loeuf_constrained) |
                     T0(num_col(v, C$pli)   >= CFG$pli_hi)]
  # NMD escape: precomputed boolean wins; else last-exon heuristic from VEP EXON cell.
  pre_nmd_raw <- v[[C$nmd_escape]]
  pre_nmd <- if (is.logical(pre_nmd_raw)) pre_nmd_raw else
             tolower(as.character(pre_nmd_raw)) %in% c("true","t","yes","y","1")
  pre_nmd[is.na(pre_nmd_raw)] <- NA
  heur <- parse_last_exon(v[[C$exon]])
  v[, nmd_escape := T0(ifelse(!is.na(pre_nmd), pre_nmd, heur))]

  loftee <- toupper(as.character(v[[C$loftee_lof]]))
  v[, loftee_hc := loftee == "HC" & !is.na(loftee)]
  v[, loftee_lc := loftee == "LC" & !is.na(loftee)]
  v[, loftee_present := loftee %in% c("HC","LC")]

  am  <- num_col(v, C$alphamissense)
  amc <- tolower(as.character(v[[C$am_class]]))
  v[, am_path   := T0(am > CFG$am_path)   | grepl("pathogenic", amc)]
  v[, am_strong := T0(am >= CFG$am_strong)]
  v[, am_ambig  := T0(am > CFG$am_benign & am <= CFG$am_path)]
  v[, am_benign := T0(am < CFG$am_benign) | grepl("benign", amc)]

  revel <- num_col(v, C$revel); esm <- num_col(v, C$esm1b)
  v[, revel_path_supp := T0(revel >= CFG$revel_pp3_supp)]
  v[, revel_path_mod  := T0(revel >= CFG$revel_pp3_mod)]
  v[, revel_benign    := T0(revel <= CFG$revel_bp4)]
  v[, esm_damaging    := T0(esm <= CFG$esm1b_damaging)]

  sa <- num_col(v, C$spliceai_ds_max)
  v[, splice_flag   := T0(sa >= CFG$spliceai_flag)]
  v[, splice_likely := T0(sa >= CFG$spliceai_likely)]
  v[, splice_high   := T0(sa >= CFG$spliceai_high)]

  cs    <- tolower(as.character(v[[C$clin_sig]]))
  stars <- num_col(v, C$clin_stars)
  v[, clinvar_plp := grepl("pathogenic", cs) & !grepl("benign", cs) &
                     (T0(stars >= 1) | is.na(stars))]
  v[, clinvar_blb := grepl("benign", cs) & !grepl("pathogenic", cs)]

  motif <- as.character(v[[C$motif_name]])
  v[, has_motif := !is.na(motif) & motif != "" & motif != "."]
  v[]
}

# ===================== causation tier (class-conditional) =====================
assign_tier_research <- function(v) {
  v[, tier := 4L]
  v[, tier_reason := "rare; no class-appropriate signal"]

  # Snapshot evidence columns as locals. data.table's `i` evaluation looks up
  # single-symbol or `!sym` filters in calling scope before column scope; using
  # locals avoids the "is.found-in-calling-scope-but-also-a-column" error and
  # also makes the row filters explicitly NA-safe (every column here is logical
  # with NA already coerced to FALSE upstream in add_evidence).
  rare           <- v$rare
  constrained    <- v$constrained
  nmd_escape     <- v$nmd_escape
  loftee_lc      <- v$loftee_lc
  loftee_hc      <- v$loftee_hc
  loftee_present <- v$loftee_present
  splice_high    <- v$splice_high
  splice_likely  <- v$splice_likely
  splice_flag    <- v$splice_flag
  am_strong      <- v$am_strong
  am_path        <- v$am_path
  am_ambig       <- v$am_ambig
  am_benign      <- v$am_benign
  has_motif      <- v$has_motif
  clinvar_plp    <- v$clinvar_plp
  clinvar_blb    <- v$clinvar_blb
  variant_class  <- v$variant_class

  # ---- pLoF ----
  is_lof <- variant_class == "pLoF"
  hc_or_unknown <- !loftee_present | loftee_hc
  v[is_lof & rare & constrained & !nmd_escape & hc_or_unknown,
    `:=`(tier = 1L, tier_reason = "pLoF in constrained gene; not NMD-escape (LOFTEE HC if present)")]
  v[is_lof & rare & (!constrained | nmd_escape | loftee_lc),
    `:=`(tier = 2L, tier_reason = "pLoF in non-constrained gene / NMD-escape / LOFTEE LC")]

  # ---- non-canonical splice ----
  is_sp <- variant_class == "splice_noncanonical"
  v[is_sp & rare & splice_high,                       `:=`(tier = 1L, tier_reason = "non-canonical splice + SpliceAI>=0.8")]
  v[is_sp & rare & splice_likely & !splice_high,      `:=`(tier = 2L, tier_reason = "non-canonical splice + SpliceAI 0.5-0.8")]
  v[is_sp & rare & splice_flag & !splice_likely,      `:=`(tier = 3L, tier_reason = "non-canonical splice + SpliceAI 0.2-0.5")]

  # ---- missense ----
  is_ms <- variant_class == "missense"
  am_num <- num_col(v, COLS$alphamissense)
  am_missing <- is.na(am_num)
  v[is_ms & rare & am_strong & constrained,
    `:=`(tier = 1L, tier_reason = "missense: AM>=0.85 in constrained gene")]
  v[is_ms & rare & am_path & !(am_strong & constrained),
    `:=`(tier = 2L, tier_reason = "missense: AM likely_pathogenic")]
  v[is_ms & rare & am_ambig & !am_path,
    `:=`(tier = 3L, tier_reason = "missense: AM ambiguous (0.34-0.564)")]
  v[is_ms & rare & am_missing & constrained,
    `:=`(tier = 3L, tier_reason = "rare missense in constrained gene; AM missing")]

  # ---- in-frame ----
  is_if <- variant_class == "inframe"
  loeuf_v <- num_col(v, COLS$loeuf)
  v[is_if & rare & (T0(loeuf_v < CFG$loeuf_constrained) | splice_likely),
    `:=`(tier = 2L, tier_reason = "in-frame in constrained gene or SpliceAI>=0.5")]
  # `tier` here is read fresh from v so it reflects any earlier := updates.
  v[is_if & rare & (constrained | splice_flag) & v$tier > 2L,
    `:=`(tier = 3L, tier_reason = "in-frame in moderately constrained or SpliceAI 0.2-0.5")]

  # ---- non-coding ----
  is_nc <- variant_class == "noncoding"
  v[is_nc & rare & splice_likely,                      `:=`(tier = 2L, tier_reason = "non-coding + SpliceAI>=0.5 (cryptic splice)")]
  v[is_nc & rare & splice_flag & !splice_likely,       `:=`(tier = 3L, tier_reason = "non-coding + SpliceAI 0.2-0.5")]

  # ---- regulatory ----
  is_reg <- variant_class == "regulatory"
  v[is_reg & rare & has_motif & splice_likely,         `:=`(tier = 2L, tier_reason = "regulatory: motif + SpliceAI>=0.5")]
  v[is_reg & rare & (has_motif | splice_flag),         `:=`(tier = 3L, tier_reason = "regulatory: motif or SpliceAI 0.2-0.5")]

  # ---- benign / common overrides (applied late) ----
  v[is_ms & am_benign,        `:=`(tier = 5L, tier_reason = "missense: AM likely_benign")]
  v[!rare,                    `:=`(tier = 5L, tier_reason = "common (failed rarity gate)")]
  v[clinvar_blb,              `:=`(tier = 5L, tier_reason = "ClinVar B/LB")]
  v[clinvar_plp,              `:=`(tier = 1L, tier_reason = "ClinVar P/LP (>=1 star)")]
  v[]
}

# =========================== ACMG-strict mode ================================
assign_class_acmg <- function(v) {
  # Snapshot evidence columns as locals (same reason as in assign_tier_research).
  am_path       <- v$am_path
  am_benign     <- v$am_benign
  revel_benign  <- v$revel_benign
  splice_likely <- v$splice_likely
  is_lof        <- v$is_lof
  constrained   <- v$constrained
  rare          <- v$rare
  nmd_escape    <- v$nmd_escape
  clinvar_plp   <- v$clinvar_plp
  clinvar_blb   <- v$clinvar_blb
  revel         <- num_col(v, COLS$revel)

  v[, pp3 := "none"]
  if (CFG$acmg_primary == "alphamissense") {
    v[am_path, pp3 := "PP3_supporting"]
  } else {
    v[T0(revel >= CFG$revel_pp3_supp),   pp3 := "PP3_supporting"]
    v[T0(revel >= CFG$revel_pp3_mod),    pp3 := "PP3_moderate"]
    v[T0(revel >= CFG$revel_pp3_strong), pp3 := "PP3_strong"]
  }
  v[splice_likely, pp3 := "PP3_splicing"]
  v[, bp4 := fifelse(am_benign | revel_benign, "BP4", "none")]
  v[, acmg_class := "VUS"]
  v[is_lof & constrained & rare & !nmd_escape, acmg_class := "LP_lof"]
  v[clinvar_plp,                              acmg_class := "P_clinvar"]
  v[clinvar_blb | !rare,                      acmg_class := "B/LB"]
  # carry a numeric "tier" alias so the gene summary path can be shared.
  v[, tier := fcase(acmg_class %in% c("P_clinvar","LP_lof"), 1L,
                    acmg_class == "VUS", 4L,
                    default = 5L)]
  v[, tier_reason := acmg_class]
  v[]
}

# ======================= modifier candidate flag =============================
add_modifier <- function(v, modifier_genes = NULL) {
  C <- COLS
  am <- num_col(v, C$alphamissense)
  sa <- num_col(v, C$spliceai_ds_max)
  af <- v$af_used

  in_modifier_band <- !is.na(af) & af > CFG$af_max & af <= CFG$af_modifier_upper
  rare_low_tier    <- v$rare & v$tier %in% c(3L, 4L)
  af_signal        <- in_modifier_band | rare_low_tier

  ms_moderate <- v$variant_class == "missense" &
                 T0(am > CFG$am_benign & am <= CFG$am_strong)
  sp_moderate <- v$variant_class %in% c("splice_noncanonical","noncoding") &
                 T0(sa >= CFG$spliceai_flag & sa < CFG$spliceai_likely)
  any_lof     <- v$variant_class == "pLoF"
  is_reg      <- v$variant_class == "regulatory"
  functional  <- ms_moderate | sp_moderate | any_lof | is_reg

  not_common  <- is.na(af) | af <= CFG$af_modifier_upper
  cand        <- af_signal & functional & not_common

  if (!is.null(modifier_genes)) {
    mg <- toupper(modifier_genes)
    in_set <- toupper(as.character(v[[C$gene]])) %in% mg |
              toupper(as.character(v[[C$gene_id]])) %in% mg
    cand <- cand & in_set
  }

  v[, modifier_candidate := cand]
  ev <- character(nrow(v))
  ev[v$modifier_candidate & ms_moderate & in_modifier_band] <-
    sprintf("missense AM=%.2f, AF=%.3g (modifier band)",
            am[v$modifier_candidate & ms_moderate & in_modifier_band],
            af[v$modifier_candidate & ms_moderate & in_modifier_band])
  ev[v$modifier_candidate & ms_moderate & rare_low_tier & ev == ""] <-
    sprintf("rare missense AM=%.2f (sub-Mendelian)",
            am[v$modifier_candidate & ms_moderate & rare_low_tier & ev == ""])
  ev[v$modifier_candidate & sp_moderate & ev == ""] <-
    sprintf("SpliceAI=%.2f sub-threshold nudge",
            sa[v$modifier_candidate & sp_moderate & ev == ""])
  ev[v$modifier_candidate & any_lof & in_modifier_band & ev == ""] <-
    sprintf("pLoF in modifier AF band (%.3g)",
            af[v$modifier_candidate & any_lof & in_modifier_band & ev == ""])
  ev[v$modifier_candidate & any_lof & ev == ""] <- "pLoF (rare, not Mendelian top tier)"
  ev[v$modifier_candidate & is_reg & ev == ""]  <- "regulatory variant + AF signal"
  v[, modifier_evidence := ev]
  v[]
}

# ============================ cohort summaries ===============================
read_genotypes <- function(path) {
  gt <- fread(path, sep = "\t", na.strings = c("","NA","./.",".|."))
  setnames(gt, tolower(names(gt)))
  needed <- c("chrom","pos","ref","alt","sample","gt")
  miss <- setdiff(needed, names(gt))
  if (length(miss)) stop(sprintf("genotype file missing columns: %s", paste(miss, collapse=", ")))
  gt[, alt_count := fcase(
    gt %in% c("0/1","1/0","0|1","1|0"), 1L,
    gt %in% c("1/1","1|1"),             2L,
    gt %in% c("0/0","0|0"),             0L,
    default = NA_integer_)]
  gt[, called := !is.na(alt_count)]
  gt[, is_het := gt %in% c("0/1","1/0","0|1","1|0")]
  gt[, is_hom := gt %in% c("1/1","1|1")]
  gt[, pos := as.integer(pos)]
  gt[, chrom := as.character(chrom)]
  gt[]
}

per_variant_cohort <- function(gt) {
  gt[, .(
    cohort_ac  = sum(alt_count, na.rm = TRUE),
    cohort_an  = 2L * sum(called, na.rm = TRUE),
    n_het      = sum(is_het, na.rm = TRUE),
    n_hom      = sum(is_hom, na.rm = TRUE),
    n_carriers = sum(T0(alt_count > 0L), na.rm = TRUE)
  ), by = .(chrom, pos, ref, alt)][, cohort_af := ifelse(cohort_an > 0L, cohort_ac / cohort_an, NA_real_)][]
}

gene_sample_burden <- function(v, gt) {
  C <- COLS
  carriers <- gt[T0(alt_count > 0L), .(chrom, pos = as.integer(pos), ref, alt, sample)]
  v_keys <- v[, .(chrom = as.character(get(C$chrom)),
                  pos   = as.integer(get(C$pos)),
                  ref   = as.character(get(C$ref)),
                  alt   = as.character(get(C$alt)),
                  gene  = as.character(get(C$gene)),
                  tier, modifier_candidate,
                  impact = as.character(get(C$impact)),
                  variant_class)]
  carriers[, chrom := as.character(chrom)]
  m <- merge(carriers, v_keys, by = c("chrom","pos","ref","alt"))
  # Default qualifying mask: Tier 1-2 in any class, OR any pLoF at Tier <= 3
  # (catches NMD-escape / LOFTEE-LC LoF that get demoted to T2/T3).
  m[, qualifying := T0(tier %in% c(1L,2L)) |
                    (variant_class == "pLoF" & T0(tier <= 3L))]
  m[, .(
    n_samples_high_impact            = uniqueN(sample[impact == "HIGH"]),
    n_samples_tier1                  = uniqueN(sample[T0(tier == 1L)]),
    n_samples_tier12                 = uniqueN(sample[T0(tier %in% c(1L,2L))]),
    n_samples_qualifying             = uniqueN(sample[qualifying]),
    n_samples_with_modifier_candidate = uniqueN(sample[T0(modifier_candidate)])
  ), by = gene]
}

# ============================ gene description ===============================
gene_description_block <- function(v) {
  C <- COLS
  v[, gene_key := as.character(get(C$gene))]
  v[, biotype_str := as.character(get(C$biotype))]
  v[, gene_id_str := as.character(get(C$gene_id))]
  v[, omim_str    := as.character(get(C$omim_phenotype))]
  v[, loeuf_n     := num_col(v, C$loeuf)]
  v[, pli_n       := num_col(v, C$pli)]
  v[, misz_n      := num_col(v, C$mis_z)]
  v[, synz_n      := num_col(v, C$syn_z)]

  desc <- v[, .(
    SYMBOL         = gene_key[1],
    Gene           = first_non_na(gene_id_str),
    BIOTYPE        = first_non_na(biotype_str),
    OMIM_phenotype = first_non_na(omim_str),
    LOEUF          = safe_min(loeuf_n),
    pLI            = safe_max(pli_n),
    mis_z          = safe_max(misz_n),
    syn_z          = safe_max(synz_n)
  ), by = gene_key]

  desc[, gene_description := sprintf(
    "%s%s%s%s%s",
    SYMBOL,
    ifelse(is.na(BIOTYPE) | BIOTYPE == "", "", sprintf(" [%s]", BIOTYPE)),
    ifelse(is.na(Gene)    | Gene == "",    "", sprintf(" %s", Gene)),
    ifelse(is.na(LOEUF), "", sprintf(" | LOEUF=%.2f", LOEUF)),
    ifelse(is.na(OMIM_phenotype) | OMIM_phenotype == "", "",
           sprintf(" | OMIM: %s", OMIM_phenotype))
  )]
  desc
}

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != "" & x != "."]
  if (length(x)) x[1] else NA_character_
}

# ============================ per-gene summary ===============================
gene_summary <- function(v_in, mode, gt = NULL) {
  C <- COLS
  v <- copy(v_in)
  v[, gene_key := as.character(get(C$gene))]
  v[, am_n     := num_col(v, C$alphamissense)]
  v[, esm_n    := num_col(v, C$esm1b)]
  v[, revel_n  := num_col(v, C$revel)]
  v[, sa_n     := num_col(v, C$spliceai_ds_max)]

  g <- v[, .(
    n_variants             = .N,
    n_lof                  = sum(variant_class == "pLoF"),
    n_splice_ge05          = sum(T0(sa_n >= CFG$spliceai_likely)),
    n_clinvar_plp          = sum(T0(clinvar_plp)),
    n_modifier_candidates  = sum(T0(modifier_candidate)),
    max_alphamissense      = safe_max(am_n),
    min_esm1b              = safe_min(esm_n),
    max_revel              = safe_max(revel_n),
    max_spliceai           = safe_max(sa_n),
    min_gnomad_af          = safe_min(af_used)
  ), by = gene_key]

  if (mode != "acmg") {
    tc <- dcast(v[, .N, by = .(gene_key, tier)], gene_key ~ tier, value.var = "N", fill = 0L)
    for (k in as.character(1:5)) if (!k %in% names(tc)) tc[, (k) := 0L]
    setnames(tc, as.character(1:5), paste0("n_tier", 1:5), skip_absent = TRUE)
    keep <- c("gene_key", paste0("n_tier", 1:5))
    tc <- tc[, ..keep]
    g <- merge(g, tc, by = "gene_key", all.x = TRUE)
    ntcols <- paste0("n_tier", 1:5)
    g[, (ntcols) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))), .SDcols = ntcols]
    g[, best_tier := apply(.SD, 1L, function(r) {
        w <- which(r > 0); if (length(w)) min(w) else NA_integer_ }), .SDcols = ntcols]
  }

  # gene description
  desc <- gene_description_block(copy(v))
  g <- merge(g, desc, by = "gene_key", all.x = TRUE)

  # cohort sample burden
  if (!is.null(gt)) {
    burden <- gene_sample_burden(v, gt)
    g <- merge(g, burden, by.x = "gene_key", by.y = "gene", all.x = TRUE)
    burden_cols <- c("n_samples_high_impact","n_samples_tier1","n_samples_tier12",
                     "n_samples_qualifying","n_samples_with_modifier_candidate")
    g[, (burden_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))),
      .SDcols = burden_cols]
  }

  setnames(g, "gene_key", C$gene)
  if (mode != "acmg") {
    setorderv(g, c("best_tier","n_tier1","n_tier2","n_modifier_candidates"),
                 c(1L,-1L,-1L,-1L))
  } else {
    setorderv(g, "n_clinvar_plp", -1L)
  }
  g[]
}

# =================================== main =====================================
main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  CFG$mode   <<- opt$mode
  CFG$af_max <<- opt$af_max
  message(sprintf("[tier_variants] mode=%s  af_max=%g  input=%s  genotypes=%s  modifier-genes=%s",
                  CFG$mode, CFG$af_max, opt$input,
                  ifelse(is.null(opt$genotypes), "<none>", opt$genotypes),
                  ifelse(is.null(opt$modifier_genes), "<none>", opt$modifier_genes)))

  dt <- fread(opt$input, sep = "\t", na.strings = c("",".","NA"))
  dt <- ensure_cols(dt)
  # Force chrom to character so merges with the genotype table don't fail on
  # chr-21 (integer) vs chrX (character) type mismatches.
  dt[, (COLS$chrom) := as.character(get(COLS$chrom))]

  vv <- collapse_to_variant(dt)
  vv <- add_evidence(vv)
  vv <- if (CFG$mode == "acmg") assign_class_acmg(vv) else assign_tier_research(vv)

  mg <- NULL
  if (!is.null(opt$modifier_genes)) {
    lines <- readLines(opt$modifier_genes)
    lines <- trimws(lines)
    mg <- lines[nzchar(lines) & !startsWith(lines, "#")]
    message(sprintf("[tier_variants] modifier-genes: %d symbols loaded", length(mg)))
  }
  vv <- add_modifier(vv, modifier_genes = mg)

  gt <- NULL
  if (!is.null(opt$genotypes)) {
    gt <- read_genotypes(opt$genotypes)
    coh <- per_variant_cohort(gt)
    vv <- merge(vv, coh, by.x = c(COLS$chrom, COLS$pos, COLS$ref, COLS$alt),
                          by.y = c("chrom","pos","ref","alt"), all.x = TRUE)
  }

  desc <- gene_description_block(copy(vv))
  vv <- merge(vv, desc[, .(gene_key, gene_description)],
              by.x = COLS$gene, by.y = "gene_key", all.x = TRUE)

  gs <- gene_summary(vv, CFG$mode, gt = gt)

  vfile <- paste0(opt$out_prefix, ".variants.tsv")
  gfile <- paste0(opt$out_prefix, ".genes.tsv")
  fwrite(vv, vfile, sep = "\t")
  fwrite(gs, gfile, sep = "\t")
  message(sprintf("[tier_variants] wrote %s (%d variants) and %s (%d genes)",
                  vfile, nrow(vv), gfile, nrow(gs)))
}

if (sys.nframe() == 0L) main()
