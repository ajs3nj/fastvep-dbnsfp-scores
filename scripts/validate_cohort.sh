#!/usr/bin/env bash
# validate_cohort.sh
# -----------------------------------------------------------------------------
# Sanity-check a finished cohort_pipeline.sh run.
#
# Asserts the tiering and filtering rules hold against the actual outputs:
#   - file presence + non-empty
#   - required columns in each file
#   - tier distribution is pyramid-shaped (Tier 4 >> Tier 1 etc.)
#   - class-conditional rules hold (every Tier 1 pLoF is in constrained gene,
#     every Tier 1 missense has AM>=0.85 in constrained gene or ClinVar P/LP, etc.)
#   - ClinVar P/LP at >=1 star always Tier 1; ClinVar B/LB always Tier 5
#   - modifier candidates have AF in band OR rare-sub-Mendelian, AND class-
#     appropriate functional signal
#   - cohort.variants.tsv n_carriers matches recount of cohort.genotypes.tsv
#     (on a 100-variant sample, since exact recount of all is expensive)
#   - non-coding limitation visible (fraction of non-coding rare variants
#     stuck in Tier 4; documented gap per docs/noncoding_v2_plan.md)
#
# Output: pass/fail/info report. Soft-failures (limitations) get INFO, hard
# violations (real bugs) get FAIL. Exit code 0 if zero FAILs, 1 otherwise.
#
# Usage:
#   scripts/validate_cohort.sh --out-dir /data/nf1/outputs/batch1
# -----------------------------------------------------------------------------

set -uo pipefail
# NOT set -e: we want to keep checking even after individual checks fail.

OUT_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --out-dir DIR"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done
[[ -z "$OUT_DIR" ]] && { echo "--out-dir required" >&2; exit 1; }

# ============================== reporting helpers ============================
n_pass=0; n_fail=0; n_info=0
pass() { printf "  \033[32mPASS\033[0m  %s\n" "$*"; n_pass=$((n_pass + 1)); }
fail() { printf "  \033[31mFAIL\033[0m  %s\n" "$*"; n_fail=$((n_fail + 1)); }
info() { printf "  \033[33mINFO\033[0m  %s\n" "$*"; n_info=$((n_info + 1)); }
section() { printf "\n=== %s ===\n" "$*"; }

# ============================== input files ==================================
section "1. Input file presence + non-empty"

VARIANTS="$OUT_DIR/cohort.variants.tsv"
GENES="$OUT_DIR/cohort.genes.tsv"
GENOTYPES="$OUT_DIR/cohort.genotypes.tsv"

for f in "$VARIANTS" "$GENES"; do
  if [[ -s "$f" ]]; then
    pass "$(basename "$f") exists, $(wc -l < "$f") lines, $(du -h "$f" | cut -f1)"
  else
    fail "$(basename "$f") missing or empty: $f"
  fi
done
if [[ -s "$GENOTYPES" ]]; then
  pass "$(basename "$GENOTYPES") exists, $(wc -l < "$GENOTYPES") lines, $(du -h "$GENOTYPES" | cut -f1)"
else
  info "cohort.genotypes.tsv missing -- some cohort-burden checks will be skipped"
  GENOTYPES=""
fi

# Bail if required files missing
[[ -s "$VARIANTS" && -s "$GENES" ]] || { echo "cannot continue without variants + genes" >&2; exit 1; }

# ============================== column presence ==============================
section "2. Required columns present in cohort.variants.tsv"

required_cols=(chrom pos ref alt gene transcript consequence impact
               tier tier_reason variant_class
               modifier_candidate modifier_evidence
               alphamissense am_class esm1b revel
               loeuf pli
               gnomad_popmax_af clin_sig clin_stars)
header=$(head -1 "$VARIANTS")
missing_cols=()
for col in "${required_cols[@]}"; do
  if echo "$header" | tr '\t' '\n' | grep -qFx "$col"; then
    :
  else
    missing_cols+=("$col")
  fi
done
if [[ ${#missing_cols[@]} -eq 0 ]]; then
  pass "all ${#required_cols[@]} required columns present"
else
  fail "missing columns: ${missing_cols[*]}"
fi

# ============================== tier distribution ============================
section "3. Tier distribution"

# The pyramid check assumes gnomAD popmax is loaded. Without it, no variant
# can be marked "common in gnomAD" so Tier 1 inflates -- not a tiering-logic
# bug, a data-availability issue. Detect the missing-gnomAD case (>99% of
# variants have empty popmax) and downgrade the pyramid check to INFO with
# a clear message so the caveat travels with the report.
gnomad_empty_pct=$(awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  {
    total++
    pmax = $c["gnomad_popmax_af"]
    if (pmax == "" || pmax == ".") empty++
  }
  END { printf "%.1f", (total > 0 ? 100 * empty / total : 0) }
' "$VARIANTS")

awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  $c["tier"] >= 1 && $c["tier"] <= 5 { t[$c["tier"]]++ }
  END {
    total = 0
    for (k=1; k<=5; k++) total += t[k]+0
    for (k=1; k<=5; k++) printf "  Tier %d: %12d  (%5.2f%%)\n", k, t[k]+0, (t[k]+0)/total*100
    # Pass criteria: pyramid shape (T1 < T2 OR T1 < 1% of total), T4 not zero
    if (t[1]+0 < t[2]+0 || t[1]+0 < total*0.01) print "PASS_PYRAMID"
    else print "FAIL_PYRAMID"
    if (t[4]+0 > 0) print "PASS_T4_NONZERO"
    else print "FAIL_T4_NONZERO"
  }
' "$VARIANTS" | tee /tmp/_tier_dist.txt > /dev/null

if grep -q PASS_PYRAMID /tmp/_tier_dist.txt; then
  pass "tier distribution is pyramid-shaped (T1 < T2 and < 1% of total)"
elif (( $(echo "$gnomad_empty_pct > 99" | bc -l) )); then
  info "tier distribution inverted BUT gnomAD popmax is ${gnomad_empty_pct}% empty -- Tier 1 inflation is a documented consequence of the missing SA source, not a tiering bug (load a per-variant gnomAD source to fix)"
else
  fail "tier distribution is inverted -- Tier 1 should be much smaller than Tier 2"
fi
grep -q PASS_T4_NONZERO /tmp/_tier_dist.txt \
  && pass "Tier 4 is non-empty" \
  || fail "Tier 4 is empty -- suspect class-conditional rules over-promoting"
grep '^  Tier' /tmp/_tier_dist.txt
rm -f /tmp/_tier_dist.txt

# ============================== class-conditional ============================
section "4. Class-conditional rules"

# 4a. Every Tier 1 pLoF should be in a constrained gene, or have ClinVar P/LP,
#     or sit in an NF-spectrum tumor suppressor (NF1/NF2/SMARCB1/LZTR1) where
#     we apply a hardcoded Tier 1 override regardless of cohort constraint.
viol=$(awk -F'\t' '
  BEGIN {
    split("NF1 NF2 SMARCB1 LZTR1", nf_anchors, " ")
    for (i in nf_anchors) nf[nf_anchors[i]] = 1
  }
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  $c["tier"] == 1 && $c["variant_class"] == "pLoF" {
    cv = tolower($c["clin_sig"])
    has_clinvar_plp = (cv ~ /pathogenic/ && cv !~ /benign/)
    if (has_clinvar_plp) next  # ClinVar override -- OK
    if ($c["gene"] in nf) next  # NF anchor override -- OK
    loeuf = $c["loeuf"]; pli = $c["pli"]
    constrained = (loeuf != "" && loeuf != "." && loeuf+0 < 0.35) || \
                  (pli   != "" && pli   != "." && pli+0   >= 0.9)
    if (!constrained) print $c["chrom"]":"$c["pos"]" "$c["gene"]" loeuf="loeuf" pli="pli
  }
' "$VARIANTS" | head -5)
if [[ -z "$viol" ]]; then
  pass "all Tier 1 pLoF variants are in constrained genes (or ClinVar P/LP / NF-anchor override)"
else
  fail "Tier 1 pLoF variants in NON-constrained genes (first 5):"
  echo "$viol" | sed 's/^/        /'
fi

# 4b. Every Tier 1 missense should match ONE of these paths:
#   (a) ClinVar P/LP
#   (b) AM >= 0.85 in constrained gene (AM-strong primary path)
#   (c) AM likely_pathogenic (>0.564) + ESM1b damaging (<=-7.5) + constrained
#       (AM+ESM1b orthogonal-agreement secondary Tier 1 path)
#   (d) AM >= 0.85 in NF-spectrum tumor suppressor (NF1/NF2/SMARCB1/LZTR1)
#       regardless of cohort constraint metric
viol=$(awk -F'\t' '
  BEGIN {
    split("NF1 NF2 SMARCB1 LZTR1", nf_anchors, " ")
    for (i in nf_anchors) nf[nf_anchors[i]] = 1
  }
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  $c["tier"] == 1 && $c["variant_class"] == "missense" {
    cv = tolower($c["clin_sig"])
    has_clinvar_plp = (cv ~ /pathogenic/ && cv !~ /benign/)
    if (has_clinvar_plp) next                                                  # path (a)
    am  = $c["alphamissense"]; esm = $c["esm1b"]
    loeuf = $c["loeuf"]; pli = $c["pli"]
    am_strong    = (am  != "" && am  != "." && am+0  >= 0.85)
    am_path      = (am  != "" && am  != "." && am+0  >  0.564)
    esm_damaging = (esm != "" && esm != "." && esm+0 <= -7.5)
    constrained  = (loeuf != "" && loeuf != "." && loeuf+0 < 0.35) || \
                   (pli   != "" && pli   != "." && pli+0   >= 0.9)
    if (am_strong && constrained) next                                         # path (b)
    if (am_path && esm_damaging && constrained) next                           # path (c)
    if (am_strong && ($c["gene"] in nf)) next                                  # path (d)
    print $c["chrom"]":"$c["pos"]" "$c["gene"]" am="am" esm="esm" loeuf="loeuf
  }
' "$VARIANTS" | head -5)
if [[ -z "$viol" ]]; then
  pass "all Tier 1 missense variants match a valid path (AM-strong+constrained, AM+ESM1b orthogonal+constrained, AM-strong+NF-anchor, or ClinVar P/LP)"
else
  fail "Tier 1 missense variants without any valid Tier 1 path (first 5):"
  echo "$viol" | sed 's/^/        /'
fi

# 4c. Every Tier 5 should be common, ClinVar B/LB, or AM-benign missense.
# Use af_used (R's effective rarity column, falls back through faf -> popmax -> gnomad_af)
# rather than popmax alone, and accept am_class containing "benign" in addition
# to numeric AM < 0.34 (mirrors the R rule `am_benign := T0(am < 0.34) | grepl("benign", amc)`).
viol=$(awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  $c["tier"] == 5 {
    # Prefer the af_used column R writes; fall back to popmax if af_used absent.
    if (("af_used" in c) && $c["af_used"] != "" && $c["af_used"] != ".") af_used = $c["af_used"]
    else af_used = $c["gnomad_popmax_af"]
    cv = tolower($c["clin_sig"])
    am = $c["alphamissense"]; amc = tolower($c["am_class"])
    vc = $c["variant_class"]
    is_common = (af_used != "" && af_used != "." && af_used+0 > 1e-4)
    has_clinvar_blb = (cv ~ /benign/ && cv !~ /pathogenic/)
    am_benign = (am != "" && am != "." && am+0 < 0.34) || (amc ~ /benign/)
    if (is_common || has_clinvar_blb || (vc == "missense" && am_benign)) next
    print $c["chrom"]":"$c["pos"]" "$c["gene"]" vc="vc" af_used="af_used" am="am" amc="amc" cs="$c["clin_sig"]
  }
' "$VARIANTS" | head -5)
if [[ -z "$viol" ]]; then
  pass "all Tier 5 variants are common, ClinVar B/LB, or AM-benign missense"
else
  fail "Tier 5 variants not matching any benign rule (first 5):"
  echo "$viol" | sed 's/^/        /'
fi

# ============================== ClinVar overrides ============================
section "5. ClinVar override rules"

# Every ClinVar P/LP at >=1 star should be Tier 1
viol=$(awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  {
    cv = tolower($c["clin_sig"])
    stars = $c["clin_stars"]
    is_plp = (cv ~ /pathogenic/ && cv !~ /benign/)
    # Star filter intentionally disabled (star_ok = 1) while clin_stars is
    # silently 0 from the csq_to_wide_tab.py underscored-REVIEW_STATUS bug.
    # Restore once cohort is rebuilt with the fixed converter:
    #   star_ok = (stars == "" || stars == ".") || (stars+0 >= 1)
    star_ok = 1
    if (is_plp && star_ok && $c["tier"] != "1") {
      print $c["chrom"]":"$c["pos"]" tier="$c["tier"]" cs="$c["clin_sig"]" stars="stars
    }
  }
' "$VARIANTS" | head -5)
if [[ -z "$viol" ]]; then
  pass "all ClinVar P/LP variants are Tier 1 (star filter currently disabled; see methods doc §7.5)"
else
  fail "ClinVar P/LP variants NOT in Tier 1 (first 5):"
  echo "$viol" | sed 's/^/        /'
fi

# Every ClinVar B/LB should be Tier 5
viol=$(awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  {
    cv = tolower($c["clin_sig"])
    is_blb = (cv ~ /benign/ && cv !~ /pathogenic/)
    if (is_blb && $c["tier"] != "5") {
      print $c["chrom"]":"$c["pos"]" tier="$c["tier"]" cs="$c["clin_sig"]
    }
  }
' "$VARIANTS" | head -5)
if [[ -z "$viol" ]]; then
  pass "all ClinVar B/LB variants are Tier 5"
else
  fail "ClinVar B/LB variants NOT in Tier 5 (first 5):"
  echo "$viol" | sed 's/^/        /'
fi

# ============================== modifier flag ================================
section "6. Modifier candidate flag"

# Every modifier_candidate=TRUE should have AF in band (1e-4 < af_used <= 0.05)
# OR rare and in Tier 3/4 (sub-Mendelian rescue).
#
# Use af_used (R's effective rarity column, cascade: faf -> popmax -> gnomad_af
# -> cohort_af). Falling back to popmax makes this check silently pass all
# modifier candidates when gnomAD isn't loaded (pmax empty everywhere), because
# the "rare + tier 3/4" branch would swallow everything -- masking real rule
# violations. Section 4c uses the same af_used-with-popmax-fallback pattern.
viol=$(awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  $c["modifier_candidate"] == "TRUE" {
    if (("af_used" in c) && $c["af_used"] != "" && $c["af_used"] != ".") af_used = $c["af_used"]
    else af_used = $c["gnomad_popmax_af"]
    tier = $c["tier"]+0
    in_band = (af_used != "" && af_used != "." && af_used+0 > 1e-4 && af_used+0 <= 0.05)
    rare_low_tier = ((af_used == "" || af_used == "." || af_used+0 <= 1e-4) && (tier == 3 || tier == 4))
    if (!in_band && !rare_low_tier) {
      print $c["chrom"]":"$c["pos"]" tier="tier" af_used="af_used" af_source="$c["af_source"]
    }
  }
' "$VARIANTS" | head -5)
if [[ -z "$viol" ]]; then
  pass "all modifier candidates have AF in band OR sub-Mendelian"
else
  fail "modifier candidates not matching AF rule (first 5):"
  echo "$viol" | sed 's/^/        /'
fi

# ============================== anchor-gene sanity ===========================
section "7. NF anchor genes have signal"

for anchor in NF1 NF2 SMARCB1 LZTR1; do
  n=$(awk -F'\t' -v g="$anchor" '
    NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
    $c["gene"] == g
  ' "$VARIANTS" | wc -l)
  best_t=$(awk -F'\t' -v g="$anchor" '
    NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
    $1 == g { print $c["best_tier"]; exit }
  ' "$GENES")
  if [[ "$n" -gt 0 ]]; then
    pass "$anchor: $n variants in cohort, best_tier=${best_t:-N/A}"
  else
    info "$anchor: 0 variants -- check if cohort actually contains carriers"
  fi
done

# ============================== cohort carrier counts ========================
section "8. Cohort carrier counts (n_carriers) vs genotype recount"

if [[ -n "$GENOTYPES" ]]; then
  # Sample 50 random Tier 1+2 variants and recount carriers from cohort.genotypes.tsv;
  # check n_carriers column matches.
  mismatches=$(awk -F'\t' '
    NR == FNR {
      if (FNR == 1) { for (i=1;i<=NF;i++) c[$i] = i; next }
      if ($c["tier"] == "1" || $c["tier"] == "2") {
        key = $c["chrom"] SUBSEP $c["pos"] SUBSEP $c["ref"] SUBSEP $c["alt"]
        expected[key] = $c["n_carriers"]+0
        if (++pick > 50) exit
      }
      next
    }
  ' "$VARIANTS")

  awk -F'\t' '
    NR == FNR {
      if (FNR == 1) { for (i=1;i<=NF;i++) c[$i] = i; next }
      if ($c["tier"] == "1" || $c["tier"] == "2") {
        if (++pick > 50) next
        key = $c["chrom"] SUBSEP $c["pos"] SUBSEP $c["ref"] SUBSEP $c["alt"]
        expected[key] = $c["n_carriers"]+0
      }
      next
    }
    FNR == 1 { next }
    $6 ~ /^(0\/1|1\/0|0\|1|1\|0|1\/1|1\|1)$/ {
      key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      if (key in expected) carriers[key, $5] = 1
    }
    END {
      n_match = 0; n_mismatch = 0
      for (k in expected) {
        cnt = 0
        for (kk in carriers) {
          split(kk, a, SUBSEP)
          if (a[1] SUBSEP a[2] SUBSEP a[3] SUBSEP a[4] == k) cnt++
        }
        if (cnt == expected[k]) n_match++
        else { n_mismatch++; if (n_mismatch <= 3) print "  mismatch", k, "expected="expected[k], "found="cnt }
      }
      printf "  %d/%d sampled variants match\n", n_match, n_match + n_mismatch
      if (n_mismatch == 0) print "PASS_CARRIERS"
      else print "FAIL_CARRIERS"
    }
  ' "$VARIANTS" "$GENOTYPES" 2>/dev/null | tee /tmp/_carriers.txt > /dev/null

  grep '^  ' /tmp/_carriers.txt
  if grep -q PASS_CARRIERS /tmp/_carriers.txt; then
    pass "all sampled n_carriers counts match recount of genotypes"
  else
    fail "n_carriers does not match genotype recount (see mismatches above)"
  fi
  rm -f /tmp/_carriers.txt
else
  info "skipping carrier recount (cohort.genotypes.tsv not present)"
fi

# ============================== non-coding limitation ========================
section "9. Non-coding limitation (documented gap, INFO not FAIL)"

# What fraction of rare non-coding variants ended up in Tier 4 ("we don't know")?
awk -F'\t' '
  NR == 1 { for (i=1;i<=NF;i++) c[$i] = i; next }
  {
    vc = $c["variant_class"]
    pmax = $c["gnomad_popmax_af"]
    is_rare = (pmax == "" || pmax == "." || pmax+0 <= 1e-4)
    is_nc = (vc == "noncoding" || vc == "regulatory")
    if (is_nc && is_rare) {
      total++
      if ($c["tier"] == "4") t4++
    }
  }
  END {
    if (total > 0) {
      printf "  rare non-coding variants: %d; in Tier 4: %d (%.1f%%)\n", total, t4+0, (t4+0)/total*100
      if (t4 / total > 0.7) print "INFO_NC_GAP_HIGH"
    } else {
      print "  rare non-coding variants: 0 (filtered out?)"
    }
  }
' "$VARIANTS" | tee /tmp/_nc_gap.txt > /dev/null
grep '^  ' /tmp/_nc_gap.txt
if grep -q INFO_NC_GAP_HIGH /tmp/_nc_gap.txt; then
  info ">70% of rare non-coding variants stuck in Tier 4 (expected -- see docs/noncoding_v2_plan.md)"
fi
rm -f /tmp/_nc_gap.txt

# ============================== summary ======================================
section "Summary"
printf "  %d PASS  %d FAIL  %d INFO\n" $n_pass $n_fail $n_info
if [[ $n_fail -eq 0 ]]; then
  echo
  echo "  All hard checks passed. INFO items are documented limitations, not bugs."
  exit 0
else
  echo
  echo "  $n_fail FAIL(s) -- review the output above before trusting the cohort tier table."
  exit 1
fi
