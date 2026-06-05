#!/usr/bin/env python3
"""
Convert a fastvep-annotated VCF into a wide-format tab table for tier_variants.R.

Reads the CSQ INFO field (per-transcript consequence + 48 VEP fields), each
FV_* projection (FV_DBNSFP, FV_GNOMAD, FV_CLINVAR, FV_OMIM, FV_GNOMAD_GENE,
SpliceAI), splits all the pipe-delimited values into named columns, and emits
one TSV row per (variant, transcript) consequence -- exactly the column set
that R/tier_variants.R's COLS contract expects.

Stdlib only. Reads gz or plain VCF. Handles gracefully when an FV_* source
isn't loaded (its columns just stay blank).

Usage:
    python3 csq_to_wide_tab.py --input annotated.vcf.gz --output annotated.tab [--sample S001]
"""
import argparse
import gzip
import io
import sys

# Output columns in order. Matches R/tier_variants.R's COLS contract names.
OUT_COLS = [
    "chrom", "pos", "ref", "alt",
    "sample",
    "gene", "gene_id", "transcript",
    "mane_select", "canonical", "mane_plus_clinical",
    "biotype", "consequence", "impact",
    "exon", "intron",
    "hgvsg", "hgvsc", "hgvsp",
    "existing_variation", "domains", "variant_class",
    "pheno", "gene_pheno", "motif_name",
    "loftee_lof",
    "gnomad_af", "gnomad_faf", "gnomad_popmax_af",
    "clin_sig", "clin_stars",
    "loeuf", "pli", "mis_z", "syn_z",
    "omim_phenotype",
    "spliceai_ds_max",
    "alphamissense", "am_class", "esm1b", "esm1b_class", "revel",
    "sift", "polyphen",
    "nmd_escape",
]

# Map fastvep CSQ field names -> our output column names. Anything not listed
# here is read from CSQ but not emitted; anything mapped to None means "carry
# its value but to a column derived in special logic (e.g. FLAGS -> canonical)".
CSQ_TO_OUT = {
    "Consequence": "consequence",
    "IMPACT": "impact",
    "SYMBOL": "gene",
    "Gene": "gene_id",
    "Feature": "transcript",
    "BIOTYPE": "biotype",
    "EXON": "exon",
    "INTRON": "intron",
    "HGVSc": "hgvsc",
    "HGVSp": "hgvsp",
    "Existing_variation": "existing_variation",
    "CANONICAL": "canonical",
    "MANE_SELECT": "mane_select",
    "MANE_PLUS_CLINICAL": "mane_plus_clinical",
    "SIFT": "sift",
    "PolyPhen": "polyphen",
    "CLIN_SIG": "clin_sig",
    "PHENO": "pheno",
    "MOTIF_NAME": "motif_name",
}

GNOMAD_POP_KEYS = [
    "AFR_AF", "AMR_AF", "ASJ_AF", "EAS_AF", "FIN_AF",
    "MID_AF", "NFE_AF", "OTH_AF", "REMAINING_AF", "SAS_AF",
]


# --------------------------------------------------------------------------- helpers
def open_vcf(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_format_field(line):
    """Extract the 'Format: F1|F2|...' field list from a ##INFO=<...> header."""
    idx = line.find("Format: ")
    if idx < 0:
        return []
    fmt = line[idx + len("Format: "):]
    # strip trailing ">\n and stray quotes
    fmt = fmt.rstrip("\n").rstrip(">").rstrip('"')
    return fmt.split("|")


def parse_info(info_str):
    """Parse VCF INFO column into dict (flags become True)."""
    d = {}
    for kv in info_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def safe_float(v):
    if v in ("", "-", ".", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def map_review_status_to_stars(rs):
    """Coarse mapping from ClinVar REVIEW_STATUS to ClinVar gold-star count."""
    if not rs:
        return ""
    rs = rs.lower()
    if "practice guideline" in rs:
        return "4"
    if "reviewed by expert panel" in rs:
        return "3"
    if "criteria provided, multiple submitters" in rs:
        return "2"
    if "criteria provided" in rs:
        return "1"
    return "0"


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--input", required=True, help="fastvep-annotated VCF (.vcf, .vcf.gz, .vcf.bgz)")
    ap.add_argument("--output", required=True, help="output TSV path")
    ap.add_argument("--sample", default="",
                    help="value to write into the 'sample' column (default empty)")
    args = ap.parse_args()

    csq_fields = []
    # info_id -> [field names] (ALLELE stripped)
    fv_fields = {}

    n_rows = 0
    with open_vcf(args.input) as fh, open(args.output, "w") as out:
        out.write("\t".join(OUT_COLS) + "\n")

        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ,"):
                    csq_fields = parse_format_field(line)
                elif line.startswith("##INFO=<ID=FV_") or line.startswith("##INFO=<ID=SpliceAI,"):
                    # extract the ID
                    head = line.split(",", 1)[0]
                    info_id = head.split("=")[-1]
                    fmt = parse_format_field(line)
                    if fmt:
                        # First entry is ALLELE for per-allele projections; skip it.
                        if fmt[0].upper() == "ALLELE":
                            fv_fields[info_id] = fmt[1:]
                        else:
                            fv_fields[info_id] = fmt
                continue
            if line.startswith("#"):
                continue
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _, ref, _alts, _, _, info_str = fields[:8]
            info = parse_info(info_str)

            csq_val = info.get("CSQ", "")
            if not csq_val:
                continue
            csq_entries = csq_val.split(",")

            # Build per-allele lookup for each FV_* projection.
            fv_by_allele = {}  # info_id -> {allele -> dict}
            for fid, fnames in fv_fields.items():
                raw = info.get(fid, "")
                if not raw:
                    continue
                d = {}
                for entry in raw.split(","):
                    parts = entry.split("|")
                    if not parts:
                        continue
                    allele = parts[0]
                    rest = parts[1:]
                    d[allele] = {
                        f: (rest[i] if i < len(rest) else "")
                        for i, f in enumerate(fnames)
                    }
                fv_by_allele[fid] = d

            for csq_entry in csq_entries:
                csq_parts = csq_entry.split("|")
                csq = {
                    f: (csq_parts[i] if i < len(csq_parts) else "")
                    for i, f in enumerate(csq_fields)
                }
                allele = csq.get("Allele", "")
                if not allele:
                    continue

                row = {col: "" for col in OUT_COLS}
                row["chrom"] = chrom
                row["pos"] = pos
                row["ref"] = ref
                row["alt"] = allele
                row["sample"] = args.sample

                # CSQ -> our columns (direct mapping).
                for csq_name, out_name in CSQ_TO_OUT.items():
                    if out_name:
                        row[out_name] = csq.get(csq_name, "")

                # FLAGS may carry "canonical"; CANONICAL field carries YES/NO too.
                if not row["canonical"]:
                    flags = csq.get("FLAGS", "").lower()
                    if "canonical" in flags:
                        row["canonical"] = "TRUE"
                else:
                    row["canonical"] = "TRUE" if row["canonical"].upper() in ("YES", "TRUE", "1") else "FALSE"

                # MANE_SELECT is "ENSTxxx" if present (truthy), blank if not. Normalize to TRUE/FALSE.
                row["mane_select"] = "TRUE" if row["mane_select"] else "FALSE"
                row["mane_plus_clinical"] = "TRUE" if row["mane_plus_clinical"] else "FALSE"

                # --- FV_DBNSFP: per-allele score block ---
                d = fv_by_allele.get("FV_DBNSFP", {}).get(allele)
                if d:
                    # NB: SIFT/POLYPHEN in dbNSFP are different from CSQ SIFT/PolyPhen;
                    # keep CSQ values (already in row) and overlay dbNSFP only for the
                    # five new fields we patched in.
                    row["alphamissense"] = d.get("ALPHAMISSENSE", "")
                    row["am_class"] = d.get("AM_CLASS", "")
                    row["esm1b"] = d.get("ESM1B", "")
                    row["esm1b_class"] = d.get("ESM1B_CLASS", "")
                    row["revel"] = d.get("REVEL", "")

                # --- SpliceAI: ds_max = max(DS_AG, DS_AL, DS_DG, DS_DL) ---
                d = fv_by_allele.get("SpliceAI", {}).get(allele)
                if d:
                    vals = [safe_float(d.get(k)) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")]
                    vals = [v for v in vals if v is not None]
                    if vals:
                        row["spliceai_ds_max"] = f"{max(vals):.4f}"

                # --- FV_GNOMAD: ALL_AF and popmax ---
                d = fv_by_allele.get("FV_GNOMAD", {}).get(allele)
                if d:
                    row["gnomad_af"] = d.get("ALL_AF", "")
                    pop_afs = [safe_float(d.get(k)) for k in GNOMAD_POP_KEYS]
                    pop_afs = [v for v in pop_afs if v is not None]
                    if pop_afs:
                        row["gnomad_popmax_af"] = f"{max(pop_afs):.6g}"

                # --- FV_CLINVAR: SIGNIFICANCE + stars from REVIEW_STATUS ---
                d = fv_by_allele.get("FV_CLINVAR", {}).get(allele)
                if d:
                    row["clin_sig"] = d.get("SIGNIFICANCE", "")
                    row["clin_stars"] = map_review_status_to_stars(d.get("REVIEW_STATUS", ""))

                # --- FV_GNOMAD_GENE: gene-level constraint (matches first entry for this allele) ---
                d = fv_by_allele.get("FV_GNOMAD_GENE", {}).get(allele)
                if d:
                    row["loeuf"] = d.get("LOEUF", "")
                    row["pli"] = d.get("PLI", "")
                    row["mis_z"] = d.get("MIS_Z", "")
                    row["syn_z"] = d.get("SYN_Z", "")

                # --- FV_OMIM: gene-level phenotype text ---
                d = fv_by_allele.get("FV_OMIM", {}).get(allele)
                if d:
                    row["omim_phenotype"] = d.get("PHENOTYPES", "")

                out.write("\t".join(row[c] for c in OUT_COLS) + "\n")
                n_rows += 1

    print(f"[csq_to_wide_tab] wrote {n_rows:,} rows -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
