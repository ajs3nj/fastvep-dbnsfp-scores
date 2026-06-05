#!/usr/bin/env python3
"""Extract AlphaMissense / ESM1b / REVEL from a dbNSFP v4.5+ TSV.

CADD is intentionally not pulled; see docs/tiering.md §1.2.

Standalone (Python 3, stdlib only). Two uses:
  1. Get the scores today, without touching fastVEP's Rust, as a tidy per-(variant, transcript) table
     you can fread() in R and join to your VEP output.
  2. A validation oracle: run this over the same dbNSFP build your fastVEP `--source dbnsfp` reads,
     and diff the numbers to confirm the Rust path agrees.

Columns are resolved by HEADER NAME (version-safe). Missing token is '.'. Per-transcript columns are
';'-split and aligned to Ensembl_transcriptid; single-cardinality columns (REVEL, CADD) are broadcast.

Usage:
    # whole file
    python3 extract_dbnsfp_scores.py dbNSFP4.7a_variant.complete.gz > scores.tsv

    # one region (requires the .gz to be tabix-indexed and `tabix` on PATH)
    python3 extract_dbnsfp_scores.py dbNSFP4.7a.gz --region 17:43044295-43125483 > brca1.tsv

    # only rows with an AlphaMissense or ESM1b score (skip non-missense noise)
    python3 extract_dbnsfp_scores.py dbNSFP4.7a.gz --require-protein-score > missense_scores.tsv
"""
import argparse
import gzip
import io
import subprocess
import sys

MISSING = "."

# dbNSFP column names (v4.5+). Confirm against your build's .readme.
COL_CHR = "#chr"
COL_POS = "pos(1-based)"
COL_REF = "ref"
COL_ALT = "alt"
COL_AAREF = "aaref"
COL_AAALT = "aaalt"
COL_TX = "Ensembl_transcriptid"
COL_PROT = "Ensembl_proteinid"
COL_GENE = "genename"
COL_REVEL = "REVEL_score"
COL_AM = "AlphaMissense_score"
COL_AM_PRED = "AlphaMissense_pred"
COL_ESM = "ESM1b_score"
COL_ESM_PRED = "ESM1b_pred"

OUT_HEADER = [
    "chr", "pos", "ref", "alt", "aaref", "aaalt", "gene",
    "transcript", "protein",
    "revel",
    "alphamissense", "alphamissense_pred", "esm1b", "esm1b_pred",
]


def clean(cell):
    cell = cell.strip()
    return None if cell == "" or cell == MISSING else cell


def split_multi(cell):
    c = clean(cell)
    return [] if c is None else [clean(e) for e in c.split(";")]


def pick(values, i):
    """Element i, with broadcast of a single-element list and graceful out-of-range."""
    if len(values) == 1:
        return values[0]
    return values[i] if i < len(values) else None


def open_stream(path, region):
    if region:
        # tabix prints data rows only (no header); we read the header separately below.
        proc = subprocess.Popen(["tabix", path, region], stdout=subprocess.PIPE, text=True)
        return proc.stdout, proc
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8"), None
    return open(path, "r", encoding="utf-8"), None


def read_header(path):
    """Read the dbNSFP header line (first line starting with '#chr')."""
    opener = gzip.open if (path.endswith(".gz") or path.endswith(".bgz")) else open
    with (io.TextIOWrapper(opener(path, "rb"), encoding="utf-8")
          if opener is gzip.open else open(path, "r", encoding="utf-8")) as fh:
        for line in fh:
            if line.startswith(COL_CHR) or line.startswith("#chr"):
                return line.rstrip("\n\r").split("\t")
    raise SystemExit("error: could not find dbNSFP header line starting with '#chr'")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dbnsfp", help="dbNSFP TSV (.gz/.bgz ok)")
    ap.add_argument("--region", help="tabix region, e.g. 17:43044295-43125483 (needs .tbi + tabix)")
    ap.add_argument("--require-protein-score", action="store_true",
                    help="only emit transcript rows that have an AlphaMissense or ESM1b score")
    args = ap.parse_args()

    header = read_header(args.dbnsfp)
    col = {name: i for i, name in enumerate(header)}

    def cidx(name):
        return col.get(name)

    needed = [COL_CHR, COL_POS, COL_REF, COL_ALT, COL_TX]
    for n in needed:
        if n not in col:
            raise SystemExit(f"error: required column '{n}' not found in header; check dbNSFP version")

    have = {k: (k in col) for k in
            [COL_AM, COL_AM_PRED, COL_ESM, COL_ESM_PRED, COL_REVEL,
             COL_PROT, COL_GENE, COL_AAREF, COL_AAALT]}
    missing_cols = [k for k, v in have.items() if not v]
    if missing_cols:
        print(f"# note: columns absent in this dbNSFP build (emitted as '.'): {missing_cols}",
              file=sys.stderr)

    stream, proc = open_stream(args.dbnsfp, args.region)
    out = sys.stdout
    out.write("\t".join(OUT_HEADER) + "\n")

    def field(parts, name):
        i = cidx(name)
        return parts[i] if (i is not None and i < len(parts)) else MISSING

    try:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")

            chrom = field(parts, COL_CHR)
            pos = field(parts, COL_POS)
            ref = field(parts, COL_REF)
            alt = field(parts, COL_ALT)
            aaref = field(parts, COL_AAREF)
            aaalt = field(parts, COL_AAALT)

            revel = clean(field(parts, COL_REVEL))

            tx = split_multi(field(parts, COL_TX))
            prot = split_multi(field(parts, COL_PROT))
            genes = split_multi(field(parts, COL_GENE))
            am = split_multi(field(parts, COL_AM))
            am_pred = split_multi(field(parts, COL_AM_PRED))
            esm = split_multi(field(parts, COL_ESM))
            esm_pred = split_multi(field(parts, COL_ESM_PRED))

            n = max(len(x) for x in (tx, prot, am, am_pred, esm, esm_pred, [None])) or 1

            for i in range(n):
                am_i = pick(am, i)
                esm_i = pick(esm, i)
                if args.require_protein_score and am_i is None and esm_i is None:
                    continue
                row = [
                    chrom, pos, ref, alt,
                    aaref if aaref else MISSING,
                    aaalt if aaalt else MISSING,
                    pick(genes, i) or MISSING,
                    pick(tx, i) or MISSING,
                    pick(prot, i) or MISSING,
                    revel or MISSING,
                    am_i or MISSING,
                    pick(am_pred, i) or MISSING,
                    esm_i or MISSING,
                    pick(esm_pred, i) or MISSING,
                ]
                out.write("\t".join(row) + "\n")
    finally:
        try:
            stream.close()
        except Exception:
            pass
        if proc is not None:
            proc.wait()


if __name__ == "__main__":
    main()
