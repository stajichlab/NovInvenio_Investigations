#!/usr/bin/env python3
"""Extract per-protein GO/Pfam/InterPro/gene-name annotation from a UniProt .dat.gz.

This is the "IPRScan summary" mechanism decided in DESIGN.md Sec 5/7: every protein in
a UniProt reference proteome already carries precomputed GO/Pfam/InterPro cross-
references (DR lines) -- no fresh InterProScan run needed for v1. This script surfaces
them as one flat per-accession TSV, which both the ORA enrichment scripts (Sec 7) and
a future per-candidate report summary consume.

Deliberately NOT using Biopython's SwissProt parser: that parses the entire record
(references, comments, sequence, everything) to get at a handful of line types. This
reads only AC/GN/OX/DR lines and resets on the `//` record separator -- a targeted,
single-pass parser, addressing the parser-performance question flagged as deferred in
DESIGN.md Sec 5/9. (If this ever turns out not fast enough for larger proteomes, the
fallback is still Biopython -- not the other direction.)

Output columns: accession, taxon_id, gene_name, description, go_ids, pfam_ids,
pfam_names, interpro_ids, ec_numbers, alphafold_id
  - go_ids / pfam_ids / pfam_names / interpro_ids / ec_numbers are "|"-separated
    (same order across pfam_ids/pfam_names -- index i of each names the same
    domain); each GO entry carries its evidence code as "GO:0016020:IEA"
    (evidence matters for ORA -- see Sec 10 addendum on TrEMBL GO being mostly IEA).
  - pfam_names is the DR Pfam line's short mnemonic name (e.g. "RmlD_sub_bind"),
    not a long description -- what nf_NovInvenio's own Pfam_Names column holds
    for an ANNOTATE_MATRIX-produced hit, so this stays display-compatible with it.
  - description is the DE line's RecName (falling back to SubName when a protein has
    no reviewed RecName, e.g. most TrEMBL entries) Full= value, semicolon/EC-number
    stripped -- the same "what is this protein" text a report's product_description
    field is meant to hold.
  - ec_numbers comes only from DE lines (the protein-name block's own EC=... value),
    not from CC comment lines that mention an EC number in prose (catalytic-activity
    descriptions among ChEBI references) -- those are far less reliable to anchor a
    regex on and are not necessarily *this* protein's own EC assignment.
  - alphafold_id is the DR AlphaFoldDB cross-reference (usually, but not assumed to
    always be, identical to `accession`) -- AlphaFold DB covers nearly all of
    UniProt, so this fires for nearly every protein and replaces the report's
    generic "structure search" fallback with a direct predicted-structure link.
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

AC_RE = re.compile(r"^AC\s+([A-Z0-9]+)")
OX_RE = re.compile(r"NCBI_TaxID=(\d+)")
GN_RE = re.compile(r"(?:Name|ORFNames)=([^;{]+)")
DE_RECNAME_RE = re.compile(r"^DE\s+RecName:\s*Full=([^;{]+)")
DE_SUBNAME_RE = re.compile(r"^DE\s+SubName:\s*Full=([^;{]+)")
DR_GO_RE = re.compile(r"^DR\s+GO;\s*(GO:\d+);\s*[A-Z]:[^;]*;\s*([A-Za-z0-9_]+):")
# Captures both the accession and the short mnemonic name (3rd field), e.g.
# "DR   Pfam; PF04321; RmlD_sub_bind; 1." -- the name is what a report should
# display; the bare accession alone (PFxxxxx) is functional but not friendly.
DR_PFAM_RE = re.compile(r"^DR\s+Pfam;\s*(PF\d+);\s*([^;]+);")
DR_INTERPRO_RE = re.compile(r"^DR\s+InterPro;\s*(IPR\d+);")
DR_ALPHAFOLD_RE = re.compile(r"^DR\s+AlphaFoldDB;\s*([A-Z0-9]+);")
# EC numbers only from DE lines (the protein-name block, e.g.
# "DE   RecName: Full=...; EC=2.7.7.7;" or a continuation "DE   EC=2.7.7.7;") --
# not from CC comment lines, which mention EC numbers in prose (catalytic-activity
# descriptions with ChEBI references) that are much less reliable to anchor on.
DE_EC_RE = re.compile(r"EC=([\d.]+(?:-)?)")


def parse_dat_gz(path: Path):
    """Yield one dict per protein entry."""
    accession = None
    taxon_id = None
    gene_name = None
    rec_description = None
    sub_description = None
    go_ids = []
    pfam_ids = []
    pfam_names = []
    interpro_ids = []
    ec_numbers = []
    alphafold_id = None

    def reset():
        nonlocal accession, taxon_id, gene_name, rec_description, sub_description
        nonlocal go_ids, pfam_ids, pfam_names, interpro_ids, ec_numbers, alphafold_id
        accession = taxon_id = gene_name = rec_description = sub_description = None
        go_ids, pfam_ids, pfam_names, interpro_ids, ec_numbers = [], [], [], [], []
        alphafold_id = None

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("//"):
                if accession:
                    yield {
                        "accession": accession,
                        "taxon_id": taxon_id or "",
                        "gene_name": gene_name or "",
                        "description": (rec_description or sub_description or "").strip(),
                        "go_ids": "|".join(go_ids),
                        "pfam_ids": "|".join(pfam_ids),
                        "pfam_names": "|".join(pfam_names),
                        "interpro_ids": "|".join(interpro_ids),
                        "ec_numbers": "|".join(ec_numbers),
                        "alphafold_id": alphafold_id or "",
                    }
                reset()
                continue
            if accession is None and line.startswith("AC"):
                m = AC_RE.match(line)
                if m:
                    accession = m.group(1)
                continue
            if gene_name is None and line.startswith("GN"):
                m = GN_RE.search(line)
                if m:
                    gene_name = m.group(1).strip()
                continue
            if line.startswith("DE"):
                if rec_description is None:
                    m = DE_RECNAME_RE.match(line)
                    if m:
                        rec_description = m.group(1).strip()
                if sub_description is None:
                    m = DE_SUBNAME_RE.match(line)
                    if m:
                        sub_description = m.group(1).strip()
                m = DE_EC_RE.search(line)
                if m:
                    ec_numbers.append(m.group(1))
                continue
            if taxon_id is None and line.startswith("OX"):
                m = OX_RE.search(line)
                if m:
                    taxon_id = m.group(1)
                continue
            if line.startswith("DR"):
                m = DR_GO_RE.match(line)
                if m:
                    go_ids.append(f"{m.group(1)}:{m.group(2)}")
                    continue
                m = DR_PFAM_RE.match(line)
                if m:
                    pfam_ids.append(m.group(1))
                    pfam_names.append(m.group(2).strip())
                    continue
                m = DR_INTERPRO_RE.match(line)
                if m:
                    interpro_ids.append(m.group(1))
                    continue
                if alphafold_id is None:
                    m = DR_ALPHAFOLD_RE.match(line)
                    if m:
                        alphafold_id = m.group(1)
                        continue


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dat-gz", required=True, type=Path, help="Path to a UniProt {Proteome}_{taxid}.dat.gz")
    ap.add_argument("--output", required=True, type=Path, help="Output TSV path")
    args = ap.parse_args()

    n = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "taxon_id", "gene_name", "description", "go_ids", "pfam_ids", "pfam_names", "interpro_ids", "ec_numbers", "alphafold_id"], delimiter="\t")
        w.writeheader()
        for rec in parse_dat_gz(args.dat_gz):
            w.writerow(rec)
            n += 1

    print(f"Wrote {args.output} ({n} proteins)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
