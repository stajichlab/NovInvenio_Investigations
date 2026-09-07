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

Output columns: accession, taxon_id, gene_name, go_ids, pfam_ids, interpro_ids
  - go_ids / pfam_ids / interpro_ids are "|"-separated; each GO entry carries its
    evidence code as "GO:0016020:IEA" (evidence matters for ORA -- see Sec 10 addendum
    on TrEMBL GO being mostly IEA).
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
DR_GO_RE = re.compile(r"^DR\s+GO;\s*(GO:\d+);\s*[A-Z]:[^;]*;\s*([A-Za-z0-9_]+):")
DR_PFAM_RE = re.compile(r"^DR\s+Pfam;\s*(PF\d+);")
DR_INTERPRO_RE = re.compile(r"^DR\s+InterPro;\s*(IPR\d+);")


def parse_dat_gz(path: Path):
    """Yield one dict per protein entry."""
    accession = None
    taxon_id = None
    gene_name = None
    go_ids = []
    pfam_ids = []
    interpro_ids = []

    def reset():
        nonlocal accession, taxon_id, gene_name, go_ids, pfam_ids, interpro_ids
        accession = taxon_id = gene_name = None
        go_ids, pfam_ids, interpro_ids = [], [], []

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("//"):
                if accession:
                    yield {
                        "accession": accession,
                        "taxon_id": taxon_id or "",
                        "gene_name": gene_name or "",
                        "go_ids": "|".join(go_ids),
                        "pfam_ids": "|".join(pfam_ids),
                        "interpro_ids": "|".join(interpro_ids),
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
                    continue
                m = DR_INTERPRO_RE.match(line)
                if m:
                    interpro_ids.append(m.group(1))
                    continue


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dat-gz", required=True, type=Path, help="Path to a UniProt {Proteome}_{taxid}.dat.gz")
    ap.add_argument("--output", required=True, type=Path, help="Output TSV path")
    args = ap.parse_args()

    n = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "taxon_id", "gene_name", "go_ids", "pfam_ids", "interpro_ids"], delimiter="\t")
        w.writeheader()
        for rec in parse_dat_gz(args.dat_gz):
            w.writerow(rec)
            n += 1

    print(f"Wrote {args.output} ({n} proteins)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
