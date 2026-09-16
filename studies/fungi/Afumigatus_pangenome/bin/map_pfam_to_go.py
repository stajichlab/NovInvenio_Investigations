#!/usr/bin/env python3
"""Adds GO-term annotations to domain_enrichment.with_urls.tsv (or any
table carrying a `pfam_accession` column) via the standard Pfam2GO
mapping -- a near-free way to add GO-level functional information on top
of the Pfam domains already found, since Pfam is itself an InterPro
member database and this mapping (InterPro's own, redistributed via the
Gene Ontology Consortium) already ties every Pfam accession to its
associated GO terms. Deliberately NOT running a fresh InterProScan (a
much heavier, multi-hour job for this many sequences) -- this only
annotates domains ALREADY found by the existing Pfam-A hmmscan, it
cannot discover a new domain hmmscan itself missed.

Mapping file: standard Pfam2GO release, e.g.
  curl -sL http://current.geneontology.org/ontology/external2go/pfam2go \\
      -o pfam2go.txt
(ephemeral input, re-fetchable -- not itself a tracked data file; see
this study's CLAUDE.md class-1 data-provenance convention.)

Usage:
  map_pfam_to_go.py --domain_enrichment domain_enrichment.with_urls.tsv \\
      --pfam2go pfam2go.txt \\
      --output domain_enrichment.with_go.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

_PFAM2GO_LINE_RE = re.compile(r"^Pfam:(PF\d+)\s+\S+\s+>\s+GO:(.+?)\s*;\s*(GO:\d+)\s*$")


def load_pfam_to_go(path: str) -> dict[str, list[tuple[str, str]]]:
    """{pfam_accession: [(go_id, go_term_name), ...]} from a standard
    pfam2go release file (lines like
    "Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930"
    -- comment lines starting with "!" are skipped)."""
    mapping: dict[str, list[tuple[str, str]]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("!") or not line.strip():
                continue
            m = _PFAM2GO_LINE_RE.match(line.strip())
            if not m:
                continue
            accession, term_name, go_id = m.groups()
            mapping.setdefault(accession, []).append((go_id, term_name))
    return mapping


def add_go_terms(rows: list[dict], pfam_to_go: dict[str, list[tuple[str, str]]]) -> list[dict]:
    annotated = []
    for row in rows:
        terms = pfam_to_go.get(row.get("pfam_accession", ""), [])
        new_row = dict(row)
        new_row["n_go_terms"] = str(len(terms))
        new_row["go_terms"] = "; ".join(f"{go_id} ({name})" for go_id, name in terms) if terms else ""
        annotated.append(new_row)
    return annotated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domain_enrichment", required=True)
    ap.add_argument("--pfam2go", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    pfam_to_go = load_pfam_to_go(args.pfam2go)
    print(f"map_pfam_to_go: {len(pfam_to_go)} Pfam accessions have >=1 GO term in this mapping",
          file=sys.stderr)

    with open(args.domain_enrichment, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    annotated = add_go_terms(rows, pfam_to_go)
    n_with_go = sum(1 for r in annotated if int(r["n_go_terms"]) > 0)
    print(f"map_pfam_to_go: {n_with_go}/{len(annotated)} domains in this table have >=1 GO term",
          file=sys.stderr)

    fieldnames = list(annotated[0].keys()) if annotated else []
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)


if __name__ == "__main__":
    main()
