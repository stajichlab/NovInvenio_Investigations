#!/usr/bin/env python3
"""Adds `pfam_accession` and `pfam_url` columns to domain_enrichment.tsv
(from summarize_island_functions.py), so every domain in the enrichment
table links directly to its canonical Pfam/InterPro entry page --
requested for the report's hotlinks.

The accession comes from the SAME hmmscan --domtblout files already
produced (column 2, target accession -- hmmscan's domain-name-keyed
output doesn't carry the accession through to summarize_island_functions
.py's own tables, so this re-derives it from the raw domtblout rather
than re-running hmmscan). Pfam is now hosted under InterPro
(pfam.xfam.org is retired) -- the canonical entry URL is
https://www.ebi.ac.uk/interpro/entry/pfam/<ACCESSION>/, accession only
(no .NN version suffix, which domtblout includes but the URL doesn't
use).

Usage:
  add_pfam_urls.py --domain_enrichment domain_enrichment.tsv \\
      --domtblout island_family_reps_vs_pfam.domtblout \\
      --domtblout extra_background_vs_pfam.domtblout \\
      --output domain_enrichment.with_urls.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys

PFAM_URL_TEMPLATE = "https://www.ebi.ac.uk/interpro/entry/pfam/{accession}/"


def load_name_to_accession(domtblout_paths: list[str]) -> dict[str, str]:
    """{domain_name: bare_accession} (version suffix stripped, e.g.
    "PF00109.33" -> "PF00109") from one or more hmmscan --domtblout
    files' target-name/target-accession columns (1 and 2)."""
    mapping: dict[str, str] = {}
    for path in domtblout_paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                name, accession = parts[0], parts[1]
                mapping[name] = accession.split(".", 1)[0]
    return mapping


def add_urls(rows: list[dict], name_to_accession: dict[str, str]) -> list[dict]:
    annotated = []
    for row in rows:
        new_row = dict(row)
        accession = name_to_accession.get(row["domain"], "")
        new_row["pfam_accession"] = accession
        new_row["pfam_url"] = PFAM_URL_TEMPLATE.format(accession=accession) if accession else ""
        annotated.append(new_row)
    return annotated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domain_enrichment", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    name_to_accession = load_name_to_accession(args.domtblout)
    with open(args.domain_enrichment, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    annotated = add_urls(rows, name_to_accession)
    n_missing = sum(1 for r in annotated if not r["pfam_accession"])
    if n_missing:
        print(f"add_pfam_urls: WARNING {n_missing}/{len(annotated)} domains had no "
              "accession found in the given --domtblout files", file=sys.stderr)

    fieldnames = list(annotated[0].keys()) if annotated else []
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)
    print(f"add_pfam_urls: wrote {len(annotated)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
