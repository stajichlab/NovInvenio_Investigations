#!/usr/bin/env python3
"""Adds a `swissprot_names` column to significant_islands.with_enrichment.tsv
(or any of the island tables sharing its `member_families` column) --
each island's real, human-readable protein names/descriptions from
family_swissprot_annotation.tsv (annotate_families_with_swissprot.py),
instead of only cryptic family IDs and Pfam domain codes.

Usage:
  add_swissprot_to_islands.py --islands significant_islands.with_enrichment.tsv \\
      --swissprot_annotation family_swissprot_annotation.tsv \\
      --output significant_islands.with_enrichment.with_swissprot.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys


def load_family_names(path: str) -> dict[str, str]:
    """{family_id: "NAME (description)"} -- description omitted when
    blank (the self-annotated case, where the UniProt name field alone
    is already informative)."""
    names: dict[str, str] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            label = row["name"]
            if row["description"]:
                label += f" ({row['description']})"
            names[row["family"]] = label
    return names


def add_swissprot_names(island_rows: list[dict], family_names: dict[str, str]) -> list[dict]:
    annotated = []
    for row in island_rows:
        members = row["member_families"].split(",")
        hits = [family_names[m] for m in members if m in family_names]
        new_row = dict(row)
        new_row["n_swissprot_hits"] = str(len(hits))
        new_row["swissprot_names"] = "; ".join(hits) if hits else "-"
        annotated.append(new_row)
    return annotated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--islands", required=True)
    ap.add_argument("--swissprot_annotation", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    family_names = load_family_names(args.swissprot_annotation)
    with open(args.islands, newline="") as fh:
        island_rows = list(csv.DictReader(fh, delimiter="\t"))

    annotated = add_swissprot_names(island_rows, family_names)
    n_with_hit = sum(1 for r in annotated if int(r["n_swissprot_hits"]) > 0)
    print(
        f"add_swissprot_to_islands: {n_with_hit}/{len(annotated)} islands have "
        ">=1 member with a SwissProt-based name", file=sys.stderr,
    )

    fieldnames = list(annotated[0].keys()) if annotated else []
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)


if __name__ == "__main__":
    main()
