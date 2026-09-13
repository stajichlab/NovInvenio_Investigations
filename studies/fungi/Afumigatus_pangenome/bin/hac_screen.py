#!/usr/bin/env python3
"""Targeted screen for the HAC (hrmA-Associated Cluster, PF11001/
IPR047092 paralog family) and the independent hacA (Afu3g04070) locus
across all 293 strains -- notes/superpowers/specs/2026-09-13-pangenome-
cluster-profile-design.md, component 7. Run only after Task 10's
benchmark scorecard has validated the clustering backend being used
here (per the spec's explicit ask to test recovery of known Starships
before trusting the method on this target).

Usage:
  hac_screen.py --matrix presence_matrix.rescued.tsv \
      --hac_family_id <tier1_rep_for_the_PF11001_family> \
      --haca_family_id <tier1_rep_for_hacA> \
      --strain_starship_map strain_to_starship.tsv \
      --output hac_screen_report.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402


def screen_family(
    matrix: PresenceMatrix, family_id: str, strain_to_starship: dict[str, str]
) -> list[dict]:
    rows = []
    for strain in matrix.strains:
        state = matrix.call(family_id, strain)
        rows.append({
            "strain": strain,
            "state": state,
            "copy_number": matrix.copy_number.get((family_id, strain), 0),
            "starship": strain_to_starship.get(strain, "unknown"),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--hac_family_id", required=True)
    ap.add_argument("--haca_family_id", required=True)
    ap.add_argument("--strain_starship_map", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    strain_to_starship: dict[str, str] = {}
    with open(args.strain_starship_map) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                strain_to_starship[parts[0]] = parts[1]

    hac_rows = screen_family(matrix, args.hac_family_id, strain_to_starship)
    haca_rows = screen_family(matrix, args.haca_family_id, strain_to_starship)

    with open(args.output, "w") as fh:
        fh.write("locus\tstrain\tstate\tcopy_number\tstarship\n")
        for row in hac_rows:
            fh.write(f"HAC\t{row['strain']}\t{row['state']}\t{row['copy_number']}\t{row['starship']}\n")
        for row in haca_rows:
            fh.write(f"hacA\t{row['strain']}\t{row['state']}\t{row['copy_number']}\t{row['starship']}\n")


if __name__ == "__main__":
    main()
