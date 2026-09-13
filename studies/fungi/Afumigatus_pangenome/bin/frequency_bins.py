#!/usr/bin/env python3
"""Core/soft-core/shell/cloud/singleton frequency binning -- notes/
superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 3.

Usage:
  frequency_bins.py --matrix presence_matrix.rescued.tsv --output frequency_table.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402


def assign_bin(
    freq: float,
    strain_count: int,
    core_cutoff: float = 0.95,
    softcore_cutoff: float = 0.90,
    shell_cutoff: float = 0.15,
) -> str:
    if strain_count <= 1:
        return "singleton"
    if freq >= core_cutoff:
        return "core"
    if freq >= softcore_cutoff:
        return "soft_core"
    if freq >= shell_cutoff:
        return "shell"
    return "cloud"


def compute_frequency_table(
    matrix: PresenceMatrix,
    core_cutoff: float = 0.95,
    softcore_cutoff: float = 0.90,
    shell_cutoff: float = 0.15,
) -> list[dict]:
    rows = []
    for fam in matrix.families:
        freq = matrix.frequency(fam)
        count = matrix.strain_count(fam)
        rows.append({
            "family": fam,
            "frequency": freq,
            "strain_count": count,
            "bin": assign_bin(freq, count, core_cutoff, softcore_cutoff, shell_cutoff),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--core_cutoff", type=float, default=0.95)
    ap.add_argument("--softcore_cutoff", type=float, default=0.90)
    ap.add_argument("--shell_cutoff", type=float, default=0.15)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    table = compute_frequency_table(
        matrix, args.core_cutoff, args.softcore_cutoff, args.shell_cutoff
    )
    with open(args.output, "w") as fh:
        fh.write("family\tfrequency\tstrain_count\tbin\n")
        for row in table:
            fh.write(f"{row['family']}\t{row['frequency']:.4f}\t{row['strain_count']}\t{row['bin']}\n")


if __name__ == "__main__":
    main()
