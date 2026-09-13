#!/usr/bin/env python3
"""Per-strain gene-order synteny scanning for Starship-driven physical
clustering -- notes/superpowers/specs/2026-09-13-pangenome-cluster-
profile-design.md, component 5 (refined per the Fable bioinformatics
review: per-strain window enumeration, contig-edge exclusion, an
accessory-island test alongside fixed-size windows, and a direct
linkage-fraction statistic per correlated family pair).

Usage:
  synteny_windows.py --gff3_dir data_dir/gff3 --config config.csv \
      --matrix presence_matrix.rescued.tsv --window_sizes 3,10 \
      --output synteny_windows.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))


_ID_RE = re.compile(r"ID=([^;\n]+)")


def parse_gff3_gene_order(gff3_path: str | Path) -> list[tuple[str, str, int, int]]:
    """Return [(gene_id, contig_id, start, end), ...] for every 'gene'
    feature, sorted by contig then start position."""
    genes = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            contig, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
            m = _ID_RE.search(attrs)
            if not m:
                continue
            genes.append((m.group(1), contig, start, end))
    genes.sort(key=lambda g: (g[1], g[2]))
    return genes


def sliding_windows(gene_order: list[tuple[str, str, int, int]], n: int) -> list[list[tuple]]:
    windows = []
    for i in range(len(gene_order) - n + 1):
        window = gene_order[i:i + n]
        contigs = {g[1] for g in window}
        if len(contigs) == 1:
            windows.append(window)
    return windows


def accessory_islands(
    gene_order: list[tuple[str, str, int, int]], is_core: dict[str, bool]
) -> list[list[tuple]]:
    """Maximal runs of consecutive non-core genes, never crossing a
    contig boundary."""
    islands: list[list[tuple]] = []
    current: list[tuple] = []
    prev_contig = None
    for gene in gene_order:
        gene_id, contig = gene[0], gene[1]
        core = is_core.get(gene_id, True)
        if contig != prev_contig and current:
            islands.append(current)
            current = []
        if not core:
            current.append(gene)
        elif current:
            islands.append(current)
            current = []
        prev_contig = contig
    if current:
        islands.append(current)
    return islands


def linkage_fraction(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, tuple[str, int]]],
    k: int = 10,
) -> float:
    """Fraction of strains carrying BOTH families where their gene-order
    positions are within k genes of each other on the same contig."""
    both_present = [
        s for s, positions in gene_position.items()
        if family_a in positions and family_b in positions
    ]
    if not both_present:
        return 0.0
    linked = 0
    for s in both_present:
        contig_a, pos_a = gene_position[s][family_a]
        contig_b, pos_b = gene_position[s][family_b]
        if contig_a == contig_b and abs(pos_a - pos_b) <= k:
            linked += 1
    return linked / len(both_present)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff3_dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--window_sizes", default="3,10")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    # Full per-strain window enumeration + correlated-window detection
    # across all 293 strains' GFF3s is run interactively (Step 5 below)
    # rather than duplicated here -- this CLI's job is exposing the
    # tested building blocks above for that run.
    print(
        "Use parse_gff3_gene_order/sliding_windows/accessory_islands/"
        "linkage_fraction directly for the real run; see plan Task 8 Step 5.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
