#!/usr/bin/env python3
"""Join gene_positions.tsv (protein_id -> genomic position, per strain)
with a tier-1 cluster TSV (protein_id -> family_id) into per-strain,
per-family gene-order RANK positions -- the input `synteny_windows.py`'s
`linkage_fraction`/`accessory_islands` and `bin/pair_classification.py`
actually consume.

Rank, not raw bp coordinate: `linkage_fraction`'s `k` parameter is a gene-
COUNT window ("within k genes of each other"), matching its original
design and tests -- a protein's rank is its index (0-based) among all
proteins on the SAME contig, sorted by start position. A family with
multiple copies in one strain (this study's own PF11001 finding: 3-9
copies per strain) gets a LIST of (contig, rank) positions per strain, not
one -- see synteny_windows.py's must-fix-M4 fix to `linkage_fraction`.

Usage:
  build_family_positions.py --gene_positions gene_positions.tsv \\
      --cluster_tsv tier1_cluster.tsv --output family_positions.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv  # noqa: E402


def build_family_positions(
    gene_position_rows: list[tuple[str, str, str, int, int]],
    member_to_rep: dict[str, str],
    id_sep: str = "|",
) -> dict[str, dict[str, list[tuple[str, int]]]]:
    """`gene_position_rows`: [(Short, protein_id, contig, start, end), ...].
    Returns {Short: {family_id: [(contig, rank), ...]}}.

    A protein_id with no cluster membership (e.g. a rare gene filtered out
    before clustering) is silently skipped -- not every protein needs to be
    in a family for this to work, only ones that are."""
    # First pass: per-strain ordered (protein_id, contig, start) list, sorted
    # by contig then start, so a running index gives contig-relative rank
    # (a jump in the running index at a contig boundary is harmless --
    # linkage_fraction only ever compares positions when contig_a == contig_b).
    per_strain_ordered: dict[str, list[tuple[str, str, int]]] = {}
    for short, protein_id, contig, start, _end in gene_position_rows:
        per_strain_ordered.setdefault(short, []).append((protein_id, contig, start))

    result: dict[str, dict[str, list[tuple[str, int]]]] = {}
    for short, ordered in per_strain_ordered.items():
        ordered.sort(key=lambda row: (row[1], row[2]))
        family_positions: dict[str, list[tuple[str, int]]] = {}
        for rank, (protein_id, contig, _start) in enumerate(ordered):
            member_id = f"{short}{id_sep}{protein_id}"
            family = member_to_rep.get(member_id)
            if family is None:
                continue
            family_positions.setdefault(family, []).append((contig, rank))
        result[short] = family_positions
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    member_to_rep = read_cluster_tsv(args.cluster_tsv)

    rows = []
    with open(args.gene_positions) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            short, protein_id, contig, start, end = parts
            rows.append((short, protein_id, contig, int(start), int(end)))

    positions = build_family_positions(rows, member_to_rep, id_sep=args.id_sep)

    n_rows = 0
    with open(args.output, "w") as out:
        out.write("Short\tfamily\tcontig\trank\n")
        for short, family_positions in positions.items():
            for family, entries in family_positions.items():
                for contig, rank in entries:
                    out.write(f"{short}\t{family}\t{contig}\t{rank}\n")
                    n_rows += 1

    print(f"build_family_positions: {n_rows} (strain, family, copy) positions written", file=sys.stderr)


if __name__ == "__main__":
    main()
