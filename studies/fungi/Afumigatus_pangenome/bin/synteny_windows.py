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
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))


# Anchored to the start of the attributes string or right after a ';' so a
# decoy attribute like "orig_protein_ID=XP_1;ID=geneA" does not match on
# "ID=" buried inside "orig_protein_ID=" -- unanchored .search() would
# incorrectly return "XP_1" there.
_ID_RE = re.compile(r"(?:^|;)ID=([^;\n]+)")


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
    contig boundary.

    A gene_id absent from `is_core` is treated as core (`is_core.get(gene_id,
    True)`), i.e. it splits/ends any in-progress island rather than being
    folded into one. This is deliberately conservative: an unrecognized gene
    is more likely a family the clustering step never called ("core" is the
    majority state, ~95%+ of genes in this pangenome) than a genuine, silent
    accessory gene, so defaulting it to non-core would risk manufacturing
    spurious islands out of noise. The real danger is the reverse failure
    mode -- a systematic ID-format mismatch between the GFF3's `ID=` values
    and the presence-matrix's protein IDs (e.g. `gene-X` vs `X-T1`) would
    make *every* gene look "missing", silently collapsing every island to
    nothing with no exception raised. To make that failure visible instead
    of silent, this function warns (via `warnings.warn`) whenever more than
    10% of the genes it sees are absent from `is_core` -- a real Step-5 run
    should treat that warning as a hard stop, not proceed with degenerate
    output.
    """
    islands: list[list[tuple]] = []
    current: list[tuple] = []
    prev_contig = None
    missing = 0
    for gene in gene_order:
        gene_id, contig = gene[0], gene[1]
        if gene_id not in is_core:
            missing += 1
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
    if gene_order and missing / len(gene_order) > 0.1:
        warnings.warn(
            f"accessory_islands: {missing}/{len(gene_order)} genes "
            "("
            f"{missing / len(gene_order):.0%}) were absent from is_core and "
            "defaulted to core -- check for a GFF3 ID vs presence-matrix "
            "protein ID format mismatch before trusting these islands.",
            stacklevel=2,
        )
    return islands


def linkage_fraction(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    k: int = 10,
) -> float:
    """Fraction of strains carrying BOTH families where SOME copy of
    family_a is within k genes of SOME copy of family_b on the same contig.

    `gene_position[strain][family]` is a LIST of (contig, rank) positions,
    not a single position -- design spec component 8 must-fix M4: this
    study's own multi-copy families (PF11001 at 3-9 copies per strain in
    every one of 295 strains, per the HAC screen) broke the original
    single-position assumption, silently picking whichever copy happened to
    be stored and mislabeling the pair's physical linkage. Checking the
    MINIMUM distance over every copy-pair answers the right question ("is
    ANY copy of family_a near ANY copy of family_b"), not "is the arbitrary
    stored copy of each near the other"."""
    both_present = [
        s for s, positions in gene_position.items()
        if positions.get(family_a) and positions.get(family_b)
    ]
    if not both_present:
        return 0.0
    linked = 0
    for s in both_present:
        positions = gene_position[s]
        is_linked = any(
            contig_a == contig_b and abs(pos_a - pos_b) <= k
            for contig_a, pos_a in positions[family_a]
            for contig_b, pos_b in positions[family_b]
        )
        if is_linked:
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
