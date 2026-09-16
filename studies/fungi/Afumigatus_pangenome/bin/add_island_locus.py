#!/usr/bin/env python3
"""Adds a stable, human-readable genomic locus identifier to each
significant island (e.g. `Asfu_Z5:JAIBVV010000123.1:104512-118930`), using
the island's `example_strain` occurrence as the reference coordinate.

Without this, an island is only identifiable by its (arbitrary-order)
comma-joined member_families list -- unwieldy to cite or look up in a
genome browser. Resolves real bp coordinates via the SAME two files the
rest of the pipeline already builds/uses, rather than re-parsing GFF3:

- `tier1_cluster.tsv` (rep\\tmember per line): inverted to find which of the
  example_strain's own proteins belongs to each member family, since a
  family ID is itself the rep protein's ID (usually from a DIFFERENT
  strain) and cannot be looked up directly against gene_positions.tsv.
- `gene_positions.tsv` (Short/protein_id/contig/start/end, from
  build_gene_positions.py's GFF3 parse): the actual bp coordinates for
  that resolved protein.

A family member with no resolvable protein in the example_strain (most
often a rescue-pass GENOME_ONLY call, which has a family-level tblastn
position but no annotated protein_id -- see extract_rescue_positions.py)
is simply excluded from the span, and counted in n_members_with_coordinates
vs. the island's full member count, rather than causing a crash or a
silently wrong span.

Usage:
  add_island_locus.py --significant_islands significant_islands.with_enrichment.tsv \\
      --cluster_tsv results/full_293run/tier1_cluster.tsv \\
      --gene_positions results/full_293run/gene_positions.tsv \\
      --output significant_islands.with_locus.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402


def load_needed_families(significant_islands_path: str) -> set[str]:
    needed = set()
    with open(significant_islands_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            needed.update(row["member_families"].split(","))
    return needed


def load_family_members_by_strain(
    cluster_tsv_path: str, needed_families: set[str]
) -> dict[str, dict[str, list[str]]]:
    """{family: {strain: [protein_id, ...]}}, restricted to families that
    actually appear in some significant island -- the full cluster TSV has
    millions of rows, most irrelevant here."""
    result: dict[str, dict[str, list[str]]] = {}
    with open_maybe_compressed(cluster_tsv_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep, member = line.split("\t", 1)
            if rep not in needed_families:
                continue
            if "|" not in member:
                continue
            strain, protein_id = member.split("|", 1)
            result.setdefault(rep, {}).setdefault(strain, []).append(protein_id)
    return result


def load_gene_positions_for_strains(
    gene_positions_path: str, needed_strains: set[str]
) -> dict[tuple[str, str], tuple[str, int, int]]:
    """{(strain, protein_id): (contig, start, end)}, restricted to strains
    that are some island's example_strain -- the full file covers all 295
    strains' full protein sets."""
    result: dict[tuple[str, str], tuple[str, int, int]] = {}
    with open(gene_positions_path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            short, protein_id, contig, start, end = parts
            if short not in needed_strains:
                continue
            result[(short, protein_id)] = (contig, int(start), int(end))
    return result


def compute_locus(
    example_strain: str,
    member_families: list[str],
    family_members_by_strain: dict[str, dict[str, list[str]]],
    gene_positions: dict[tuple[str, str], tuple[str, int, int]],
) -> dict:
    """Resolves each member family to its example_strain occurrence's bp
    position (if any), then reports the min-start/max-end span. Reports
    n_contigs > 1 rather than silently spanning across contigs -- a real
    accessory island is a single-contig consecutive run by construction
    (see synteny_windows.accessory_islands()), so >1 here signals a data
    problem to investigate, not a valid wide locus."""
    positions = []
    for family in member_families:
        for protein_id in family_members_by_strain.get(family, {}).get(example_strain, []):
            pos = gene_positions.get((example_strain, protein_id))
            if pos is not None:
                positions.append(pos)

    n_total = len(member_families)
    n_resolved = len(positions)
    if not positions:
        return {
            "locus_id": "", "contig": "", "start": None, "end": None,
            "n_resolved": n_resolved, "n_total": n_total, "n_contigs": 0,
        }

    contigs = {p[0] for p in positions}
    contig = positions[0][0]
    start = min(p[1] for p in positions)
    end = max(p[2] for p in positions)
    return {
        "locus_id": f"{example_strain}:{contig}:{start}-{end}",
        "contig": contig, "start": start, "end": end,
        "n_resolved": n_resolved, "n_total": n_total, "n_contigs": len(contigs),
    }


def annotate_islands(
    rows: list[dict],
    family_members_by_strain: dict[str, dict[str, list[str]]],
    gene_positions: dict[tuple[str, str], tuple[str, int, int]],
) -> list[dict]:
    annotated = []
    for row in rows:
        member_families = row["member_families"].split(",")
        locus = compute_locus(row["example_strain"], member_families, family_members_by_strain, gene_positions)
        new_row = dict(row)
        new_row["locus_id"] = locus["locus_id"]
        new_row["locus_contig"] = locus["contig"]
        new_row["locus_start"] = "" if locus["start"] is None else str(locus["start"])
        new_row["locus_end"] = "" if locus["end"] is None else str(locus["end"])
        new_row["n_members_with_coordinates"] = str(locus["n_resolved"])
        new_row["n_contigs_in_locus"] = str(locus["n_contigs"])
        annotated.append(new_row)
    return annotated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_islands", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    needed_families = load_needed_families(args.significant_islands)
    print(f"add_island_locus: {len(needed_families)} distinct member families across all islands",
          file=sys.stderr)

    with open(args.significant_islands, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    needed_strains = {row["example_strain"] for row in rows}

    family_members_by_strain = load_family_members_by_strain(args.cluster_tsv, needed_families)
    gene_positions = load_gene_positions_for_strains(args.gene_positions, needed_strains)
    print(f"add_island_locus: resolved cluster membership for {len(family_members_by_strain)} families, "
          f"gene positions for {len(needed_strains)} example strains", file=sys.stderr)

    annotated = annotate_islands(rows, family_members_by_strain, gene_positions)
    n_full = sum(
        1 for r in annotated
        if r["n_members_with_coordinates"] and r["locus_id"]
        and int(r["n_members_with_coordinates"]) == len(r["member_families"].split(","))
    )
    n_multi_contig = sum(1 for r in annotated if r["n_contigs_in_locus"] not in ("0", "1"))
    print(f"add_island_locus: {n_full}/{len(annotated)} islands fully resolved (all members have "
          f"coordinates); {n_multi_contig} flagged with members spanning >1 contig", file=sys.stderr)

    fieldnames = list(annotated[0].keys()) if annotated else []
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)


if __name__ == "__main__":
    main()
