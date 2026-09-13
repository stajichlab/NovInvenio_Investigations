#!/usr/bin/env python3
"""Strain inventory: per-strain assembly-quality proxy (N50/contig count)
and mash-based dereplication, feeding every later step of the pangenome
cluster-profile analysis (see notes/superpowers/specs/
2026-09-13-pangenome-cluster-profile-design.md, component 1b).

Usage:
  dereplicate_strains.py --config config.csv --data_dir data_dir \
      --output strain_inventory.tsv
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def compute_assembly_stats(fasta_path: str | Path) -> dict:
    """Return {n_contigs, n50, total_length} for a DNA FASTA, without
    external tools -- a fast, dependency-free fragmentation proxy."""
    lengths = []
    current = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
        if current:
            lengths.append(current)
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break
    return {"n_contigs": len(lengths), "n50": n50, "total_length": total}


def parse_mash_dist(lines: list[str], threshold: float) -> list[set[str]]:
    """Parse `mash dist -t` output (a query x reference distance matrix,
    tab-separated, first row is '#query\\tref1\\tref2\\t...') into
    dereplication groups: strains whose pairwise mash distance is below
    `threshold` are grouped together (union-find over the threshold graph)."""
    if not lines:
        return []
    header = lines[0].lstrip("#").split("\t")
    names = header[1:]
    parent = {name: name for name in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for row in lines[1:]:
        parts = row.split("\t")
        query = parts[0]
        for ref, dist_str in zip(names, parts[1:]):
            if query == ref:
                continue
            if float(dist_str) < threshold:
                union(query, ref)

    groups: dict[str, set[str]] = {}
    for name in names:
        root = find(name)
        groups.setdefault(root, set()).add(name)
    return list(groups.values())


def choose_representatives(
    dedup_groups: list[set[str]], assembly_stats: dict[str, dict]
) -> dict[str, str]:
    """For every strain in every group, map it to the group's chosen
    representative (the member with the highest N50 -- the least
    fragmented assembly, per the review's fragmentation-risk concern)."""
    result: dict[str, str] = {}
    for group in dedup_groups:
        # Use deterministic tie-break on equal N50 (lexically last Short ID)
        rep = max(group, key=lambda s: (assembly_stats[s]["n50"], s))
        for member in group:
            result[member] = rep
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--mash_threshold", type=float, default=0.001)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    # Import config_parser only when main() is called (not at module import time)
    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    NOVINVENIO_LIB = Path(__file__).resolve().parents[4] / "NovInvenio" / "lib"
    if NOVINVENIO_LIB.exists():
        sys.path.insert(0, str(NOVINVENIO_LIB))
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(args.config)
    data_dir = Path(args.data_dir)
    dna_paths = {s.short: data_dir / "dna" / s.dna for s in samples if s.dna}

    stats = {short: compute_assembly_stats(p) for short, p in dna_paths.items()}

    sketch_prefix = "strain_sketches"
    subprocess.run(
        ["mash", "sketch", "-o", sketch_prefix] + [str(p) for p in dna_paths.values()],
        check=True,
    )
    dist_out = subprocess.run(
        ["mash", "dist", "-t", f"{sketch_prefix}.msh", f"{sketch_prefix}.msh"],
        check=True, capture_output=True, text=True,
    ).stdout
    groups = parse_mash_dist(dist_out.splitlines(), args.mash_threshold)

    # Translate groups from full paths (as reported by mash) to Short IDs
    path_to_short = {str(p): short for short, p in dna_paths.items()}
    groups = [{path_to_short[p] for p in group} for group in groups]

    reps = choose_representatives(groups, stats)
    dedup_group_id = {}
    for i, group in enumerate(groups):
        for member in group:
            dedup_group_id[member] = i

    with open(args.output, "w") as fh:
        fh.write("Short\tn_contigs\tn50\ttotal_length\tdedup_group\tis_representative\n")
        for short in sorted(dna_paths):
            s = stats[short]
            fh.write(
                f"{short}\t{s['n_contigs']}\t{s['n50']}\t{s['total_length']}\t"
                f"{dedup_group_id[short]}\t{int(reps[short] == short)}\n"
            )


if __name__ == "__main__":
    main()
