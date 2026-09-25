#!/usr/bin/env python3
"""Collapse the `trans` (physically unlinked, statistically significant)
co-occurrence pairs from a pangenome study's pair_classification.tsv into
interpretable modules via Leiden community detection.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/detect_trans_modules.py) for the coccidioides_pangenome study. Only
change from the original: plain-file I/O (no compressed_io dependency --
this study's pair_classification.tsv outputs aren't gzipped), everything
else (graph construction, Leiden call, module renumbering) is unchanged.
See that script's docstring for the Leiden-vs-Louvain rationale (Traag,
Waltman & van Eck 2019, Sci Rep 9:5233) and the resolution-tuning caveat --
no published pangenome-specific precedent exists for this parameter at
this scale; `--resolution` is exposed for empirical tuning, not trusted
blindly at its default.

Usage:
  detect_trans_modules.py --pair_classification pair_classification.tsv \\
      --resolution 1.0 --seed 0 \\
      --output_families family_modules.tsv \\
      --output_modules module_summary.tsv
"""
from __future__ import annotations

import argparse
import math
import sys

import igraph as ig
import leidenalg

MAX_WEIGHT = 300.0  # clip -log10(fdr_q) for fdr_q underflowing to 0.0


def load_trans_edges(path: str) -> list[tuple[str, str, float]]:
    """Returns [(family_a, family_b, weight), ...] for classification=="trans"
    rows, weight = -log10(fdr_q) (clipped) -- a stronger statistical signal
    (smaller q) gets a heavier edge."""
    edges = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[idx["classification"]] != "trans":
                continue
            fam_a, fam_b = parts[idx["family_a"]], parts[idx["family_b"]]
            q = float(parts[idx["fdr_q"]])
            weight = MAX_WEIGHT if q <= 0 else min(MAX_WEIGHT, -math.log10(q))
            edges.append((fam_a, fam_b, weight))
    return edges


def build_graph(edges: list[tuple[str, str, float]]) -> ig.Graph:
    nodes = sorted({n for a, b, _ in edges for n in (a, b)})
    node_index = {n: i for i, n in enumerate(nodes)}
    g = ig.Graph()
    g.add_vertices(len(nodes))
    g.vs["name"] = nodes
    g.add_edges([(node_index[a], node_index[b]) for a, b, _ in edges])
    g.es["weight"] = [w for _, _, w in edges]
    return g


def detect_modules(
    g: ig.Graph, resolution: float = 1.0, seed: int = 0
) -> dict[str, int]:
    """Returns {family: module_id}, module IDs assigned in descending
    order of module size (module 0 is the largest)."""
    partition = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution, seed=seed,
    )
    sizes = [(i, len(members)) for i, members in enumerate(partition)]
    sizes.sort(key=lambda x: -x[1])
    old_to_new = {old: new for new, (old, _) in enumerate(sizes)}
    family_module: dict[str, int] = {}
    for old_id, members in enumerate(partition):
        new_id = old_to_new[old_id]
        for v in members:
            family_module[g.vs[v]["name"]] = new_id
    return family_module


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output_families", required=True)
    ap.add_argument("--output_modules", required=True)
    args = ap.parse_args()

    print("Loading trans edges...", file=sys.stderr)
    edges = load_trans_edges(args.pair_classification)
    print(f"{len(edges)} trans edges loaded", file=sys.stderr)
    if not edges:
        print("No trans edges found -- nothing to cluster.", file=sys.stderr)
        sys.exit(1)

    g = build_graph(edges)
    print(f"Graph: {g.vcount()} nodes, {g.ecount()} edges", file=sys.stderr)

    print(f"Running Leiden (resolution={args.resolution}, seed={args.seed})...", file=sys.stderr)
    family_module = detect_modules(g, resolution=args.resolution, seed=args.seed)

    module_sizes: dict[int, int] = {}
    for module_id in family_module.values():
        module_sizes[module_id] = module_sizes.get(module_id, 0) + 1

    with open(args.output_families, "w") as out:
        out.write("family\tmodule_id\tmodule_size\n")
        for family, module_id in sorted(family_module.items(), key=lambda kv: (kv[1], kv[0])):
            out.write(f"{family}\t{module_id}\t{module_sizes[module_id]}\n")

    with open(args.output_modules, "w") as out:
        out.write("module_id\tsize\n")
        for module_id, size in sorted(module_sizes.items(), key=lambda kv: -kv[1]):
            out.write(f"{module_id}\t{size}\n")

    n_modules = len(module_sizes)
    largest = max(module_sizes.values())
    singletons = sum(1 for s in module_sizes.values() if s == 1)
    print(
        f"detect_trans_modules: {n_modules} modules, largest={largest}, "
        f"singleton modules={singletons}", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
