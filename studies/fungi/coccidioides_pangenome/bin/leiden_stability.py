#!/usr/bin/env python3
"""Empirically test whether a Leiden trans-co-occurrence module resolution
has a defensible value: at each candidate resolution, rerun Leiden with
many different random seeds and measure how much the resulting partition
agrees across seeds via Adjusted Mutual Information (AMI). A resolution
where independent seeds converge on nearly the same partition (AMI near 1)
reflects real, reproducible community structure; low AMI means the
partition is close to arbitrary tie-breaking among near-degenerate
solutions, not a stable answer.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/leiden_stability.py) for the coccidioides_pangenome study, unchanged
except importing detect_trans_modules from this study's own bin/.

Usage:
  leiden_stability.py --trans_edges trans_edges.tsv \\
      --resolutions 0.5,1.0,2.0,5.0,10.0 --n_seeds 10 \\
      --output leiden_stability.tsv

`trans_edges.tsv` is a pre-extracted (family_a, family_b, fdr_q) cache,
e.g.:
  awk -F'\\t' 'NR==1 || $3=="trans"{print $1"\\t"$2"\\t"$7}' pair_classification.tsv \\
      > trans_edges.tsv
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from itertools import combinations
from pathlib import Path

from sklearn.metrics import adjusted_mutual_info_score

sys.path.insert(0, str(Path(__file__).parent))
from detect_trans_modules import build_graph, detect_modules  # noqa: E402


def load_cached_trans_edges(path: str) -> list[tuple[str, str, float]]:
    edges = []
    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("family_a"):
            fh.seek(0)
        for line in fh:
            fam_a, fam_b, fdr_q = line.rstrip("\n").split("\t")
            q = float(fdr_q)
            weight = 300.0 if q <= 0 else min(300.0, -math.log10(q))
            edges.append((fam_a, fam_b, weight))
    return edges


def align_labels(a: dict, b: dict, families: list[str]) -> tuple[list, list]:
    """Two partitions' label *values* are arbitrary (Leiden doesn't number
    modules consistently run-to-run) -- AMI/ARI only care about which
    families share a label, so callers must present both partitions'
    labels in the same family order, not compare raw label dicts."""
    return [a[f] for f in families], [b[f] for f in families]


def pairwise_ami(partitions: list[dict], families: list[str]) -> tuple[float, float]:
    scores = []
    for p1, p2 in combinations(partitions, 2):
        la, lb = align_labels(p1, p2, families)
        scores.append(adjusted_mutual_info_score(la, lb))
    mean_ami = statistics.mean(scores)
    std_ami = statistics.pstdev(scores) if len(scores) > 1 else 0.0
    return mean_ami, std_ami


def singleton_fraction(partition: dict) -> float:
    sizes: dict = {}
    for module_id in partition.values():
        sizes[module_id] = sizes.get(module_id, 0) + 1
    n_singleton = sum(1 for s in sizes.values() if s == 1)
    return n_singleton / len(sizes)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trans_edges", required=True, help="pre-extracted family_a/family_b/fdr_q cache")
    ap.add_argument("--resolutions", default="0.5,1.0,2.0,5.0,10.0")
    ap.add_argument("--n_seeds", type=int, default=10)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    print("Loading cached trans edges...", file=sys.stderr)
    edges = load_cached_trans_edges(args.trans_edges)
    print(f"{len(edges)} trans edges loaded", file=sys.stderr)

    g = build_graph(edges)
    families = g.vs["name"]
    print(f"Graph: {g.vcount()} nodes, {g.ecount()} edges", file=sys.stderr)

    resolutions = [float(r) for r in args.resolutions.split(",")]
    rows = []
    for resolution in resolutions:
        partitions = []
        module_counts = []
        singleton_fracs = []
        for seed in range(args.n_seeds):
            partition = detect_modules(g, resolution=resolution, seed=seed)
            partitions.append(partition)
            module_counts.append(len(set(partition.values())))
            singleton_fracs.append(singleton_fraction(partition))
        mean_ami, std_ami = pairwise_ami(partitions, families)
        row = {
            "resolution": resolution,
            "n_seeds": args.n_seeds,
            "mean_pairwise_ami": f"{mean_ami:.4f}",
            "std_pairwise_ami": f"{std_ami:.4f}",
            "median_module_count": statistics.median(module_counts),
            "median_singleton_fraction": f"{statistics.median(singleton_fracs):.4f}",
        }
        rows.append(row)
        print(f"  resolution={resolution}: mean_AMI={mean_ami:.4f} (std={std_ami:.4f}), "
              f"median_modules={row['median_module_count']}, "
              f"median_singleton_frac={row['median_singleton_fraction']}", file=sys.stderr)

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
