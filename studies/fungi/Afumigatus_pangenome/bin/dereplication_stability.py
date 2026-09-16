#!/usr/bin/env python3
"""Empirically test whether the Mash dereplication threshold has a
defensible value, per PANGENOME_CLUSTER_PROFILE_NOTES.md's 2026-09-15
"extreme sensitivity, no stable plateau" finding.

Two independent checks, both against real data already in this study:

1. External validation: score each threshold's dereplication grouping
   against the reference paper's own published population clusters
   (`results/clade_assignment/s21_matches.tsv`, S21_groupID -- 254/295
   strains, 86% coverage) via Adjusted Rand Index. A threshold that merges
   two distinct published clusters into one dereplication group, or that
   fragments one published cluster into many singletons, scores worse.

2. Linkage-method sensitivity: `bin/dereplicate_strains.py` uses
   single-linkage clustering (union-find over the threshold graph) via
   `parse_mash_dist()`, which is known to "chain" -- a run of pairwise-close
   strains can transitively merge into one group even though its endpoints
   are far apart. This script reruns the same threshold sweep with
   complete-linkage (`scipy.cluster.hierarchy`, which requires every pair
   within a cluster, not just a chain, to be below threshold) to check
   whether the instability documented on 2026-09-15 is intrinsic to the
   data or an artifact of single-linkage chaining specifically.

Usage:
  dereplication_stability.py --mash_dist all_mash_dist.tsv \\
      --ground_truth results/clade_assignment/s21_matches.tsv \\
      --config config.csv --data_dir data_dir \\
      --thresholds 0.0002,0.0004,...,0.003 \\
      --output dereplication_stability.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score


def load_pairwise_mash(path: str) -> list[tuple[str, str, float]]:
    """Parses `mash dist sketches.msh sketches.msh` all-vs-all output
    (query \\t ref \\t distance \\t pvalue \\t shared-hashes), skipping any
    stray non-data lines (e.g. pixi's cache-redirect WARN banner) and
    self-pairs."""
    pairs = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5:
                continue
            query, ref, dist_str = parts[0], parts[1], parts[2]
            if query == ref:
                continue
            try:
                dist = float(dist_str)
            except ValueError:
                continue
            pairs.append((query, ref, dist))
    return pairs


def load_ground_truth_clusters(path: str) -> dict[str, str]:
    """{Short: S21_groupID} from results/clade_assignment/s21_matches.tsv."""
    truth = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            truth[row["Short"]] = row["S21_groupID"]
    return truth


def load_path_to_short(config_path: str, data_dir: str) -> dict[str, str]:
    """{data_dir/dna/<DNA filename>: Short} from config.csv, matching how
    dereplicate_strains.py/assign_clades.py build their own mash sketch
    input paths."""
    mapping = {}
    with open(config_path) as fh:
        for row in csv.DictReader(fh):
            if row.get("DNA"):
                mapping[f"{data_dir}/dna/{row['DNA']}"] = row["Short"]
    return mapping


def _union_find(nodes: set[str]):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    return find, union


def single_linkage_groups(
    pairs: list[tuple[str, str, float]], threshold: float, nodes: set[str]
) -> list[set[str]]:
    """Same algorithm as dereplicate_strains.py's parse_mash_dist(): merge
    any two nodes whose distance is below threshold, transitively."""
    find, union = _union_find(nodes)
    for a, b, dist in pairs:
        if dist < threshold:
            union(a, b)
    groups: dict[str, set[str]] = {}
    for n in nodes:
        groups.setdefault(find(n), set()).add(n)
    return list(groups.values())


def hierarchical_groups(
    pairs: list[tuple[str, str, float]],
    threshold: float,
    nodes: set[str],
    method: str,
) -> list[set[str]]:
    """Complete- or average-linkage clustering (scipy) on the same distance
    data, cut at the same threshold -- unlike single-linkage, these require
    every (complete) or the average (average) pairwise distance within a
    cluster to be below threshold, not just a chain of neighbors."""
    ordered = sorted(nodes)
    idx = {n: i for i, n in enumerate(ordered)}
    n = len(ordered)
    dist_matrix = np.zeros((n, n))
    max_dist = max((d for _, _, d in pairs), default=1.0) * 2 + 1.0
    dist_matrix[:] = max_dist
    np.fill_diagonal(dist_matrix, 0.0)
    for a, b, dist in pairs:
        if a in idx and b in idx:
            i, j = idx[a], idx[b]
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    condensed = squareform(dist_matrix, checks=False)
    if n < 2:
        return [set(ordered)]
    Z = linkage(condensed, method=method)
    labels = fcluster(Z, t=threshold, criterion="distance")
    groups: dict[int, set[str]] = {}
    for name, label in zip(ordered, labels):
        groups.setdefault(label, set()).add(name)
    return list(groups.values())


def groups_to_labels(groups: list[set[str]]) -> dict[str, int]:
    labels = {}
    for i, group in enumerate(groups):
        for member in group:
            labels[member] = i
    return labels


def ari_against_ground_truth(
    pred_labels: dict[str, int], truth: dict[str, str]
) -> float | None:
    """Adjusted Rand Index between predicted group labels and the ground
    truth, restricted to strains present in both. Returns None (rather
    than raising or silently returning 0) when there is no overlap, since
    that's a data problem, not a score of 0 agreement."""
    common = sorted(set(pred_labels) & set(truth))
    if not common:
        return None
    return float(
        adjusted_rand_score([pred_labels[c] for c in common], [truth[c] for c in common])
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mash_dist", required=True, help="mash dist all-vs-all pairwise output")
    ap.add_argument("--ground_truth", required=True, help="results/clade_assignment/s21_matches.tsv")
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--thresholds", default="0.0002,0.0004,0.0006,0.0008,0.001,0.0012,0.0015,0.002,0.003,0.005")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    path_to_short = load_path_to_short(args.config, args.data_dir)
    raw_pairs = load_pairwise_mash(args.mash_dist)
    pairs = [
        (path_to_short[a], path_to_short[b], d)
        for a, b, d in raw_pairs
        if a in path_to_short and b in path_to_short
    ]
    nodes = set(path_to_short.values())
    truth = load_ground_truth_clusters(args.ground_truth)
    print(f"dereplication_stability: {len(nodes)} strains, {len(pairs)} pairwise distances, "
          f"{len(truth)}/{len(nodes)} with ground-truth cluster labels", file=sys.stderr)

    thresholds = [float(t) for t in args.thresholds.split(",")]
    rows = []
    for threshold in thresholds:
        sl_groups = single_linkage_groups(pairs, threshold, nodes)
        sl_ari = ari_against_ground_truth(groups_to_labels(sl_groups), truth)
        for method in ("complete", "average"):
            h_groups = hierarchical_groups(pairs, threshold, nodes, method)
            h_ari = ari_against_ground_truth(groups_to_labels(h_groups), truth)
            rows.append({
                "threshold": threshold,
                "method": method,
                "n_representative_groups": len(h_groups),
                "ari_vs_s21_clusters": "" if h_ari is None else f"{h_ari:.4f}",
            })
        rows.append({
            "threshold": threshold,
            "method": "single",
            "n_representative_groups": len(sl_groups),
            "ari_vs_s21_clusters": "" if sl_ari is None else f"{sl_ari:.4f}",
        })
        print(f"  threshold={threshold}: single={len(sl_groups)} groups (ARI={sl_ari}), "
              f"complete/average computed", file=sys.stderr)

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["threshold", "method", "n_representative_groups", "ari_vs_s21_clusters"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
