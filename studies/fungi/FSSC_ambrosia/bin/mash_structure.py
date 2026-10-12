#!/usr/bin/env python3
"""Describe the mash structure of a study's selected strains.

Usage: <NI pixi python> mash_structure.py --study-dir <study>
Reads selection/selected_strains.tsv and selection/selected_mash_matrix.tsv
(bin/select_fssc_subclade.py). Writes selection/ingroup_mash_clades.tsv:
UPGMA (average linkage) clusters of the ingroup cut at several heights, and
prints within/between summaries. Needs scipy.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
from scipy.spatial.distance import squareform

CUTS = (0.005, 0.01, 0.02, 0.04)

ap = argparse.ArgumentParser()
ap.add_argument("--study-dir", required=True)
sel = Path(ap.parse_args().study_dir) / "selection"
rows = {r["asmid"]: r for r in csv.DictReader(open(sel / "selected_strains.tsv"), delimiter="\t")}
with open(sel / "selected_mash_matrix.tsv") as fh:
    ids = fh.readline().rstrip("\n").split("\t")[1:]
    M = np.array([[float(x) for x in line.rstrip("\n").split("\t")[1:]] for line in fh])
ing = [i for i, a in enumerate(ids) if rows[a]["GROUP"] == "IN"]
out = [i for i, a in enumerate(ids) if rows[a]["GROUP"] == "OUT"]
Mi = M[np.ix_(ing, ing)]
Z = linkage(squareform(Mi, checks=False), method="average")
labs = {c: fcluster(Z, t=c, criterion="distance") for c in CUTS}
order = dendrogram(Z, no_plot=True)["leaves"]
with open(sel / "ingroup_mash_clades.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["dendrogram_order", "asmid", "species", "strain"] + [f"upgma_cut_{c}" for c in CUTS])
    for k, j in enumerate(order):
        a = ids[ing[j]]
        w.writerow([k, a, rows[a]["species"], rows[a]["strain"]] + [int(labs[c][j]) for c in CUTS])
iu = np.triu_indices(len(ing), 1)
print(f"ingroup n={len(ing)} pairwise mash: min {Mi[iu].min():.4f} median {np.median(Mi[iu]):.4f} max {Mi[iu].max():.4f}")
for c in CUTS:
    print(f"UPGMA cut {c}: {len(set(labs[c]))} clusters, sizes {sorted(np.bincount(labs[c])[1:].tolist(), reverse=True)}")
for i in out:
    di = M[i, ing]
    print(f"OUT {rows[ids[i]]['species']} {rows[ids[i]]['strain']}: to ingroup min {di.min():.4f} median {np.median(di):.4f} max {di.max():.4f}")
print("dendrogram order (cut 0.02 / 0.04):")
for k, j in enumerate(order):
    a = ids[ing[j]]
    print(f"  {labs[0.02][j]:2d} {labs[0.04][j]:2d}  {rows[a]['species']} {rows[a]['strain']}")
