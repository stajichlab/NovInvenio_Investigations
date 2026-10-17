#!/usr/bin/env python3
"""Figures for the shared-domain context analysis (16:9, >= 14 pt, colourblind-safe).

fig1_island_cooccurrence_<count_type>.png/.pdf  log2(observed/expected) per label pair, Af / v5 / v2
fig2_island_tier1_counts.png/.pdf               islands with 1-4 tier-1 labels per run
fig3_duf3435_tree.png/.pdf                       DUF3435 domain tree, tips coloured by group
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
D = NII / "studies/fungi/coccidioides_pangenome/analysis/shared_domains_context"
T = D / "tables"
F = D / "figures"
LAB6 = ["APH", "DUF3435", "DUF3723", "DUF3723_N", "PF11001", "HrmA"]
RUNS = [("af", "A. fumigatus full_v072"), ("v5", "C. immitis v5"), ("v2", "C. posadasii v2")]
# Okabe-Ito
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00",
      "sky": "#56B4E9", "purple": "#CC79A7", "grey": "#999999", "black": "#000000"}
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14})


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fig1(count_type, title):
    rows = [r for r in read(T / "island_pair_cooccurrence.tsv") if r["count_type"] == count_type]
    fig, axes = plt.subplots(1, 3, figsize=(16, 9), constrained_layout=True)
    vmax = 5
    for ax, (tag, name) in zip(axes, RUNS):
        M = np.full((6, 6), np.nan)
        A = [["" for _ in range(6)] for _ in range(6)]
        for r in rows:
            if r["tag"] != tag or r["label_a"] not in LAB6 or r["label_b"] not in LAB6:
                continue
            i, j = LAB6.index(r["label_a"]), LAB6.index(r["label_b"])
            o, e = float(r["observed_islands"]), float(r["perm_mean"])
            v = np.log2((o + 0.5) / (e + 0.5))
            star = "*" if r.get("q_upper_bh") and float(r["q_upper_bh"]) < 0.05 else ""
            for (a, b) in ((i, j), (j, i)):
                M[a, b] = v
                A[a][b] = f"{int(o)}\n{e:.0f}{star}"
        im = ax.imshow(M, cmap="PuOr_r", vmin=-vmax, vmax=vmax)
        for a in range(6):
            for b in range(6):
                if A[a][b]:
                    ax.text(b, a, A[a][b], ha="center", va="center", fontsize=14,
                            color="white" if abs(M[a, b]) > 3 else "black")
        ax.set_xticks(range(6))
        ax.set_xticklabels(LAB6, rotation=60, ha="right")
        ax.set_yticks(range(6))
        ax.set_yticklabels(LAB6 if tag == "af" else [""] * 6)
        ax.set_title(name)
    cb = fig.colorbar(im, ax=axes, shrink=0.6, location="bottom")
    cb.set_label("log2((observed + 0.5) / (permutation mean + 0.5))")
    fig.suptitle(title, fontsize=16)
    cb.ax.set_title("Cell text: observed / permutation mean (1,000 label shuffles among island-resident "
             "families). * = BH q < 0.05 (upper tail) within run.", fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(F / f"fig1_island_cooccurrence_{count_type}.{ext}", dpi=150)
    plt.close(fig)


def fig2():
    rows = read(T / "island_tier1_count_summary.tsv")
    fig, ax = plt.subplots(figsize=(16, 9), constrained_layout=True)
    w = 0.26
    cols = [OI["blue"], OI["orange"], OI["green"]]
    for k, ((tag, name), c) in enumerate(zip(RUNS, cols)):
        v = [int(r["n_islands"]) for r in rows if r["tag"] == tag]
        tot = [int(r["n_islands_total"]) for r in rows if r["tag"] == tag][0]
        x = np.arange(1, 5) + (k - 1) * w
        ax.bar(x, v, width=w, color=c, label=f"{name} ({tot:,} islands)")
        for xi, vi in zip(x, v):
            ax.text(xi, vi + 20, f"{vi:,}", ha="center", va="bottom", fontsize=13)
    ax.set_xticks(range(1, 5))
    ax.set_xlabel("Number of tier-1 labels in the island (APH, DUF3435, DUF3723, DUF3723_N)")
    ax.set_ylabel("Islands")
    ax.set_title("Islands with >= 1 family carrying a tier-1 domain (family representative, hmmsearch --cut_ga)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(F / f"fig2_island_tier1_counts.{ext}", dpi=150)
    plt.close(fig)


def fig3():
    tips = {r["tip"]: r for r in read(D / "tree/duf3435_tips.tsv")}
    tree = Phylo.read(D / "tree/duf3435.fasttree.midpoint.nwk", "newick")
    col = {"A. fumigatus Af293": OI["vermillion"], "Coccidioides": OI["blue"], "other Dikarya": OI["grey"]}
    term = tree.get_terminals()
    ypos = {id(t): i for i, t in enumerate(term)}
    depth = tree.depths()
    if not max(depth.values()):
        depth = tree.depths(unit_branch_lengths=True)

    def y(cl):
        if id(cl) not in ypos:
            ypos[id(cl)] = np.mean([y(c) for c in cl.clades])
        return ypos[id(cl)]

    fig, ax = plt.subplots(figsize=(16, 9), constrained_layout=True)
    for cl in tree.find_clades(order="preorder"):
        x0 = depth[cl]
        for c in cl.clades:
            ax.plot([x0, x0, depth[c]], [y(cl), y(c), y(c)], color="black", lw=0.9)
        if not cl.is_terminal() and cl.confidence is not None and cl.confidence >= 0.9 and cl.clades:
            ax.plot(x0, y(cl), "o", ms=4, color=OI["black"])
    xmax = max(depth.values())
    for t in term:
        g = tips[t.name]["group"]
        ax.plot(depth[t], y(t), "s", ms=8 if g != "other Dikarya" else 6, color=col[g])
    for g, c in col.items():
        n = sum(1 for r in tips.values() if r["group"] == g)
        ax.plot([], [], "s", color=c, label=f"{g}\n({n} domains)")
    ax.plot([], [], "o", color="black", ms=4, label="node SH-like\nsupport >= 0.9")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_yticks([])
    ax.set_xlabel("substitutions per site (FastTree -lg -gamma, midpoint root)")
    ax.set_title("DUF3435 (PF11917) domain envelopes: Af293, RS, Silveira and 12 other Dikarya proteomes")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlim(-0.02 * xmax, xmax * 1.03)
    for ext in ("png", "pdf"):
        fig.savefig(F / f"fig3_duf3435_tree.{ext}", dpi=150)
    plt.close(fig)


def main():
    F.mkdir(exist_ok=True)
    fig1("distinct_families", "Island co-occurrence of domain labels: islands with the two labels on different families")
    fig1("family_pairs", "Island co-occurrence of domain labels: family pairs that share >= 1 island")
    fig2()
    fig3()


if __name__ == "__main__":
    main()
