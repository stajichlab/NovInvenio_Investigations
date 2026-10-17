#!/usr/bin/env python3
"""Slide figure (16:9, >= 14 pt): Af293 antiSMASH SM clusters per type, with vs without
Starship overlap (rule A: >= 1 cluster gene midpoint inside a Table S7 element).
Reads <dir>/sm_by_type.tsv and <dir>/sm_permutation.tsv; writes <prefix>.png/.svg/.pdf.
Colours: Okabe-Ito vermillion (#D55E00) for overlap, blue (#0072B2) for no overlap.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OVER, NOOV, INK, INK2 = "#D55E00", "#0072B2", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 15, "axes.titlesize": 18, "axes.labelsize": 15,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "axes.edgecolor": INK2, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2})


def tsv(p: Path):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    a = ap.parse_args()
    rows = tsv(a.dir / "sm_by_type.tsv")
    perm = {r["statistic"]: r for r in tsv(a.dir / "sm_permutation.tsv")}
    allrow = next(r for r in rows if r["type"] == "ALL")
    rows = [r for r in rows if r["type"] != "ALL"][::-1]
    labels = [r["type"] for r in rows]
    ov = [int(r["n_overlap_ge1gene"]) for r in rows]
    tot = [int(r["n_clusters"]) for r in rows]
    no = [t - o for t, o in zip(tot, ov)]
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    y = range(len(rows))
    ax.barh(y, no, color=NOOV, label="No gene inside a Starship")
    ax.barh(y, ov, left=no, color=OVER, label=">= 1 gene inside a Starship")
    xmax = max(tot)
    for i, (o, t, r) in enumerate(zip(ov, tot, rows)):
        b = int(r["n_overlap_ge50pct"])
        ax.text(t + xmax * 0.015, i, f"{o}/{t} ({100 * o / t:.0f}%)" + (f"; >= 50%: {b}" if b else ""),
                va="center", fontsize=14, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.set_xlim(0, xmax * 1.45)
    ax.set_xlabel("antiSMASH regions (SM clusters) in Af293")
    pa = perm["clusters_overlap_ge1gene"]
    n_all, o_all, b_all = int(allrow["n_clusters"]), int(allrow["n_overlap_ge1gene"]), int(allrow["n_overlap_ge50pct"])
    fig.suptitle(f"A. fumigatus Af293: {o_all} of {n_all} SM clusters overlap a Starship\n"
                 f"(>= 1 gene inside: {o_all}; >= 50% of genes inside: {b_all})",
                 x=0.01, ha="left", fontsize=18)
    ax.text(0.99, 0.04, f"Element-shift permutation (n = {pa['n_perm']}): observed {pa['observed']}, "
                        f"null mean {float(pa['null_mean']):.2f}, p(>=) = {float(pa['p_perm_greater']):.3g}\n"
                        "antiSMASH 8.0.4 regions; Starships: Gluck-Thaler 2025 Table S7 (8 elements)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=14, color=INK2)
    ax.legend(loc="center right", frameon=False)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{a.out_prefix}.{ext}", dpi=150 if ext == "png" else None)
    print(f"wrote {a.out_prefix}.png/.svg/.pdf")


if __name__ == "__main__":
    main()
