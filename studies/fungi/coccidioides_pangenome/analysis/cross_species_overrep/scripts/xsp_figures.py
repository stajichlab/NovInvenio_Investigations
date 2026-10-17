#!/usr/bin/env python3
"""Slide figures (16:9, >= 14 pt) for the cross-species over-representation analysis.

fig1_tiers_A_B.png / fig1_tiers_C.png: heat plot of -log10 q per run for tier 1 and 2
  terms; cell text = k/K. Sequential single-hue blue; white text on dark cells.
fig2_upset_A.png: counts of island-enriched Pfam domains (definition A, q < 0.05)
  per combination of runs (UpSet style).
fig3_family_island.png: fraction of Af families whose Coccidioides RBH partner is an
  island family, Af island vs non-island, per Coccidioides run and pair set.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

OUT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/"
           "coccidioides_pangenome/analysis/cross_species_overrep")
FIG = OUT / "figures"
plt.rcParams.update({"font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15,
                     "font.family": "DejaVu Sans"})
BLUE = LinearSegmentedColormap.from_list(
    "blue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])
C_ISL, C_NON = "#2a78d6", "#eb6834"
RUNLAB = {"af": "A. fum.\nfull_v072", "v5": "C. imm.\nv5", "v2": "C. pos.\nv2"}
DEFLAB = {
    "A_island_pfam": "A: island (Pfam)", "A_island_clan": "A: island (clan)",
    "B_shell_soft_core_gain_pfam": "B: gain, shell+soft-core (Pfam)",
    "B_shell_soft_core_gain_clan": "B: gain, shell+soft-core (clan)",
    "B_all_accessory_gain_pfam": "B: gain, all accessory (Pfam)",
    "B_all_accessory_gain_clan": "B: gain, all accessory (clan)",
    "C_accessory_vs_all_pfam": "C: accessory vs all",
    "C_shellcloud_vs_core_pfam": "C: shell+cloud vs core",
}


def heat(ax, d, title, vmax):
    runs = ["af", "v5", "v2"]
    m = np.array([[-math.log10(max(r[f"{k}_q"], 1e-300)) if pd.notna(r[f"{k}_q"]) else np.nan
                   for k in runs] for _, r in d.iterrows()])
    ax.imshow(np.clip(m, 0, vmax), cmap=BLUE, vmin=0, vmax=vmax, aspect="auto")
    for i, (_, r) in enumerate(d.iterrows()):
        for j, k in enumerate(runs):
            if pd.isna(r[f"{k}_q"]):
                txt = "n.t."
            else:
                txt = f"{int(r[f'{k}_k'])}/{int(r[f'{k}_K'])}"
            v = m[i, j]
            col = "white" if (not np.isnan(v) and v > vmax * 0.45) else "#1a1a19"
            star = "" if (np.isnan(v) or r[f"{k}_q"] >= 0.05) else " *"
            ax.text(j, i, txt + star, ha="center", va="center", fontsize=13, color=col)
    ax.set_xticks(range(3), [RUNLAB[k] for k in runs], fontsize=13)
    ax.xaxis.tick_top()
    labels = [f"T{r['tier']}  {r['name']}".replace("AFUB_07903_YDR124W_hel", "PF11001 YDR124W_hel")
              for _, r in d.iterrows()]
    ax.set_yticks(range(len(d)), labels, fontsize=14)
    ax.set_title(title, loc="left", pad=48)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, 3), minor=True)
    ax.set_yticks(np.arange(-0.5, len(d)), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)


def fig1():
    t = pd.read_csv(OUT / "tiers.tsv", sep="\t")
    t = t[t["tier"].isin([1, 2])].copy()
    vmax = 12
    groups = {
        "fig1_tiers_A_B": [["A_island_pfam", "A_island_clan"],
                           ["B_shell_soft_core_gain_pfam", "B_shell_soft_core_gain_clan",
                            "B_all_accessory_gain_pfam", "B_all_accessory_gain_clan"]],
        "fig1_tiers_C": [["C_accessory_vs_all_pfam"], ["C_shellcloud_vs_core_pfam"]],
    }
    for fname, panels in groups.items():
        fig = plt.figure(figsize=(16, 9))
        axes = [fig.add_axes([0.235, 0.09, 0.2, 0.71]), fig.add_axes([0.685, 0.09, 0.2, 0.71])]
        cax = fig.add_axes([0.925, 0.2, 0.012, 0.5])
        for ax, defs in zip(axes, panels):
            parts = []
            for dfn in defs:
                s = t[t["definition"] == dfn].sort_values(["tier", "af_q"]).copy()
                if fname == "fig1_tiers_A_B":
                    scope = ("ssc " if "shell_soft" in dfn
                             else "acc " if "all_acc" in dfn else "")
                    kind = "clan" if dfn.endswith("clan") else "Pfam"
                    s["name"] = s["name"] + f"  [{scope}{kind}]"
                parts.append(s)
            d = pd.concat(parts)
            title = " / ".join(DEFLAB[x].split(" (")[0] for x in defs[:1]) if fname == "fig1_tiers_C" \
                else ("A: islands vs shell+cloud" if defs[0].startswith("A") else "B: gained families")
            heat(ax, d.reset_index(drop=True), title, vmax)
        sm = plt.cm.ScalarMappable(cmap=BLUE, norm=plt.Normalize(0, vmax))
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(f"-log10 q (capped at {vmax})")
        fig.text(0.01, 0.01, "Cell: k/K families with term (k in test set, K in universe); "
                 "* q < 0.05; n.t. not tested.\nT1: q < 0.05 in all 3 runs; T2: Af + one Coccidioides run. "
                 "ssc = shell+soft-core scope; acc = all accessory scope.",
                 fontsize=13, color="#52514d")
        fig.savefig(FIG / f"{fname}.png", dpi=150)
        fig.savefig(FIG / f"{fname}.pdf")
        plt.close(fig)


def fig2():
    d = pd.read_csv(OUT / "domains_A_island_pfam.tsv", sep="\t")
    # include terms tested in only one run too: read per-run sig sets from all rows
    sig = {k: set(d.loc[d[f"{k}_q"] < 0.05, "term"]) for k in ("af", "v5", "v2")}
    # terms tested in only one run are not in domains_*.tsv; count them separately
    import itertools
    combos = []
    for r in (3, 2, 1):
        for c in itertools.combinations(("af", "v5", "v2"), r):
            combos.append(c)
    allsig = set().union(*sig.values())
    counts = []
    for c in combos:
        n = sum(all(t in sig[k] for k in c) and not any(t in sig[k] for k in ("af", "v5", "v2")
                                                          if k not in c) for t in allsig)
        counts.append(n)
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.25, 0.42, 0.7, 0.48])
    am = fig.add_axes([0.25, 0.1, 0.7, 0.28], sharex=ax)
    x = np.arange(len(combos))
    ax.bar(x, counts, width=0.6, color="#2a78d6", edgecolor="white", linewidth=2)
    for xi, n in zip(x, counts):
        ax.text(xi, n + 0.2, str(n), ha="center", va="bottom", fontsize=15)
    ax.set_ylabel("Pfam domains with q < 0.05")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0, labelbottom=False)
    ax.set_title("Definition A: island-enriched Pfam domains, overlap between runs "
                 "(domains tested in >= 2 runs)", loc="left")
    names = ["af", "v5", "v2"]
    for xi, c in zip(x, combos):
        ys = [2 - names.index(k) for k in c]
        am.scatter([xi] * 3, [0, 1, 2], s=140, color="#d6d5cf", zorder=1)
        am.scatter([xi] * len(ys), ys, s=160, color="#1a1a19", zorder=2)
        if len(ys) > 1:
            am.plot([xi, xi], [min(ys), max(ys)], color="#1a1a19", linewidth=2.5, zorder=1)
    summ = pd.read_csv(OUT / "summary_domains.tsv", sep="\t").set_index("definition").loc["A_island_pfam"]
    am.set_yticks([2, 1, 0], [RUNLAB[k].replace("\n", " ") + f" ({len(sig[k])} of {int(summ[k + '_n_q_lt_0.05'])})"
                              for k in names])
    ax.set_ylim(0, max(counts) * 1.15)
    am.set_ylim(-0.6, 2.6)
    am.tick_params(length=0)
    am.set_xticks([])
    for s in am.spines.values():
        s.set_visible(False)
    fig.text(0.01, 0.02, "Counts are exclusive intersections. (x of y): x significant domains tested in >= 2 runs, "
             "y significant domains in that run in total.", fontsize=13, color="#52514d")
    fig.savefig(FIG / "fig2_upset_A.png", dpi=150)
    fig.savefig(FIG / "fig2_upset_A.pdf")
    plt.close(fig)


def fig3():
    t = pd.read_csv(OUT / "family_tests.tsv", sep="\t")
    t = t[t["test"].isin(["island|all_rbh_pairs", "island|both_shell_or_cloud"])]
    fig, ax = plt.subplots(figsize=(16, 9))
    labs, xi, w = [], 0, 0.36
    for _, r in t.iterrows():
        run = "C. immitis v5" if "immitis_in" in r["cocci_run"] else "C. posadasii v2"
        ps = "all RBH pairs" if r["test"].endswith("all_rbh_pairs") else "both shell/cloud"
        n1 = r["af_yes_partner_yes"] + r["af_yes_partner_no"]
        n0 = r["af_no_partner_yes"] + r["af_no_partner_no"]
        f1, f0 = r["frac_partner_yes_given_af_yes"], r["frac_partner_yes_given_af_no"]
        ax.bar(xi - w / 2 - 0.01, f1, w, color=C_ISL, label="Af island family" if xi == 0 else None)
        ax.bar(xi + w / 2 + 0.01, f0, w, color=C_NON, label="Af non-island family" if xi == 0 else None)
        ax.text(xi - w / 2, f1 + 0.01, f"{f1:.2f}\n{r['af_yes_partner_yes']}/{n1}", ha="center",
                va="bottom", fontsize=13)
        ax.text(xi + w / 2, f0 + 0.01, f"{f0:.2f}\n{r['af_no_partner_yes']}/{n0}", ha="center",
                va="bottom", fontsize=13)
        labs.append(f"{run}\n{ps}\nOR {r['odds_ratio']} p {float(r['p_two_sided']):.1e}")
        xi += 1
    ax.set_xticks(range(len(labs)), labs, fontsize=13)
    ax.set_ylabel("Fraction with RBH partner in a\nCoccidioides island family")
    ax.set_ylim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("A. fumigatus families with a Coccidioides reciprocal best hit: "
                 "is the partner an island family?", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_family_island.png", dpi=150)
    fig.savefig(FIG / "fig3_family_island.pdf")
    plt.close(fig)


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["1", "2", "3"]
    if "1" in which:
        fig1()
    if "2" in which:
        fig2()
    if "3" in which:
        fig3()
