#!/usr/bin/env python3
"""Slide figures (16:9, >= 14 pt) from bin/starship_island_overlap.py outputs.

  --overlap-dir   <run> output dir (enrichment.tsv, per_starship.tsv, sm_*.tsv)
  --out-prefix    figures/<run>  -> <prefix>_starship_islands.png, <prefix>_sm.png
  --title         run label for the figure title
Colours: reference categorical palette (blue #2a78d6, orange #eb6834) + grey baseline.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#9a9893", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
NAMES = {"starship_all": "Starship families\n(Af293, 8 elements)",
         "starship_high": "Starship families\n(Af293, 7 high-conf.)",
         "captain": "Captain families\n(DUF3435 hit)"}


def tsv(p: Path):
    return list(csv.DictReader(open(p), delimiter="\t"))


def panel_fraction(ax, enr, universe):
    rows = [r for r in enr if r["universe"] == universe]
    labels, vals, txt, cols = [], [], [], []
    base = rows[0]
    labels.append("All non-core\nfamilies (baseline)" if universe == "noncore" else "Universe (baseline)")
    vals.append(float(base["frac_universe_in_island"]))
    txt.append(f"{int(base['n_island']):,} / {int(base['N_universe']):,}")
    cols.append(GREY)
    for r in rows:
        labels.append(NAMES.get(r["annotation"], r["annotation"]))
        k, K = int(r["k_both"]), int(r["K_annotated"])
        vals.append(k / K if K else 0)
        txt.append(f"{k:,} / {K:,}")
        cols.append(BLUE)
    y = range(len(labels))[::-1]
    ax.barh(list(y), vals, color=cols, height=0.6)
    for yi, v, t in zip(y, vals, txt):
        ax.text(v + 0.01, yi, f"{100 * v:.1f}%  ({t})", va="center", fontsize=14, color=INK)
    ax.set_yticks(list(y), labels)
    ax.set_xlim(0, 1.15)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Families in >= 1 island (non-core families only)")
    ax.set_title("Which families fall in pangenome islands", loc="left")


def panel_per_starship(ax, ships):
    ships = sorted(ships, key=lambda r: r["element_id"])
    labels = [f'{r["label"].replace(" (likely)", "?")}' for r in ships]
    fr = [float(r["frac_families_in_island"] or 0) for r in ships]
    y = list(range(len(ships)))[::-1]
    for yi, r, v in zip(y, ships, fr):
        hi = r["confidence"] == "high"
        ax.barh(yi, v, color=ORANGE if hi else "white", edgecolor=ORANGE, hatch=None if hi else "//", height=0.6)
        ax.text(v + 0.01, yi, f'{r["n_families_in_island"]}/{r["n_families"]}', va="center", fontsize=14)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1.15)
    ax.set_xticks([0, 0.5, 1.0], ["0", "50%", "100%"])
    ax.set_xlabel("In >= 1 island (all bins)")
    ax.set_title("Per Af293 Starship\n(hatched: unknown conf.)", loc="left", fontsize=15)


def fig_islands(d: Path, out: Path, title: str):
    enr = tsv(d / "enrichment.tsv")
    has_ships = (d / "per_starship.tsv").exists()
    fig = plt.figure(figsize=(16, 9))
    if has_ships:
        ax1 = fig.add_axes([0.20, 0.16, 0.40, 0.66])
        ax2 = fig.add_axes([0.78, 0.16, 0.17, 0.66])
        panel_per_starship(ax2, tsv(d / "per_starship.tsv"))
    else:
        ax1 = fig.add_axes([0.22, 0.14, 0.68, 0.70])
    panel_fraction(ax1, enr, "noncore")
    fig.suptitle(title, x=0.02, ha="left", fontsize=20, fontweight="bold")
    short = {"starship_all": "Starship (8)", "starship_high": "Starship (7 high)", "captain": "captain"}
    ps = "; ".join(f'{short.get(r["annotation"], r["annotation"])} p = {r["p_hypergeom_greater"]}'
                   for r in enr if r["universe"] == "noncore")
    fig.text(0.02, 0.055, "Bars: share of each non-core family set in >= 1 island (significant_islands.tsv); grey = all non-core families. "
             "Labels: in island / total.", fontsize=13, color=INK2)
    fig.text(0.02, 0.02, "Hypergeometric, universe = non-core families: " + ps, fontsize=13, color=INK2)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_sm_cocci(d: Path, out: Path, title: str):
    rows = [r for r in tsv(d / "sm_by_type.tsv") if r["species"] == "all"]
    rows = sorted(rows, key=lambda r: -int(r["n_clusters"]))
    top = rows[:10]
    other = rows[10:]
    if other:
        top.append({"sm_type": f"other ({len(other)} types)", "n_clusters": sum(int(r["n_clusters"]) for r in other),
                    "n_with_captain_gene": sum(int(r["n_with_captain_gene"]) for r in other)})
    fig, ax = plt.subplots(figsize=(16, 9))
    fig.subplots_adjust(left=0.25, right=0.95, top=0.85, bottom=0.14)
    y = list(range(len(top)))[::-1]
    for yi, r in zip(y, top):
        n, c = int(r["n_clusters"]), int(r["n_with_captain_gene"])
        ax.barh(yi, n - c, color=GREY, height=0.6)
        ax.barh(yi, c, left=n - c, color=ORANGE, height=0.6)
        ax.text(n * 1.01 + 5, yi, f"{c} of {n:,} with captain gene ({100 * c / n:.1f}%)", va="center", fontsize=14)
    ax.set_yticks(y, [r["sm_type"] for r in top])
    ax.set_xlim(0, max(int(r["n_clusters"]) for r in top) * 1.55)
    ax.set_xlabel("antiSMASH SM regions (all strains with antiSMASH output)")
    fig.suptitle(title, x=0.02, ha="left", fontsize=20, fontweight="bold")
    ax.set_title("Orange = region contains a DUF3435 (captain-type) gene. No Starship calls exist for Coccidioides.",
                 loc="left", fontsize=15, color=INK2)
    m = {r["metric"]: r["value"] for r in tsv(d / "sm_starship_overlap.tsv")}
    fig.text(0.02, 0.03, f'Genes: {m["captain_genes_in_sm_clusters"]} captain genes in SM regions, expected '
             f'{m["expected_captain_genes_in_sm_clusters"]} (hypergeometric depletion p = {m["p_hypergeom_less"]}; '
             f'{m["strains_with_antismash"]} strains, {int(m["genes_total"]):,} genes; antiSMASH 8.0.1 annotation).',
             fontsize=13, color=INK2)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_sm_af(d: Path, out: Path, title: str):
    import ast
    rows = tsv(d / "sm_by_type.tsv")
    met = {r["metric"]: r["value"] for r in tsv(d / "sm_starship_overlap.tsv")}
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.76, bottom=0.20, wspace=0.40)
    labels = [r["sm_type"] for r in rows]
    y = list(range(len(rows)))[::-1]
    ax = axes[0]
    for yi, r in zip(y, rows):
        n, c = int(r["af293_genes"]), int(r["af293_genes_in_starship"])
        ax.barh(yi, n - c, color=GREY, height=0.45)
        ax.barh(yi, c, left=n - c, color=ORANGE, height=0.45)
        ax.text(n + 0.4, yi, f"{c} of {n} inside a Starship", va="center", fontsize=14)
    ax.set_yticks(y, labels)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(0, max(int(r["af293_genes"]) for r in rows) * 2.0)
    ax.set_xlabel("Af293 backbone genes")
    ax.set_title("Af293 genes vs Table S7 elements", loc="left")
    ax = axes[1]
    h = 0.25
    for yi, r in zip(y, rows):
        n, c, i = int(r["families_all_strains"]), int(r["families_starship"]), int(r["families_in_island"])
        for off, v, col, txt in ((h, n, GREY, f"{n} families"), (0, i, BLUE, f"{i} in >= 1 island"),
                                 (-h, c, ORANGE, f"{c} Starship families")):
            ax.barh(yi + off, v, color=col, height=h * 0.9)
            ax.text(v + 3, yi + off, txt, va="center", fontsize=14)
    ax.set_yticks(y, labels)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(0, max(int(r["families_all_strains"]) for r in rows) * 1.7)
    ax.set_xlabel("Pangenome families with a backbone hit")
    ax.set_title("Families (all 295 strains)", loc="left")
    fig.suptitle(title, x=0.02, ha="left", fontsize=20, fontweight="bold")
    fig.text(0.02, 0.86, "SM proxy: hmmsearch PF00109 (PKS KS) + PF00668 (NRPS C), E <= 1e-5. "
             "No antiSMASH clusters exist for this annotation.", fontsize=14, color=INK2)
    hg = ast.literal_eval(met["hypergeom_af293_genes"])
    hf = ast.literal_eval(met["hypergeom_family_island_noncore"])
    fig.text(0.02, 0.07, f"Af293: {hg['k']} of {hg['n']} backbone genes in Starships, expected {hg['expected']} "
             f"(hypergeometric p = {hg['p_greater']}, universe {hg['N']:,} Af293 genes).", fontsize=14, color=INK2)
    fig.text(0.02, 0.025, f"Families: {hf['k']} of {hf['K']} non-core backbone families in islands, expected {hf['expected']} "
             f"(p = {hf['p_greater']}, universe {hf['N']:,} non-core families).", fontsize=14, color=INK2)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlap-dir", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument("--title", required=True)
    ap.add_argument("--sm", choices=["none", "af", "cocci"], default="none")
    a = ap.parse_args()
    a.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    suffix = (": Starship families and pangenome islands" if (a.overlap_dir / "per_starship.tsv").exists()
              else ": captain (DUF3435) families and islands (no Starship calls)")
    fig_islands(a.overlap_dir, Path(f"{a.out_prefix}_starship_islands.png"), a.title + suffix)
    if a.sm == "af":
        fig_sm_af(a.overlap_dir, Path(f"{a.out_prefix}_sm.png"), a.title + ": SM backbone genes and Starships")
    elif a.sm == "cocci":
        fig_sm_cocci(a.overlap_dir, Path(f"{a.out_prefix}_sm.png"), a.title + ": SM regions and captain genes")


if __name__ == "__main__":
    main()
