#!/usr/bin/env python3
"""PF11001 figures (16:9, >= 14 pt, Okabe-Ito colours).

fig1_island_fraction.{png,pdf}: per run, fraction of PF11001 shell+cloud
    families that are island members (with 95% Clopper-Pearson CI) vs the
    background fraction, n, p and q (pipeline annotation).
fig2_af293_chromosome_map.{png,pdf}: Af293 chromosomes, PF11001 genes,
    Table S7 Starship elements, and Af293 instances of islands that contain a
    PF11001 family.
"""
import csv
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import beta

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
D = AF / "analysis/pf11001"
FIG = D / "figures"
FIG.mkdir(parents=True, exist_ok=True)
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          sky="#56B4E9", purple="#CC79A7", grey="#999999", black="#000000")
plt.rcParams.update({"font.size": 15, "axes.titlesize": 17, "axes.labelsize": 16,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14})
LABEL = {"full_v072": "A. fumigatus\nfull_v072",
         "immitis_in_posadasii_out_v5": "C. immitis\nv5",
         "posadasii_in_immitis_out_v2": "C. posadasii\nv2"}


def cp_ci(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return lo, hi


def fig1():
    rows = []
    with open(D / "tables/island_enrichment_summary.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["annotation"].startswith("pipeline"):
                rows.append(r)
    fig, ax = plt.subplots(figsize=(16, 9))
    for i, r in enumerate(rows):
        a, K = int(r["a"]), int(r["K"])
        exp = float(r["expected_frac"])
        obs = a / K
        lo, hi = cp_ci(a, K)
        ax.bar(i - 0.18, exp, width=0.34, color=OI["grey"], label="all shell+cloud families (background)" if i == 0 else None)
        ax.bar(i + 0.18, obs, width=0.34, color=OI["blue"], label="PF11001 shell+cloud families" if i == 0 else None)
        ax.errorbar(i + 0.18, obs, yerr=[[obs - lo], [hi - obs]], color=OI["black"], capsize=8, lw=2)
        ax.text(i - 0.18, exp + 0.02, f"{int(r['n_island'])}/{int(r['n_bg'])}\n= {exp:.2f}", ha="center", va="bottom", fontsize=14)
        ax.text(i + 0.18, 1.04, f"{a}/{K} = {obs:.2f}", ha="center", va="bottom", fontsize=14)
        sig = "q < 0.05" if float(r["q"]) < 0.05 else "q >= 0.05 (not significant)"
        ax.text(i, -0.2, f"p = {float(r['p']):.1e}\nq = {float(r['q']):.3f}\n{sig}\nrank {r['pf11001_rank_by_p']} of {r['n_domains_tested']} domains",
                ha="center", va="top", fontsize=14)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([LABEL[r["run"]] for r in rows])
    ax.tick_params(axis="x", pad=8)
    ax.set_ylim(0, 1.3)
    ax.set_ylabel("fraction of families in >= 1 island")
    ax.set_title("PF11001 (AFUB_07903_YDR124W_hel): island membership of shell+cloud families\n"
                 "pipeline test: one-sided Fisher, BH-FDR over all domains tested; error bar = 95% Clopper-Pearson CI")
    ax.legend(loc="upper center", ncol=2, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.subplots_adjust(bottom=0.3, top=0.86, left=0.08, right=0.98)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"fig1_island_fraction.{ext}", dpi=150)
    plt.close(fig)


def fig2():
    genes = []
    with open(D / "tables/af293_pf11001_proteins.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.append(r)
    chrom = {}
    with open(AF / "data_dir/dna/UP000002530_330879.dna.fa") as fh:
        name = None
        n = 0
        for line in fh:
            if line.startswith(">"):
                if name:
                    chrom[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
        chrom[name] = n
    order = sorted(chrom)
    els = []
    with open(AF / "analysis/starship_island_overlap/inputs/af293_starships_table_s7.elements.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            els.append(r)
    # island instances in Af293 for islands that contain a PF11001 family
    pf_islands = set()
    with open(D / "tables/pf11001_islands.full_v072.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pf_islands.add(r["island_idx"])
    inst = []
    with gzip.open(AF / "analysis/starship_island_overlap/full_v072/per_island.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["island_idx"] in pf_islands and r["instances_in_annotated_strains"]:
                for x in r["instances_in_annotated_strains"].split(";"):
                    short, c, se = x.split(":")
                    s, e = se.split("-")
                    inst.append((c, int(s), int(e)))
    fig, ax = plt.subplots(figsize=(16, 9))
    for y, c in enumerate(order):
        yy = len(order) - 1 - y
        ax.plot([0, chrom[c] / 1e6], [yy, yy], color=OI["grey"], lw=6, solid_capstyle="round", zorder=1)
        ax.text(-0.15, yy, f"chr{y + 1}", ha="right", va="center", fontsize=15)
        for e in els:
            if e["contig"] == c:
                s, en = int(e["start"]) / 1e6, int(e["end"]) / 1e6
                ax.add_patch(Rectangle((s, yy - 0.32), max(en - s, 0.02), 0.64, color=OI["orange"], alpha=0.6, zorder=2))
                ax.text((s + en) / 2, yy - 0.36, e["label"], ha="center", va="top", fontsize=13, color="#8a5a00")
        for (ic, s, en) in inst:
            if ic == c:
                ax.add_patch(Rectangle((s / 1e6, yy - 0.18), max((en - s) / 1e6, 0.03), 0.36,
                                       facecolor="none", edgecolor=OI["green"], lw=2, zorder=3))
        for g in genes:
            if g.get("contig") == c and g.get("start"):
                x = int(g["start"]) / 1e6
                core = g.get("family_bin") in ("core", "soft_core")
                ax.scatter([x], [yy], marker="v" if core else "o", s=180, zorder=4,
                           color=OI["black"] if core else OI["blue"], edgecolor="white", lw=1)
                ax.text(x + 0.04, yy + 0.12, g.get("locus_tag") or g["uniprot"], ha="left", va="bottom", fontsize=13)
    unplaced = [g for g in genes if not g.get("start")]
    ax.scatter([], [], marker="o", s=180, color=OI["blue"], label="PF11001 gene, family not core")
    ax.scatter([], [], marker="v", s=180, color=OI["black"], label="PF11001 gene, core family")
    ax.add_patch(Rectangle((0, 0), 0, 0, color=OI["orange"], alpha=0.6, label="Starship element (Table S7)"))
    ax.add_patch(Rectangle((0, 0), 0, 0, facecolor="none", edgecolor=OI["green"], lw=2,
                           label="Af293 instance of an island with a PF11001 family"))
    ax.legend(loc="lower right", frameon=False, fontsize=14)
    ax.set_xlim(-0.9, max(chrom.values()) / 1e6 + 0.2)
    ax.set_ylim(-0.9, len(order) - 0.2)
    ax.set_yticks([])
    ax.set_xlabel("position (Mb)")
    t = f"Af293: {len(genes)} PF11001 proteins (hmmsearch --cut_ga)"
    if unplaced:
        t += f"; {len(unplaced)} without coordinates not drawn"
    ax.set_title(t)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.92, bottom=0.1)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"fig2_af293_chromosome_map.{ext}", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    fig1()
    fig2()
