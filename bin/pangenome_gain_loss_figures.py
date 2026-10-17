#!/usr/bin/env python3
"""Slide figures (16:9) from one or more bin/pangenome_gain_loss.py output dirs.

Usage:
  bin/pangenome_gain_loss_figures.py --out-dir <dir> \\
      --study "A. fumigatus=studies/fungi/Afumigatus_pangenome/analysis/gain_loss/full_v072" \\
      --study "C. immitis=studies/fungi/coccidioides_pangenome/analysis/gain_loss/immitis_in_posadasii_out_v5"

Writes fig1_gain_loss_by_study, fig2_gain_loss_by_class, fig3_size_distribution
(.png, .svg, .pdf) and figure_data.tsv (the numbers drawn). Counts only; no
biological interpretation.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import yaml  # noqa: E402

# Categorical slots 1-3 of the dataviz reference palette (validated: CVD dE >= 9.2).
COLORS = {"gain": "#2a78d6", "loss": "#eb6834", "ambiguous": "#1baf7a"}
DIRS = ("gain", "loss", "ambiguous")
CLASSES = ("soft_core", "shell", "cloud", "singleton")
SIZE_ORDER = ("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 16, "xtick.labelsize": 14,
    "ytick.labelsize": 14, "legend.fontsize": 14, "figure.titlesize": 18,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "svg.hashsalt": "gain_loss",
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def read_tsv(p: Path) -> list[dict]:
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(study_dir: Path) -> dict:
    prov = yaml.safe_load(open(study_dir / "PROVENANCE.yaml"))
    scopes = {r["scope"]: r for r in read_tsv(study_dir / "summary_scopes.tsv")}
    sel = next(k for k in scopes if k.startswith("selected:"))
    bycls = {r["class"]: r for r in read_tsv(study_dir / "summary_direction_by_class.tsv")}
    size = read_tsv(study_dir / "summary_size_distribution.tsv")
    return {"prov": prov, "scopes": scopes, "sel": sel, "bycls": bycls, "size": size}


def subtitle(name: str, d: dict) -> str:
    n = d["prov"]["strain_counts"]
    return (f"{name}: {n['in_rep'] + n['in_nonrep']} ingroup genomes "
            f"({n['in_rep']} representatives), outgroup {n['out_all']} genomes")


def rule_text(d: dict) -> str:
    p = d["prov"]["parameters"]
    return (f"outgroup-frequency rule: loss if >= {p['loss_min_frac']:g} of outgroup genomes carry the "
            f"family, gain if <= {p['gain_max_frac']:g}, else ambiguous")


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def stacked(ax, xs, values: dict, width, label_min_frac=0.04, total_labels=True):
    """values: dir -> list per x. Stacked bars with a 2px surface gap and direct labels."""
    bottoms = [0.0] * len(xs)
    tops = [sum(values[d][i] for d in DIRS) for i in range(len(xs))]
    ymax = max(tops) if tops else 1
    for d in DIRS:
        v = values[d]
        ax.bar(xs, v, width, bottom=bottoms, color=COLORS[d], edgecolor=SURF, linewidth=2, label=d)
        for i, x in enumerate(xs):
            if v[i] and v[i] / ymax >= label_min_frac:
                ax.text(x, bottoms[i] + v[i] / 2, f"{v[i]:,}", ha="center", va="center",
                        color="white", fontsize=13, fontweight="bold")
            elif v[i]:
                # small segment: label beside the bar in text ink, never in series colour
                off = {"gain": -0.02, "loss": 0.0, "ambiguous": 0.02}[d] * ymax
                ax.text(x + width / 2 + 0.03, bottoms[i] + v[i] / 2 + off, f"{d} {v[i]:,}",
                        ha="left", va="center", color=INK, fontsize=12)
        bottoms = [b + y for b, y in zip(bottoms, v)]
    ax.set_ylim(0, ymax * 1.1)
    if total_labels:
        for x, t in zip(xs, tops):
            ax.text(x, t + ymax * 0.01, f"n={t:,}", ha="center", va="bottom", color=INK, fontsize=13)


def save(fig, out: Path, stem: str):
    for ext in ("png", "svg", "pdf"):
        kw = {"dpi": 200} if ext == "png" else {}
        md = {"Date": None} if ext == "svg" else ({"CreationDate": None} if ext == "pdf" else {})
        fig.savefig(out / f"{stem}.{ext}", facecolor=SURF, metadata=md, **kw)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", action="append", required=True, help="NAME=gain_loss output dir")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args(argv)
    studies = []
    for s in args.study:
        name, path = s.split("=", 1)
        studies.append((name, load(Path(path))))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data_rows = []
    rule = rule_text(studies[0][1])
    subs = "\n".join(subtitle(n, d) for n, d in studies)

    # fig 1: all accessory vs pair-only, per study
    fig, ax = plt.subplots(figsize=(13.333, 7.5))
    xs, labels, vals = [], [], {d: [] for d in DIRS}
    for i, (name, d) in enumerate(studies):
        for j, (scope, lab) in enumerate(((d["sel"], "all accessory\nfamilies"),
                                          ("pair_family_a", "co-occurring-pair\nfamilies only"))):
            r = d["scopes"][scope]
            xs.append(i * 2.6 + j * 1.05)
            labels.append(f"{lab}")
            for dd in DIRS:
                vals[dd].append(int(r[f"freq_{dd}"]))
                data_rows.append(["fig1", name, scope, dd, r[f"freq_{dd}"]])
    stacked(ax, xs, vals, 0.9)
    ax.set_xlim(xs[0] - 0.6, xs[-1] + 1.1)
    ax.set_xticks(xs, labels)
    for i, (name, _) in enumerate(studies):
        ax.text(i * 2.6 + 0.525, -0.17, name, transform=ax.get_xaxis_transform(), ha="center",
                va="top", fontsize=16, fontstyle="italic", color=INK)
    ax.set_ylabel("gene families")
    style(ax)
    cls = "+".join(CLASSES)
    ax.set_title(f"Gain / loss labels: all accessory families ({cls}) vs. families in co-occurring pairs\n"
                 f"{rule}", loc="left", fontsize=15)
    ax.legend(frameon=False, loc="upper right", ncol=3)
    fig.text(0.01, 0.01, subs, fontsize=13, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.27)
    save(fig, args.out_dir, "fig1_gain_loss_by_study")

    # fig 2: by pangenome class, one panel per study
    fig, axes = plt.subplots(1, len(studies), figsize=(13.333, 7.5), sharey=False)
    axes = axes if len(studies) > 1 else [axes]
    for ax, (name, d) in zip(axes, studies):
        v = {dd: [int(d["bycls"][c][f"direction_{dd}"]) for c in CLASSES] for dd in DIRS}
        for c in CLASSES:
            for dd in DIRS:
                data_rows.append(["fig2", name, c, dd, d["bycls"][c][f"direction_{dd}"]])
        # grouped (not stacked) so small loss counts stay readable; value on each bar
        w = 0.27
        ymax = max(max(v[dd]) for dd in DIRS)
        for k, dd in enumerate(DIRS):
            x = [i + (k - 1) * (w + 0.02) for i in range(len(CLASSES))]
            ax.bar(x, v[dd], w, color=COLORS[dd], edgecolor=SURF, linewidth=2, label=dd)
            for xi, yi in zip(x, v[dd]):
                ax.text(xi, yi + ymax * 0.01, f"{yi:,}", ha="center", va="bottom",
                        fontsize=12, color=INK, rotation=90)
        ax.set_ylim(0, ymax * 1.18)
        ax.set_xticks(range(len(CLASSES)), [c.replace("_", "-") for c in CLASSES])
        ax.set_title(name, fontstyle="italic")
        ax.set_ylabel("gene families")
        style(ax)
    axes[-1].legend(frameon=False, loc="upper left")
    fig.suptitle(f"Gain / loss labels by ingroup pangenome class (core excluded)\n{rule}",
                 x=0.01, ha="left", fontsize=15)
    fig.text(0.01, 0.01, subs, fontsize=13, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.2, wspace=0.25)
    save(fig, args.out_dir, "fig2_gain_loss_by_class")

    # fig 3: size distribution (100% stacked; n on top)
    fig, axes = plt.subplots(1, len(studies), figsize=(13.333, 7.5), sharey=True)
    axes = axes if len(studies) > 1 else [axes]
    for ax, (name, d) in zip(axes, studies):
        agg = {b: {dd: 0 for dd in DIRS} for b in SIZE_ORDER}
        for r in d["size"]:
            if r["class"] in CLASSES:
                for dd in DIRS:
                    agg[r["ingroup_rep_strains_bin"]][dd] += int(r[dd])
        bins = [b for b in SIZE_ORDER if sum(agg[b].values())]
        tot = [sum(agg[b].values()) for b in bins]
        bottoms = [0.0] * len(bins)
        for dd in DIRS:
            frac = [100 * agg[b][dd] / t for b, t in zip(bins, tot)]
            ax.bar(range(len(bins)), frac, 0.8, bottom=bottoms, color=COLORS[dd], edgecolor=SURF,
                   linewidth=2, label=dd)
            bottoms = [x + y for x, y in zip(bottoms, frac)]
            for b in bins:
                data_rows.append(["fig3", name, b, dd, agg[b][dd]])
        for i, t in enumerate(tot):
            ax.text(i, 101, f"{t:,}", ha="center", va="bottom", fontsize=12, color=INK, rotation=90)
        ax.set_xticks(range(len(bins)), bins, rotation=45)
        ax.set_xlabel("ingroup representative genomes carrying the family")
        ax.set_title(name, fontstyle="italic", pad=48)
        ax.set_ylim(0, 100)
        style(ax)
    axes[0].set_ylabel("% of families in size bin")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, loc="lower right", ncol=3, bbox_to_anchor=(0.99, 0.01))
    fig.suptitle(f"Gain / loss share by family size ({cls}; numbers = families per bin)\n{rule}",
                 x=0.01, ha="left", fontsize=15)
    fig.text(0.01, 0.01, subs, fontsize=13, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.76, bottom=0.27, wspace=0.12)
    save(fig, args.out_dir, "fig3_size_distribution")

    with open(args.out_dir / "figure_data.tsv", "w") as fh:
        fh.write("figure\tstudy\tgroup\tdirection\tn_families\n")
        for r in data_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote figures to {args.out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
