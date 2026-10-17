#!/usr/bin/env python3
"""Slide figures (16:9) from one or more bin/pangenome_gain_loss.py output dirs.

Usage:
  bin/pangenome_gain_loss_figures.py --out-dir <dir> \\
      --study "A. fumigatus=studies/fungi/Afumigatus_pangenome/analysis/gain_loss/full_v072" \\
      --study "C. immitis=studies/fungi/coccidioides_pangenome/analysis/gain_loss/immitis_in_posadasii_out_v5"

Writes fig1_gain_loss_by_study, fig2_gain_loss_by_class, fig3_size_distribution,
fig4_gain_loss_shell_softcore (.png, .svg, .pdf), figure_data.tsv (the numbers
drawn in fig1-3) and figure_data_shell_softcore.tsv (the numbers drawn in fig4).
Counts only; no biological interpretation.

fig4 keeps the persistent accessory classes (shell + soft-core) apart from the
rare classes (cloud + singleton), counts and 100%-stacked shares.
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
# fig4 class groups (order = panel order). Never mixed in one bar.
CLASS_GROUPS = (("shell+soft_core", ("shell", "soft_core")),
                ("cloud+singleton", ("cloud", "singleton")))
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


# Ingroup frequency-bin cutoffs (fraction of ingroup representative genomes) of
# nf_NovInvenio bin/pangenome_frequency_bins.py. Not stored in the gain_loss
# outputs; the default values below match the full_v072 and
# immitis_in_posadasii_out_v5 Nextflow run reports. Used for title text only.
BIN_CUTOFFS = {"core": 0.95, "soft_core": 0.90, "shell": 0.15}
GROUP_TEXT = {
    "shell+soft_core": (f"shell + soft-core families\nin {BIN_CUTOFFS['shell']:.0%} to "
                        f"<{BIN_CUTOFFS['core']:.0%} of ingroup rep. genomes"),
    "cloud+singleton": (f"cloud + singleton families\nin <{BIN_CUTOFFS['shell']:.0%} of "
                        f"ingroup rep. genomes (singleton: 1 genome)"),
}


def class_group_counts(bycls: dict, groups=CLASS_GROUPS) -> dict:
    """bycls: class -> summary_direction_by_class.tsv row. Returns
    {group: {direction: int}} summed over the group's classes (gain/loss/ambiguous)."""
    out = {}
    for g, members in groups:
        out[g] = {dd: sum(int(bycls[c][f"direction_{dd}"]) for c in members) for dd in DIRS}
    return out


def class_group_rows(studies, groups=CLASS_GROUPS) -> list[list]:
    """Tidy rows [study, class_group, direction, count, fraction] (fraction of the
    study's class-group total; 0 when the group is empty)."""
    rows = []
    for name, d in studies:
        for g, v in class_group_counts(d["bycls"], groups).items():
            tot = sum(v.values())
            for dd in DIRS:
                rows.append([name, g, dd, v[dd], f"{(v[dd] / tot if tot else 0.0):.6f}"])
    return rows


def _spread(ys: list[float], gap: float) -> list[float]:
    """Push label y positions upward so neighbours are at least `gap` apart."""
    out = []
    for y in ys:
        out.append(max(y, out[-1] + gap) if out else y)
    return out


def _stack_panel(ax, names, vals: list[dict], pct: bool, width=0.5, min_frac=0.10):
    """One bar per study, stacked gain/loss/ambiguous. Segment labels inside when the
    segment is >= min_frac of the axis height, else beside the bar: a square in the
    series colour plus the value in text ink, spread apart vertically."""
    tots = [sum(v.values()) for v in vals]
    ymax = 100.0 if pct else max(tots + [1])
    for i, (v, tot) in enumerate(zip(vals, tots)):
        bottom, side = 0.0, []
        for dd in DIRS:
            h = (100.0 * v[dd] / tot if tot else 0.0) if pct else v[dd]
            ax.bar(i, h, width, bottom=bottom, color=COLORS[dd], edgecolor=SURF, linewidth=2,
                   label=dd if i == 0 else None)
            txt = f"{100.0 * v[dd] / tot:.1f}%" if pct else f"{v[dd]:,}"
            if v[dd] and h / ymax >= min_frac:
                ax.text(i, bottom + h / 2, txt, ha="center", va="center", color="white",
                        fontsize=14, fontweight="bold")
            elif v[dd]:
                side.append((bottom + h / 2, dd, txt))
            bottom += h
        ys = _spread([y for y, _, _ in side], 0.12 * ymax)
        if ys and ys[-1] > 0.97 * ymax:   # keep side labels inside the axes: shift down
            ys = [y - (ys[-1] - 0.97 * ymax) for y in ys]
        for y, (_, dd, t) in zip(ys, side):
            ax.plot([i + width / 2 + 0.08], [y], marker="s", markersize=10, color=COLORS[dd],
                    linestyle="none", clip_on=False)
            ax.text(i + width / 2 + 0.15, y, t, ha="left", va="center", color=INK, fontsize=14)
    ax.set_xlim(-0.45, len(names) - 0.5 + 0.35)
    ax.set_ylim(0, 100 if pct else ymax * 1.12)
    ax.set_xticks(range(len(names)), [f"{n}\nn={t:,}" for n, t in zip(names, tots)])
    for lab in ax.get_xticklabels():
        lab.set_fontstyle("italic")
    style(ax)


def limitation_text(studies) -> str:
    """One line per small outgroup size (< 10 genomes: the 0.1 / 0.9 rule cannot
    resolve finer than 1/n)."""
    by_n: dict[int, list[str]] = {}
    for name, d in studies:
        n = d["prov"]["strain_counts"]["out_all"]
        if n < 10:
            by_n.setdefault(n, []).append(name)
    lines = []
    for n, names in sorted(by_n.items()):
        vals = ", ".join(f"{k}/{n}" for k in range(n + 1))
        lines.append(f"Limitation: {' and '.join(names)}: outgroup = {n} genomes, so outgroup "
                     f"frequency is only {vals}.\n    gain = in no outgroup genome, "
                     f"loss = in all {n}, ambiguous = in some but not all.")
    return "\n".join(lines)


def fig4(studies, out_dir: Path, rule: str, subs: str):
    names = [n for n, _ in studies]
    counts = [class_group_counts(d["bycls"]) for _, d in studies]
    fig, axes = plt.subplots(2, 2, figsize=(13.333, 7.5), sharex="col")
    for col, (g, _) in enumerate(CLASS_GROUPS):
        bw = 0.5 if len(names) <= 2 else 0.56   # wider bars when 3 studies share a panel
        vals = [c[g] for c in counts]
        _stack_panel(axes[0][col], names, vals, pct=False, width=bw)
        _stack_panel(axes[1][col], names, vals, pct=True, width=bw)
        axes[0][col].set_title(GROUP_TEXT[g], loc="left", fontsize=14)
        axes[0][col].tick_params(labelbottom=False)
    axes[0][0].set_ylabel("gene families")
    axes[1][0].set_ylabel("% of bar")
    h, lab = axes[0][0].get_legend_handles_labels()
    fig.suptitle("Gain / loss labels: shell + soft-core families shown apart from cloud + "
                 f"singleton families\n{rule}", x=0.01, y=0.99, ha="left", va="top", fontsize=15)
    lim = limitation_text(studies)
    foot = subs + ("\n" + lim if lim else "")
    fig.text(0.01, 0.01, foot, fontsize=13, color=INK2, va="bottom")
    nfoot = foot.count("\n") + 1
    bottom = 0.105 + 0.03 * nfoot
    fig.subplots_adjust(left=0.08, right=0.98, top=0.80, bottom=bottom, hspace=0.1, wspace=0.22)
    fig.legend(h, lab, frameon=False, loc="lower right", ncol=3, bbox_to_anchor=(0.99, 0.87),
               fontsize=14, handlelength=1.2, columnspacing=1.2)
    save(fig, out_dir, "fig4_gain_loss_shell_softcore")
    with open(out_dir / "figure_data_shell_softcore.tsv", "w") as fh:
        fh.write("study\tclass_group\tdirection\tcount\tfraction\n")
        for r in class_group_rows(studies):
            fh.write("\t".join(str(x) for x in r) + "\n")


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
    fig4(studies, args.out_dir, rule, subs)
    print(f"wrote figures to {args.out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
