#!/usr/bin/env python3
"""Slide figure (16:9) for one bin/pangenome_gain_loss_enrichment.py output dir.

Two panels: gain test set (left) and loss test set (right). Each panel shows the top
--top Pfam terms by p value (ties: term), as bars of -log10(q) (BH FDR). Bars with
q >= 0.05 are grey. Each bar is labelled "k/N  fold x": k = test-set families with the
term, N = annotated test-set families. A dashed line marks q = 0.05.
The title gives the run, genome counts (from the gain/loss PROVENANCE.yaml), the rule,
and the family scope. Counts and statistics only; no biological interpretation.

Usage:
  bin/pangenome_gain_loss_enrichment_figures.py --name "C. immitis (in) vs C. posadasii (out)" \\
      --enrichment-dir studies/fungi/<s>/analysis/gain_loss_enrichment/<run>/shell_soft_core \\
      --gain-loss-dir studies/fungi/<s>/analysis/gain_loss/<run> \\
      --out-dir studies/fungi/<s>/analysis/gain_loss_enrichment/figures --stem <run>_pfam_top
Writes <stem>.png/.svg/.pdf and <stem>.tsv (the numbers drawn).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import yaml  # noqa: E402

# Same palette and style as bin/pangenome_gain_loss_figures.py.
COLORS = {"gain": "#2a78d6", "loss": "#eb6834"}
NS_COLOR = "#b7b5ae"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SCOPE_TEXT = {"shell+soft_core": "shell + soft-core families",
              "soft_core+shell+cloud+singleton": "all accessory families (soft-core, shell, cloud, singleton)",
              "cloud+singleton": "cloud + singleton families"}

plt.rcParams.update({
    "font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15, "xtick.labelsize": 14,
    "ytick.labelsize": 14, "legend.fontsize": 14, "figure.titlesize": 17,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "svg.hashsalt": "gain_loss_enrich",
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def read_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def top_terms(rows: list[dict], test_set: str, vocab: str, n: int) -> list[dict]:
    sel = [r for r in rows if r["test_set"] == test_set and r["vocab"] == vocab]
    sel.sort(key=lambda r: (float(r["p"]), r["term"]))
    return sel[:n]


def short_label(r: dict, width: int = 34) -> str:
    name = r["name"] or r["term"]
    pid, _, de = name.partition(": ")
    text = f"{pid} ({r['term']})" if de else name
    return text if len(text) <= width else text[: width - 1] + "…"


def neglog(q: float) -> float:
    return -math.log10(max(q, 1e-300))


def panel(ax, rows: list[dict], test_set: str, n_set: int, n_all: int, xmax: float):
    ys = list(range(len(rows)))[::-1]
    for y, r in zip(ys, rows):
        q = float(r["q"])
        color = COLORS[test_set] if q < 0.05 else NS_COLOR
        v = neglog(q)
        ax.barh(y, v, height=0.72, color=color, edgecolor=SURF, linewidth=1.5)
        fold = r["fold"]
        ftxt = "NA" if fold == "NA" else f"{float(fold):.1f}x"
        ax.text(v + xmax * 0.015, y, f"{r['k']}/{r['n_set']}  {ftxt}  q={q:.2g}", va="center", ha="left",
                color=INK, fontsize=13, bbox=dict(facecolor=SURF, edgecolor="none", pad=1.0))
    ax.set_yticks(ys)
    ax.set_yticklabels([short_label(r) for r in rows])
    ax.axvline(neglog(0.05), color=INK2, linestyle="--", linewidth=1.2)
    ax.set_xlim(0, xmax * 1.6)
    ax.set_xlabel("-log10 q (BH)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    n_sig = sum(1 for r in rows if float(r["q"]) < 0.05)
    ax.set_title(f"{test_set}: {n_set:,} of {n_all:,} families have >= 1 Pfam\n"
                 f"top {len(rows)} by p; {n_sig} with q < 0.05", loc="left", fontsize=15)
    if not rows:
        ax.text(0.5, 0.5, "no term tested", transform=ax.transAxes, ha="center")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="run name for the title")
    ap.add_argument("--enrichment-dir", required=True, type=Path)
    ap.add_argument("--gain-loss-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--stem", required=True)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--vocab", default="pfam")
    args = ap.parse_args(argv)

    rows = read_rows(args.enrichment_dir / "enrichment_domains.tsv.gz")
    gprov = yaml.safe_load(open(args.gain_loss_dir / "PROVENANCE.yaml"))
    eprov = yaml.safe_load(open(args.enrichment_dir / "PROVENANCE.yaml"))
    with open(args.enrichment_dir / "coverage.tsv") as fh:
        cov = {r["direction"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(args.enrichment_dir / "tests.tsv") as fh:
        tests = {(r["test_set"], r["vocab"]): r for r in csv.DictReader(fh, delimiter="\t")}
    n = gprov["strain_counts"]
    par = gprov["parameters"]
    scope = "+".join(eprov["parameters"]["classes"])
    test_sets = eprov["parameters"]["test_sets"]

    tops = {ts: top_terms(rows, ts, args.vocab, args.top) for ts in test_sets}
    vals = [neglog(float(r["q"])) for t in tops.values() for r in t] + [neglog(0.05)]
    xmax = max(vals)
    fig, axes = plt.subplots(1, len(test_sets), figsize=(16, 9))
    if len(test_sets) == 1:
        axes = [axes]
    for ax, ts in zip(axes, test_sets):
        t = tests[(ts, args.vocab)]
        panel(ax, tops[ts], ts, int(t["N_set_annotated"]), int(cov[ts]["n_families"]) if ts in cov else 0,
              xmax)
    fig.suptitle(
        f"{args.name}: {n['in_rep'] + n['in_nonrep']} ingroup genomes "
        f"({n['in_rep']} representatives), outgroup {n['out_all']} genomes\n"
        f"Pfam enrichment of {SCOPE_TEXT.get(scope, scope)}; gain = outgroup frequency <= "
        f"{par['gain_max_frac']:g}, loss = >= {par['loss_min_frac']:g}\n"
        f"Run {eprov['label']}; universe: {scope.replace('+', ' + ')} families with >= 1 Pfam "
        f"(M = {int(tests[(test_sets[0], args.vocab)]['M_universe_annotated']):,})",
        x=0.01, ha="left", fontsize=16)
    fig.text(0.01, 0.01, "Bar label: k/N = test-set families with the domain / annotated test-set families; "
             "fold = k / expected.\nGrey bar: q >= 0.05. Dashed line: q = 0.05. "
             "Annotation = family representative protein only (Pfam-A 38.2, hmmsearch --cut_ga).", fontsize=13, color=INK2)
    fig.subplots_adjust(top=0.79, bottom=0.17, left=0.185, right=0.985)
    fig.subplots_adjust(wspace=0.85)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        kw = {"dpi": 200} if ext == "png" else {}
        md = {"Date": None} if ext == "svg" else ({"CreationDate": None} if ext == "pdf" else {})
        fig.savefig(args.out_dir / f"{args.stem}.{ext}", facecolor=SURF, metadata=md, **kw)
    plt.close(fig)
    with open(args.out_dir / f"{args.stem}.tsv", "w") as fh:
        cols = ["test_set", "vocab", "term", "name", "k", "n_set", "K", "M", "fold", "p", "q"]
        fh.write("\t".join(["label", "scope"] + cols) + "\n")
        for ts in test_sets:
            for r in tops[ts]:
                fh.write("\t".join([r["label"], r["scope"]] + [r[c] for c in cols]) + "\n")
    print(f"wrote {args.out_dir / args.stem}.png/.svg/.pdf/.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
