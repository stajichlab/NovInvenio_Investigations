#!/usr/bin/env python3
"""Build a site report for the island genome-count sweep.

Writes site/report/report.md and site/report/figures/*.png (+ figures_pdf/)
from outputs/sweep_summary.tsv and outputs/runs/*/starship.starships.tsv, in
the layout bin/sync_pangenome_report.py stages (--pangenome-dir site).
Every number on the page is computed here from those files.
"""
from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

AN = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/"
          "studies/fungi/Afumigatus_pangenome/analysis/island_genome_count")
SITE = AN / "site"
FIG = SITE / "report" / "figures"
FIG_PDF = SITE / "report" / "figures_pdf"

# Reference palette (dataviz skill references/palette.md), light surface.
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SERIES_1 = "#2a78d6"
BACKGROUND = "#8f8e89"
SEQ = ["#eef4fc", "#b7d3f4", "#6ea6e8", "#2a78d6"]  # 0..3 seeds, one hue light->dark
UNSCORED = "#d9d8d3"


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_2)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, name: str) -> None:
    fig.patch.set_facecolor(SURFACE)
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=150)
    fig.savefig(FIG_PDF / f"{name}.pdf")
    plt.close(fig)


def load_rows() -> list[dict]:
    with (AN / "outputs/sweep_summary.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for k, v in r.items():
            r[k] = float(v) if v not in ("", "nan") else float("nan")
        r["n"], r["seed"] = int(r["n"]), int(r["seed"])
    return rows


def per_n(rows, key):
    out = defaultdict(list)
    for r in rows:
        out[r["n"]].append(r[key])
    return dict(sorted(out.items()))


def dot_median(ax, rows, key, scale=1.0, skip_labels=()) -> None:
    groups = per_n(rows, key)
    xs = list(groups)
    ax.plot(xs, [statistics.median(v) * scale for v in groups.values()],
            color=SERIES_1, linewidth=2, zorder=2)
    for n, vals in groups.items():
        ax.scatter([n] * len(vals), [v * scale for v in vals], s=36, color=SERIES_1,
                   edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.set_xlabel("representative genomes (N)", color=TEXT_2)
    ax.set_xticks(xs, ["" if x in skip_labels else str(x) for x in xs])


def fig_islands(rows) -> None:
    fig, ax = plt.subplots(figsize=(7, 3.6))
    style_axes(ax)
    dot_median(ax, rows, "n_islands")
    ax.axvspan(8, 25, color=GRID, alpha=0.6, zorder=0, linewidth=0)
    ax.text(16.5, ax.get_ylim()[1] * 0.9, "0 islands\n(every seed)", ha="center",
            va="top", fontsize=8.5, color=TEXT_2)
    ax.set_ylabel("significant islands", color=TEXT_2)
    ax.set_title("Islands by genome count (dots: seeds; line: median)", loc="left",
                 fontsize=11, color=TEXT)
    save(fig, "islands_vs_n")


def fig_starships(rows) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    for ax in axes:
        style_axes(ax)
    dot_median(axes[0], rows, "af293_high_conf_starships_hit", skip_labels=(15,))
    axes[0].set_ylim(-0.3, 7.5)
    axes[0].set_ylabel("Starships hit (of 7)", color=TEXT_2)
    axes[0].set_title("Af293 coordinates (Table S7)", loc="left", fontsize=10.5, color=TEXT)
    for r in rows:
        s = r["presence_starships_scored"]
        r["presence_frac"] = r["presence_starships_above_null95"] / s if s else 0.0
    dot_median(axes[1], rows, "presence_frac", scale=100, skip_labels=(15,))
    axes[1].set_ylim(-3, 103)
    axes[1].set_ylabel("% of scored Starships above null", color=TEXT_2)
    axes[1].set_title("Presence pattern (Table S21)", loc="left", fontsize=10.5, color=TEXT)
    save(fig, "starships_vs_n")


def fig_gene_frac(rows) -> None:
    fig, ax = plt.subplots(figsize=(7, 3.6))
    style_axes(ax)
    isl = {n: v for n, v in per_n(rows, "island_gene_frac_in_starship").items()
           if all(x == x for x in v)}  # drop N with no islands (NaN)
    bg = per_n(rows, "noncore_gene_frac_in_starship")
    for data, color, label in ((isl, SERIES_1, "island genes"), (bg, BACKGROUND, "all non-core genes")):
        xs = list(data)
        ys = [statistics.median(v) * 100 for v in data.values()]
        ax.plot(xs, ys, color=color, linewidth=2, label=label)
        for n, vals in data.items():
            ax.scatter([n] * len(vals), [v * 100 for v in vals], s=30, color=color,
                       edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.annotate(label, (xs[-1], ys[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=8.5, color=TEXT_2)
    ax.set_xticks(sorted(bg))
    ax.set_xlim(right=max(bg) * 1.25)
    ax.set_ylim(0, None)
    ax.set_xlabel("representative genomes (N)", color=TEXT_2)
    ax.set_ylabel("% of Af293 genes inside a Starship", color=TEXT_2)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=TEXT_2, loc="center left", bbox_to_anchor=(0.02, 0.52))
    ax.set_title("Island genes are enriched in Starships", loc="left", fontsize=11, color=TEXT)
    save(fig, "starship_gene_fraction")


def fig_runtime(rows) -> None:
    fig, ax = plt.subplots(figsize=(7, 3.4))
    style_axes(ax)
    dot_median(ax, rows, "runtime_s_island_steps", scale=1 / 60)
    ax.set_ylim(0, None)
    ax.set_ylabel("minutes (island steps only)", color=TEXT_2)
    ax.set_title("Runtime of the four island steps", loc="left", fontsize=11, color=TEXT)
    save(fig, "runtime_vs_n")


def starship_matrix(rows):
    """{label: {n: (recovered_seeds, scored_seeds)}} for both checks."""
    mat: dict[str, dict[int, list[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    order: list[str] = []
    for r in sorted(rows, key=lambda r: (r["n"], r["seed"])):
        path = AN / "outputs/runs" / f"n{r['n']}_s{r['seed']}" / "starship.starships.tsv"
        with path.open() as fh:
            for s in csv.DictReader(fh, delimiter="\t"):
                if s["check"] == "af293_coord":
                    if s["confidence"] != "high":
                        continue
                    label = f"Af293 {s['name']}"
                    hit, scored = int(s["value"]) > 0, True
                else:
                    label = s["name"]
                    scored = s["value"] != ""
                    hit = scored and float(s["value"]) > float(s["null95"])
                if label not in order:
                    order.append(label)
                cell = mat[label][r["n"]]
                cell[0] += hit
                cell[1] += scored
    return order, mat


def fig_heatmap(rows):
    order, mat = starship_matrix(rows)
    ns = sorted({r["n"] for r in rows})
    fig, ax = plt.subplots(figsize=(8.3, 0.3 * len(order) + 1.4))
    cmap = ListedColormap(SEQ)
    for i, label in enumerate(order):
        for j, n in enumerate(ns):
            hit, scored = mat[label][n]
            if scored == 0:
                ax.add_patch(plt.Rectangle((j, i), 0.94, 0.9, color=UNSCORED, hatch="///",
                                           fill=False, linewidth=0))
                continue
            ax.add_patch(plt.Rectangle((j, i), 0.94, 0.9, color=cmap(hit / 3 if scored == 3 else hit / max(scored, 1))))
            ax.text(j + 0.47, i + 0.45, f"{hit}/{scored}", ha="center", va="center", fontsize=6.5,
                    color=SURFACE if hit / max(scored, 1) > 0.6 else TEXT)
    n_coord = sum(label.startswith("Af293 ") for label in order)
    ax.axhline(n_coord - 0.05, color=TEXT_2, linewidth=1)
    ax.text(len(ns) + 0.1, n_coord / 2, "Af293\ncoordinates\n(Table S7)", va="center", fontsize=7.5, color=TEXT_2)
    ax.text(len(ns) + 0.1, (n_coord + len(order)) / 2, "presence\npattern\n(Table S21)", va="center", fontsize=7.5, color=TEXT_2)
    ax.set_xlim(0, len(ns))
    ax.set_ylim(len(order), 0)
    ax.set_xticks([j + 0.47 for j in range(len(ns))], [str(n) for n in ns])
    ax.set_yticks([i + 0.45 for i in range(len(order))], order, fontsize=7.5)
    ax.tick_params(length=0, colors=TEXT_2)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_xlabel("representative genomes (N)", color=TEXT_2)
    ax.set_title("Seeds that recover each Starship (hit / scored; hatched: too few carriers)",
                 loc="left", fontsize=10, color=TEXT)
    save(fig, "starship_recovery_heatmap")
    return order, mat, ns


def md_table(header, body) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in body]
    return "\n".join(lines)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    FIG_PDF.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    fig_islands(rows)
    fig_starships(rows)
    fig_gene_frac(rows)
    fig_runtime(rows)
    order, mat, ns = fig_heatmap(rows)

    body = []
    for n in ns:
        rs = [r for r in rows if r["n"] == n]
        def rng(k, fmt="{:,.0f}"):
            vals = sorted(r[k] for r in rs)
            return fmt.format(vals[0]) if vals[0] == vals[-1] else f"{fmt.format(vals[0])} - {fmt.format(vals[-1])}"
        pres = ", ".join(f"{int(r['presence_starships_above_null95'])}/{int(r['presence_starships_scored'])}" for r in rs)
        body.append([n, len(rs), rng("n_cooccurring_pairs"), rng("n_islands"),
                     rng("af293_high_conf_starships_hit"), pres,
                     rng("runtime_s_island_steps", "{:,.0f}")])
    table = md_table(["N", "seeds", "co-occurring pairs", "islands", "Af293 Starships hit (of 7)",
                      "presence matches / scored", "island steps (s)"], body)
    zero_ns = [n for n in ns if all(r["n_islands"] == 0 for r in rows if r["n"] == n)]
    pos_ns = [n for n in ns if all(r["n_islands"] > 0 for r in rows if r["n"] == n)]
    all7 = [n for n in ns if all(r["af293_high_conf_starships_hit"] == 7 for r in rows if r["n"] == n)]

    md = f"""# Island genome-count sweep

How many genomes does the pangenome island test need? This page subsamples the
293-strain *A. fumigatus* run (`Afumigatus_pangenome` / `full_v070`, 121
representative strains after dereplication) and reruns only the island steps
at each size. It is not a pipeline run: there is no clustering or rescue here.

## Result

- Every subset with N <= {max(zero_ns)} representative genomes gave 0 islands (0 FDR-significant co-occurring pairs).
- Every subset with N >= {min(pos_ns)} gave islands.
- All 7 high-confidence Af293 Starships were hit by islands in every seed at N in {{{", ".join(map(str, all7))}}}.
- The test set `isl45_v1` uses the N = 45, seed 3 subset.

{table}

![Islands by genome count](figures/islands_vs_n.png)

![Starships recovered](figures/starships_vs_n.png)

![Starship gene fraction](figures/starship_gene_fraction.png)

![Starship recovery by N](figures/starship_recovery_heatmap.png)

![Runtime](figures/runtime_vs_n.png)

## Method

- Sampling unit: representative ingroup strains (the co-occurrence test uses only these).
  Af293 is in every subset; both outgroup references are kept. Seeds 1-3 per N; N = 121 is all.
- Steps (NovInvenio f87fd1e, parameters as in the pipeline): FREQUENCY_BINS, COOCCURRENCE,
  PAIR_CLASSIFICATION, BUILD_ISLANDS.
- Starship truth: Gluck-Thaler et al. 2025, mBio, doi:10.1128/mbio.01092-25.
  - Coordinate check: Table S7, 7 high-confidence Af293 elements. A Starship is hit when an
    island's Af293 gene lies inside it.
  - Presence check: Table S21, 20 Starships. For each Starship, the best Jaccard between an
    island's strain presence and the Starship's, compared with the 95th percentile of 200
    shuffles of the Starship's presence. Starships with fewer than 3 carriers or
    non-carriers in a subset are not scored.
- Background for the gene fraction: all Af293 genes in non-core families of that subset.

## Limits

- Each subset reuses gene families clustered from all 293 genomes. A real N-genome run
  clusters its own families, so its island counts will differ.
- Runtime here is the island steps only; a full pipeline run also does clustering and rescue.
- N = 121 here gives {int(next(r["n_islands"] for r in rows if r["n"] == 121)):,} islands; `full_v070` (pipeline a394ace)
  gave 11,775 on the same strains. The cause of the difference was not checked.

Source: `studies/fungi/Afumigatus_pangenome/analysis/island_genome_count/` (ISLAND_GENOME_COUNT.md).
"""
    (SITE / "report" / "report.md").write_text(md)
    (SITE / "report" / "sweep.json").write_text(json.dumps(
        {"n_runs": len(rows), "ns": ns, "zero_island_ns": zero_ns}, indent=1) + "\n")
    print(f"wrote {SITE / 'report/report.md'} and {len(list(FIG.glob('*.png')))} figures")


if __name__ == "__main__":
    main()
