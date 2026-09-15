#!/usr/bin/env python3
"""Standard Roary/Panaroo-style static pangenome summary figures --
design spec component 9's precursor: per the 2026-09-14 conversation
("do the common/existing-tool visualization first"), these are the fast,
low-risk figures that don't require the full interactive canvas-based
report, generated directly from data this study already has for real.

Produces, from `frequency_table.tsv` (+ optionally `presence_matrix.tsv`
and `pair_classification.tsv`):
  1. frequency_distribution.png  -- family-frequency histogram, colored by
     core/soft-core/shell/cloud/singleton band (the classic Roary
     `plot_roary`-style frequency plot).
  2. core_shell_cloud_pie.png    -- pie chart of family counts per band.
  3. presence_absence_matrix.png -- families (rows, frequency-sorted) x
     strains (columns) raster, the classic Roary/Phandango-style sorted
     presence/absence matrix. Rasterized via matplotlib imshow, not one
     drawn element per cell, so it stays fast even at this study's real
     scale (tens of thousands of families).
  4. accumulation_curve.png      -- pangenome/core-genome size vs. number
     of strains sampled, averaged over random strain-order permutations
     (the classic Roary/Panaroo rarefaction-curve figure).
  5. pair_classification_summary.png -- bar chart of component 8's six
     classification labels (only if --pair_classification is given).

Usage:
  plot_pangenome_summary.py --frequency_table frequency_table.tsv \\
      --matrix presence_matrix.tsv \\
      --pair_classification pair_classification.tsv \\
      --out_dir figures/
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402

BAND_ORDER = ["core", "soft_core", "shell", "cloud", "singleton"]
BAND_COLORS = {
    "core": "#2c7fb8",
    "soft_core": "#7fcdbb",
    "shell": "#fed976",
    "cloud": "#fd8d3c",
    "singleton": "#bdbdbd",
}


def read_frequency_table(path: str) -> list[dict]:
    with open_maybe_compressed(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def band_counts(frequency_table: list[dict]) -> dict[str, int]:
    counts = {b: 0 for b in BAND_ORDER}
    for row in frequency_table:
        counts[row["bin"]] = counts.get(row["bin"], 0) + 1
    return counts


def matrix_to_binary_array(
    matrix: PresenceMatrix, family_order: list[str], strain_order: list[str]
) -> np.ndarray:
    """Boolean [n_families x n_strains] array in the given row/column order."""
    arr = np.zeros((len(family_order), len(strain_order)), dtype=bool)
    for i, fam in enumerate(family_order):
        for j, strain in enumerate(strain_order):
            arr[i, j] = matrix.is_present(fam, strain)
    return arr


def accumulation_curve(
    arr: np.ndarray, n_permutations: int = 20, seed: int = 0
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """[families x strains] boolean array -> (pan_mean, pan_std, core_mean,
    core_std), each length n_strains, averaged over `n_permutations` random
    strain orderings. `pan[k]` = number of families present in at least one
    of the first k+1 (randomly ordered) strains; `core[k]` = number present
    in ALL of the first k+1."""
    n_families, n_strains = arr.shape
    rng = random.Random(seed)
    pan_runs = np.zeros((n_permutations, n_strains), dtype=int)
    core_runs = np.zeros((n_permutations, n_strains), dtype=int)
    strain_indices = list(range(n_strains))
    for p in range(n_permutations):
        order = strain_indices[:]
        rng.shuffle(order)
        cum_any = np.zeros(n_families, dtype=bool)
        cum_all = np.ones(n_families, dtype=bool)
        for k, s in enumerate(order):
            cum_any |= arr[:, s]
            cum_all &= arr[:, s]
            pan_runs[p, k] = cum_any.sum()
            core_runs[p, k] = cum_all.sum()
    return (
        pan_runs.mean(axis=0), pan_runs.std(axis=0),
        core_runs.mean(axis=0), core_runs.std(axis=0),
    )


def plot_frequency_distribution(frequency_table: list[dict], out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for band in BAND_ORDER:
        freqs = [float(r["frequency"]) for r in frequency_table if r["bin"] == band]
        if freqs:
            ax.hist(freqs, bins=40, range=(0, 1), color=BAND_COLORS[band],
                     label=f"{band} (n={len(freqs)})", alpha=0.85)
    ax.set_xlabel("Fraction of ingroup strains carrying the family")
    ax.set_ylabel("Number of gene families")
    ax.set_title("Afumigatus_pangenome: family-frequency distribution")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_core_shell_cloud_pie(counts: dict[str, int], out_path: str) -> None:
    labels = [b for b in BAND_ORDER if counts.get(b, 0) > 0]
    sizes = [counts[b] for b in labels]
    colors = [BAND_COLORS[b] for b in labels]
    fig, ax = plt.subplots(figsize=(6, 6))
    total = sum(sizes)
    ax.pie(
        sizes, labels=[f"{b}\n({n}, {n/total:.1%})" for b, n in zip(labels, sizes)],
        colors=colors, startangle=90,
    )
    ax.set_title(f"Afumigatus_pangenome: {total} gene families")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_presence_absence_matrix(arr: np.ndarray, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 10))
    ax.imshow(arr, aspect="auto", cmap="Greys", interpolation="nearest")
    ax.set_xlabel(f"Strains (n={arr.shape[1]})")
    ax.set_ylabel(f"Gene families (n={arr.shape[0]}), sorted by frequency")
    ax.set_title("Afumigatus_pangenome: presence/absence matrix")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_accumulation_curve(
    pan_mean, pan_std, core_mean, core_std, out_path: str
) -> None:
    n = len(pan_mean)
    x = np.arange(1, n + 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, pan_mean, color="#2c7fb8", label="Pangenome (union)")
    ax.fill_between(x, pan_mean - pan_std, pan_mean + pan_std, color="#2c7fb8", alpha=0.2)
    ax.plot(x, core_mean, color="#d95f02", label="Core (intersection)")
    ax.fill_between(x, core_mean - core_std, core_mean + core_std, color="#d95f02", alpha=0.2)
    ax.set_xlabel("Number of strains sampled")
    ax.set_ylabel("Number of gene families")
    ax.set_title("Afumigatus_pangenome: accumulation curve")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_pair_classification_summary(counts: dict[str, int], out_path: str) -> None:
    order = sorted(counts, key=lambda k: -counts[k])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(order, [counts[k] for k in order], color="#4d4d4d")
    ax.set_ylabel("Number of family pairs")
    ax.set_title("Afumigatus_pangenome: co-occurring pair classification (component 8)")
    ax.tick_params(axis="x", rotation=30)
    for i, k in enumerate(order):
        ax.text(i, counts[k], f"{counts[k]:,}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--matrix")
    ap.add_argument("--pair_classification")
    ap.add_argument("--n_permutations", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading frequency table...", file=sys.stderr)
    frequency_table = read_frequency_table(args.frequency_table)
    counts = band_counts(frequency_table)
    print(f"Band counts: {counts}", file=sys.stderr)

    plot_frequency_distribution(frequency_table, str(out_dir / "frequency_distribution.png"))
    plot_core_shell_cloud_pie(counts, str(out_dir / "core_shell_cloud_pie.png"))
    print("Wrote frequency_distribution.png, core_shell_cloud_pie.png", file=sys.stderr)

    if args.matrix:
        print("Loading presence matrix...", file=sys.stderr)
        matrix = PresenceMatrix.from_tsv(args.matrix)
        freq_by_family = {r["family"]: float(r["frequency"]) for r in frequency_table}
        family_order = sorted(
            (f for f in matrix.families if f in freq_by_family),
            key=lambda f: -freq_by_family[f],
        )
        strain_order = list(matrix.strains)
        print(f"Building {len(family_order)}x{len(strain_order)} binary array...", file=sys.stderr)
        arr = matrix_to_binary_array(matrix, family_order, strain_order)
        plot_presence_absence_matrix(arr, str(out_dir / "presence_absence_matrix.png"))
        print("Wrote presence_absence_matrix.png", file=sys.stderr)

        print(f"Computing accumulation curve ({args.n_permutations} permutations)...", file=sys.stderr)
        pan_mean, pan_std, core_mean, core_std = accumulation_curve(
            arr, n_permutations=args.n_permutations, seed=args.seed,
        )
        plot_accumulation_curve(
            pan_mean, pan_std, core_mean, core_std, str(out_dir / "accumulation_curve.png"),
        )
        print("Wrote accumulation_curve.png", file=sys.stderr)

    if args.pair_classification:
        print("Reading pair classification...", file=sys.stderr)
        class_counts: dict[str, int] = {}
        with open_maybe_compressed(args.pair_classification) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = header.index("classification")
            for line in fh:
                label = line.rstrip("\n").split("\t")[idx]
                class_counts[label] = class_counts.get(label, 0) + 1
        plot_pair_classification_summary(class_counts, str(out_dir / "pair_classification_summary.png"))
        print("Wrote pair_classification_summary.png", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
