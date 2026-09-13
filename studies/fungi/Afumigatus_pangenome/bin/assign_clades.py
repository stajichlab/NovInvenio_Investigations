#!/usr/bin/env python3
"""Mash-distance-based clade assignment -- a cheap, alignment-free stand-in
for the paper's DAPC clades, to fill config.csv's currently-empty
`TaxonGroup` column (see PANGENOME_CLUSTER_PROFILE_NOTES.md's "Open items").

Method: whole-genome Mash sketch/distance (same tool `dereplicate_strains.py`
already uses for dedup) -> classical PCoA (metric MDS on the distance matrix,
via double-centering + eigendecomposition -- the "simplified PCA" analog of
DAPC's PCA-then-discriminant-analysis step, skipping the discriminant step
since we have no prior group labels to discriminate) -> k-means on the PCoA
coordinates, with k chosen by silhouette score over a small candidate range.

This is NOT a substitute for a real SNP/allele-based DAPC (see
`parsnp_snps.py` for that path) -- Mash distance is a whole-genome k-mer
sketch distance, coarser than an actual variant-based genetic distance, so
treat its clusters as a working stratification variable to unblock
cooccurrence.py's clade-permutation null, not as a validated population
structure result.

Ingroup-only by default (`--groups IN`): running this across IN+OUT on the
real 295-strain config found k=2 "best" by silhouette, but the split was
294 vs. 1 -- the single *A. lentulus* outgroup strain, whose Mash distance
to every *A. fumigatus* strain (different species) dwarfs any within-species
distance and swamps the k-means/silhouette selection entirely. That is a
real, reproducible artifact of mixing an outgroup species into a same-scale
distance clustering, not a bug -- exclude outgroup strains before assigning
within-species clades.

Usage:
  assign_clades.py --config config.csv --data_dir data_dir \
      --k_range 2,10 --output clade_assignments.tsv
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import squareform, pdist


def parse_mash_dist_matrix(lines: list[str]) -> tuple[list[str], np.ndarray]:
    """Parse `mash dist -t` output (full query x reference distance matrix)
    into (names, symmetric distance matrix). Unlike
    dereplicate_strains.parse_mash_dist (which thresholds into dedup
    groups), this keeps the full continuous distances for PCoA."""
    if not lines:
        return [], np.zeros((0, 0))
    header = lines[0].lstrip("#").split("\t")
    names = header[1:]
    n = len(names)
    mat = np.zeros((n, n))
    for i, row in enumerate(lines[1:]):
        parts = row.split("\t")
        for j, dist_str in enumerate(parts[1:]):
            mat[i, j] = float(dist_str)
    # Symmetrize (mash's own rounding can leave mat[i,j] != mat[j,i] very
    # slightly) and zero the diagonal explicitly.
    mat = (mat + mat.T) / 2.0
    np.fill_diagonal(mat, 0.0)
    return names, mat


def classical_pcoa(dist_matrix: np.ndarray, n_components: int = 4) -> np.ndarray:
    """Classical (metric) multidimensional scaling / PCoA.

    Double-centers the squared distance matrix (the Gower transform) and
    eigendecomposes it. Mash distance is not a true Euclidean metric, so
    small negative eigenvalues are expected and simply dropped rather than
    treated as an error -- keep only the positive-eigenvalue components
    actually requested.
    """
    n = dist_matrix.shape[0]
    d2 = dist_matrix ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    keep = max(0, min(n_components, int(np.sum(eigvals > 1e-10))))
    if keep == 0:
        raise ValueError(
            "classical_pcoa: no positive eigenvalues -- the distance matrix "
            "carries no separable structure (e.g. all distances ~identical)"
        )
    coords = eigvecs[:, :keep] * np.sqrt(eigvals[:keep])
    return coords


def silhouette_score(coords: np.ndarray, labels: np.ndarray) -> float:
    """Mean silhouette coefficient over all points (no sklearn dependency).
    Returns -1.0 (worst possible) for a degenerate single-cluster labeling,
    so it never wins a k-selection sweep."""
    unique_labels = np.unique(labels)
    if len(unique_labels) < 2:
        return -1.0
    dist = squareform(pdist(coords))
    scores = []
    for i in range(len(coords)):
        own = labels[i]
        a_mask = (labels == own)
        a_mask[i] = False
        if not a_mask.any():
            scores.append(0.0)
            continue
        a = dist[i, a_mask].mean()
        b = min(
            dist[i, labels == other].mean()
            for other in unique_labels if other != own
        )
        scores.append((b - a) / max(a, b))
    return float(np.mean(scores))


def choose_k_and_cluster(
    coords: np.ndarray, k_min: int, k_max: int, seed: int = 0
) -> tuple[np.ndarray, int, dict[int, float]]:
    """Try k in [k_min, k_max], pick the k with the best silhouette score.
    Returns (labels_for_best_k, best_k, {k: silhouette_score})."""
    n = coords.shape[0]
    k_max = min(k_max, n - 1)
    scores: dict[int, float] = {}
    best_k, best_labels, best_score = None, None, -2.0
    rng = np.random.default_rng(seed)
    for k in range(k_min, k_max + 1):
        _, labels = kmeans2(coords, k, seed=rng, minit="++")
        score = silhouette_score(coords, labels)
        scores[k] = score
        if score > best_score:
            best_score, best_k, best_labels = score, k, labels
    return best_labels, best_k, scores


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument(
        "--groups", default="IN",
        help="comma-separated config.csv GROUP values to include (default: IN "
        "only -- an outgroup species mixed in dominates the Mash distance "
        "structure and swamps the within-species clade signal; see this "
        "script's module docstring)",
    )
    ap.add_argument("--k_range", default="2,10", help="min,max clusters to try")
    ap.add_argument(
        "--k_fixed", type=int, default=None,
        help="skip the silhouette sweep and cluster at exactly this k "
        "(e.g. to match a known number of ground-truth clusters)",
    )
    ap.add_argument("--n_pcoa_components", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from novinvenio_path import add_novinvenio_lib_to_path  # noqa: E402

    add_novinvenio_lib_to_path()
    from config_parser import parse_config  # noqa: E402

    wanted_groups = set(args.groups.split(","))
    samples = [s for s in parse_config(args.config) if s.group in wanted_groups]
    data_dir = Path(args.data_dir)
    dna_paths = {s.short: data_dir / "dna" / s.dna for s in samples if s.dna}

    sketch_prefix = "clade_sketches"
    subprocess.run(
        ["mash", "sketch", "-o", sketch_prefix] + [str(p) for p in dna_paths.values()],
        check=True,
    )
    dist_out = subprocess.run(
        ["mash", "dist", "-t", f"{sketch_prefix}.msh", f"{sketch_prefix}.msh"],
        check=True, capture_output=True, text=True,
    ).stdout
    path_names, dist_matrix = parse_mash_dist_matrix(dist_out.splitlines())

    path_to_short = {str(p): short for short, p in dna_paths.items()}
    shorts = [path_to_short[p] for p in path_names]

    coords = classical_pcoa(dist_matrix, n_components=args.n_pcoa_components)
    if args.k_fixed is not None:
        rng = np.random.default_rng(args.seed)
        _, labels = kmeans2(coords, args.k_fixed, seed=rng, minit="++")
        best_k = args.k_fixed
        print(f"Fixed k={best_k} (silhouette sweep skipped)", file=sys.stderr)
    else:
        k_min, k_max = (int(x) for x in args.k_range.split(","))
        labels, best_k, scores = choose_k_and_cluster(coords, k_min, k_max, seed=args.seed)
        print(
            f"Selected k={best_k} (silhouette={scores[best_k]:.3f}); "
            f"scores by k: {scores}",
            file=sys.stderr,
        )

    with open(args.output, "w") as fh:
        header = ["Short", "mash_clade"] + [f"pcoa{i+1}" for i in range(coords.shape[1])]
        fh.write("\t".join(header) + "\n")
        for short, label, coord in zip(shorts, labels, coords):
            row = [short, f"clade_{label}"] + [f"{v:.6f}" for v in coord]
            fh.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
