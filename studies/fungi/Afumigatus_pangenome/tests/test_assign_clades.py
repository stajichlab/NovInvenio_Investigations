import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from assign_clades import (
    parse_mash_dist_matrix,
    classical_pcoa,
    silhouette_score,
    choose_k_and_cluster,
)


def test_parse_mash_dist_matrix_symmetrizes_and_zeroes_diagonal():
    lines = [
        "#query\ts1\ts2",
        "s1\t0.0000\t0.0500",
        "s2\t0.0502\t0.0000",
    ]
    names, mat = parse_mash_dist_matrix(lines)
    assert names == ["s1", "s2"]
    assert mat[0, 0] == 0.0 and mat[1, 1] == 0.0
    assert mat[0, 1] == mat[1, 0]
    assert abs(mat[0, 1] - 0.0501) < 1e-9


def test_parse_mash_dist_matrix_empty_input():
    names, mat = parse_mash_dist_matrix([])
    assert names == []
    assert mat.shape == (0, 0)


def test_classical_pcoa_recovers_two_well_separated_clusters():
    # Two tight clusters far apart in a toy Euclidean distance matrix
    points = np.array([[0, 0], [0.1, 0], [10, 10], [10.1, 10]])
    dist = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(-1))
    coords = classical_pcoa(dist, n_components=2)
    # Points 0,1 should be close in PCoA space; 2,3 close; the two pairs far apart
    d01 = np.linalg.norm(coords[0] - coords[1])
    d23 = np.linalg.norm(coords[2] - coords[3])
    d02 = np.linalg.norm(coords[0] - coords[2])
    assert d01 < d02 and d23 < d02


def test_silhouette_score_degenerate_single_cluster_returns_worst():
    coords = np.array([[0, 0], [1, 1], [2, 2]])
    labels = np.array([0, 0, 0])
    assert silhouette_score(coords, labels) == -1.0


def test_choose_k_and_cluster_prefers_true_cluster_count():
    rng = np.random.default_rng(0)
    cluster_a = rng.normal(loc=[0, 0], scale=0.05, size=(10, 2))
    cluster_b = rng.normal(loc=[5, 5], scale=0.05, size=(10, 2))
    coords = np.vstack([cluster_a, cluster_b])
    labels, best_k, scores = choose_k_and_cluster(coords, k_min=2, k_max=5, seed=0)
    assert best_k == 2
    assert len(np.unique(labels)) == 2
