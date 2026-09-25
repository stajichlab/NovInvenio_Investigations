import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from k_sensitivity_sweep import (
    linkage_fraction_from_distances,
    per_strain_min_distances,
    reclassify,
)


def test_per_strain_min_distances_picks_minimum_across_multi_copy_pairs():
    gene_position = {
        "s1": {"famA": [("c1", 0), ("c1", 100)], "famB": [("c1", 5)]},
        "s2": {"famA": [("c1", 0)], "famB": [("c2", 0)]},  # different contig -> None
        "s3": {"famA": [("c1", 0)]},  # missing famB -> not counted at all
    }
    dists = per_strain_min_distances("famA", "famB", gene_position)
    assert sorted(dists, key=lambda d: (d is None, d)) == [5, None]


def test_linkage_fraction_from_distances_counts_within_window_only():
    dists = [3, 15, None, 8]
    assert linkage_fraction_from_distances(dists, k=10) == 0.5  # 3 and 8 qualify, 15 and None don't
    assert linkage_fraction_from_distances(dists, k=20) == 0.75  # 3, 15, 8 qualify


def test_linkage_fraction_from_distances_empty_is_zero():
    assert linkage_fraction_from_distances([], k=10) == 0.0


def test_reclassify_trans_stays_trans_below_trans_threshold():
    label = reclassify(
        frac=0.0, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position={}, captain_families={}, k=10,
    )
    assert label == "trans"


def test_reclassify_trans_becomes_unexplained_physical_once_frac_crosses_threshold():
    # frac >= DEFAULT_PHYSICAL_THRESHOLD (0.5) and no captain-gene evidence
    label = reclassify(
        frac=0.6, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position={}, captain_families={}, k=50,
    )
    assert label == "unexplained_physical"


def test_reclassify_becomes_starship_explained_when_captain_evidence_present():
    gene_position = {"s1": {"famA": [("c1", 5)], "famCaptain": [("c1", 8)]}}
    captain_families = {"s1": {"famCaptain"}}
    label = reclassify(
        frac=0.6, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position=gene_position, captain_families=captain_families, k=10,
    )
    assert label == "starship_explained"


def test_reclassify_uses_swept_k_for_captain_check_when_captain_k_not_given():
    # captain gene is 15 genes away -- outside k=10 but inside k=20
    gene_position = {"s1": {"famA": [("c1", 5)], "famCaptain": [("c1", 20)]}}
    captain_families = {"s1": {"famCaptain"}}
    at_k10 = reclassify(
        frac=0.6, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position=gene_position, captain_families=captain_families, k=10,
    )
    at_k20 = reclassify(
        frac=0.6, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position=gene_position, captain_families=captain_families, k=20,
    )
    assert at_k10 == "unexplained_physical"
    assert at_k20 == "starship_explained"


def test_reclassify_fixed_captain_k_does_not_grow_with_swept_k():
    # same captain gene 15 genes away, but captain_k pinned to 10 regardless
    # of how large the swept k gets -- must stay unexplained_physical
    gene_position = {"s1": {"famA": [("c1", 5)], "famCaptain": [("c1", 20)]}}
    captain_families = {"s1": {"famCaptain"}}
    label = reclassify(
        frac=0.6, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position=gene_position, captain_families=captain_families,
        k=200, captain_k=10,
    )
    assert label == "unexplained_physical"


def test_reclassify_ambiguous_linkage_band():
    label = reclassify(
        frac=0.2, permutation_p=0.001, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position={}, captain_families={}, k=10,
    )
    assert label == "ambiguous_linkage"


def test_reclassify_trans_unconfirmed_when_permutation_fails():
    label = reclassify(
        frac=0.0, permutation_p=0.5, n_clades=3,
        family_a="famA", family_b="famB",
        gene_position={}, captain_families={}, k=10,
    )
    assert label == "trans_unconfirmed"


def test_reclassify_trans_unconfirmed_when_too_few_clades():
    label = reclassify(
        frac=0.0, permutation_p=0.001, n_clades=1,
        family_a="famA", family_b="famB",
        gene_position={}, captain_families={}, k=10,
    )
    assert label == "trans_unconfirmed"
