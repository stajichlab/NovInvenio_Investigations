import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT
from cooccurrence import (
    jaccard, fisher_pvalue, benjamini_hochberg, polarize_direction,
    clade_composition, permutation_null_pvalue, find_cooccurring_pairs,
)


def test_jaccard_identical_and_disjoint():
    assert jaccard([True, True, False], [True, True, False]) == 1.0
    assert jaccard([True, False], [False, True]) == 0.0


def test_fisher_pvalue_significant_for_perfect_cooccurrence():
    a = [True, True, True, True, False, False, False, False]
    b = [True, True, True, True, False, False, False, False]
    p = fisher_pvalue(a, b)
    assert p < 0.05


def test_benjamini_hochberg_orders_correctly():
    pvals = [0.01, 0.02, 0.03, 0.9]
    qvals = benjamini_hochberg(pvals)
    assert len(qvals) == 4
    assert qvals[0] <= qvals[1] <= qvals[2]
    assert qvals[3] >= qvals[2]


def test_polarize_direction():
    assert polarize_direction(present_in_outgroup=True, freq_in_ingroup=0.6) == "loss"
    assert polarize_direction(present_in_outgroup=False, freq_in_ingroup=0.6) == "gain"


def test_clade_composition_counts_per_clade():
    clade_of_strain = {"s1": "Clade_1", "s2": "Clade_1", "s3": "Clade_2"}
    counts = clade_composition(["s1", "s2", "s3"], clade_of_strain)
    assert counts == {"Clade_1": 2, "Clade_2": 1}


def test_permutation_null_pvalue_high_when_confound_is_purely_clade():
    # a and b are both simply "is this strain in Clade_1" -- co-occurrence
    # is entirely explained by clade membership, so the within-clade
    # permutation null should NOT find this significant.
    clade_of_strain = {f"s{i}": ("Clade_1" if i < 4 else "Clade_2") for i in range(8)}
    a = [clade_of_strain[f"s{i}"] == "Clade_1" for i in range(8)]
    b = list(a)
    rng = random.Random(0)
    p = permutation_null_pvalue(a, b, [clade_of_strain[f"s{i}"] for i in range(8)], n_perms=200, rng=rng)
    assert p > 0.05


def test_find_cooccurring_pairs_applies_frequency_floor_and_fdr():
    strains = [f"s{i}" for i in range(10)]
    pm = PresenceMatrix(families=["famA", "famB", "famRare1", "famRare2"], strains=strains)
    # famA/famB co-occur in strains 0-5 (6 strains, clears the floor)
    for s in strains[:6]:
        pm.set_call("famA", s, PRESENT)
        pm.set_call("famB", s, PRESENT)
    # famRare1/famRare2 co-occur only in strains 0-1 (below the floor of 5)
    pm.set_call("famRare1", "s0", PRESENT)
    pm.set_call("famRare1", "s1", PRESENT)
    pm.set_call("famRare2", "s0", PRESENT)
    pm.set_call("famRare2", "s1", PRESENT)

    frequency_table = [
        {"family": "famA", "bin": "shell"}, {"family": "famB", "bin": "shell"},
        {"family": "famRare1", "bin": "cloud"}, {"family": "famRare2", "bin": "cloud"},
    ]
    clade_of_strain = {s: "Clade_1" for s in strains}
    outgroup_presence = {"famA": False, "famB": False, "famRare1": False, "famRare2": False}

    pairs = find_cooccurring_pairs(
        pm, frequency_table, clade_of_strain, outgroup_presence,
        min_strain_count=5, fdr_alpha=0.05, n_perms=100, seed=0,
    )
    reported = {(p["family_a"], p["family_b"]) for p in pairs}
    assert ("famA", "famB") in reported or ("famB", "famA") in reported
    assert ("famRare1", "famRare2") not in reported and ("famRare2", "famRare1") not in reported
