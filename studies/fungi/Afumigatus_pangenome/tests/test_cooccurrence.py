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
    assert polarize_direction(outgroup_present_count=2, outgroup_total=2) == "loss"
    assert polarize_direction(outgroup_present_count=0, outgroup_total=2) == "gain"


def test_polarize_direction_ambiguous_when_outgroup_strains_disagree():
    # One outgroup strain carries the family, the other doesn't -- the
    # outgroup itself is split, so polarity can't be called either way.
    assert polarize_direction(outgroup_present_count=1, outgroup_total=2) == "ambiguous"
    # No outgroup strains to check at all is likewise unpolarizable.
    assert polarize_direction(outgroup_present_count=0, outgroup_total=0) == "ambiguous"


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


def test_permutation_null_pvalue_never_reports_exactly_zero():
    # Phipson & Smyth (2010) correction: even a perfectly co-occurring pair
    # that beats every single permutation must report (0+1)/(n_perms+1),
    # never a false-precision 0.0.
    a = [True, True, True, True, False, False, False, False]
    b = list(a)
    clades = ["Clade_1"] * 8
    rng = random.Random(0)
    p = permutation_null_pvalue(a, b, clades, n_perms=10, rng=rng)
    assert p == 1 / 11
    assert p > 0.0


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
    outgroup_presence = {
        "famA": (0, 0), "famB": (0, 0), "famRare1": (0, 0), "famRare2": (0, 0),
    }

    pairs = find_cooccurring_pairs(
        pm, frequency_table, clade_of_strain, outgroup_presence,
        min_strain_count=5, fdr_alpha=0.05, n_perms=100, seed=0,
    )
    reported = {(p["family_a"], p["family_b"]) for p in pairs}
    assert ("famA", "famB") in reported or ("famB", "famA") in reported
    assert ("famRare1", "famRare2") not in reported and ("famRare2", "famRare1") not in reported


def test_find_cooccurring_pairs_ignores_outgroup_strains_for_the_floor():
    # 2 ingroup strains + 1 outgroup strain. famX is present in only 1
    # ingroup strain plus the outgroup strain (3 total presences in the
    # full matrix), but must be excluded by an ingroup-only floor of 2 --
    # the outgroup strain's presence must not count toward the floor.
    strains = ["in1", "in2", "out1"]
    pm = PresenceMatrix(families=["famX", "famY"], strains=strains)
    pm.set_call("famX", "in1", PRESENT)
    pm.set_call("famX", "out1", PRESENT)
    # famY clears the ingroup-only floor of 2 on its own (both ingroup
    # strains carry it) so there's a valid partner for famX to be tested
    # against if it wrongly qualified.
    pm.set_call("famY", "in1", PRESENT)
    pm.set_call("famY", "in2", PRESENT)

    frequency_table = [
        {"family": "famX", "bin": "shell"},
        {"family": "famY", "bin": "shell"},
    ]
    clade_of_strain = {s: "Clade_1" for s in strains}
    outgroup_presence = {"famX": (1, 1), "famY": (0, 1)}

    pairs = find_cooccurring_pairs(
        pm, frequency_table, clade_of_strain, outgroup_presence,
        min_strain_count=2, fdr_alpha=0.05, n_perms=50, seed=0,
        strains=["in1", "in2"],
    )
    reported_families = {fam for p in pairs for fam in (p["family_a"], p["family_b"])}
    assert "famX" not in reported_families


def test_warn_if_unstratified_fires_on_mostly_empty_taxon_groups(capsys):
    from cooccurrence import warn_if_unstratified, unlabelled_clade_fraction

    strains = ["s1", "s2", "s3", "s4"]
    clades = {"s1": "cladeA", "s2": "", "s3": "", "s4": ""}
    assert unlabelled_clade_fraction(strains, clades) == 0.75
    assert warn_if_unstratified(strains, clades) is True
    err = capsys.readouterr().err
    assert "UNSTRATIFIED" in err


def test_warn_if_unstratified_silent_when_clades_are_labelled(capsys):
    from cooccurrence import warn_if_unstratified

    clades = {"s1": "cladeA", "s2": "cladeB", "s3": "cladeA"}
    assert warn_if_unstratified(["s1", "s2", "s3"], clades) is False
    assert capsys.readouterr().err == ""
