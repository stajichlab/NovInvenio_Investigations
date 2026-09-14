import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT
from cooccurrence import (
    jaccard, fisher_pvalue, benjamini_hochberg, polarize_direction,
    clade_composition, permutation_null_pvalue, find_cooccurring_pairs,
    presence_bitset, fisher_counts_from_bits, fisher_pvalue_from_counts,
    quick_screen_pvalue, benjamini_hochberg_sparse,
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


def test_presence_bitset_matches_presence_vector_within():
    strains = [f"s{i}" for i in range(6)]
    pm = PresenceMatrix(families=["famA"], strains=strains)
    for s in ["s1", "s3", "s5"]:
        pm.set_call("famA", s, PRESENT)
    from cooccurrence import presence_vector_within

    bits = presence_bitset(pm, "famA", strains)
    vec = presence_vector_within(pm, "famA", strains)
    assert [bool((bits >> i) & 1) for i in range(len(strains))] == vec


def test_fisher_counts_from_bits_matches_naive_list_counts():
    # a: strains 0,1,2 present; b: strains 1,2,3 present (over 5 strains)
    a_bits = 0b00111
    b_bits = 0b01110
    both, a_only, b_only, neither = fisher_counts_from_bits(a_bits, b_bits, n=5)
    # both = {1,2} -> 2; a_only = {0} -> 1; b_only = {3} -> 1; neither = {4} -> 1
    assert (both, a_only, b_only, neither) == (2, 1, 1, 1)


def test_fisher_pvalue_from_counts_matches_list_based_fisher_pvalue():
    a = [True, True, True, True, False, False, False, False]
    b = [True, True, True, True, False, False, False, False]
    both, a_only, b_only, neither = fisher_counts_from_bits(0b00001111, 0b00001111, n=8)
    p_from_counts = fisher_pvalue_from_counts(both, a_only, b_only, neither)
    p_from_lists = fisher_pvalue(a, b)
    assert p_from_counts == p_from_lists


def test_quick_screen_pvalue_never_screens_out_a_clearly_significant_table():
    # Perfect co-occurrence in 6 of 10 strains -- the real find_cooccurring_pairs
    # regression fixture below; must clear a generous screen_alpha comfortably.
    both, a_only, b_only, neither = 6, 0, 0, 4
    screen_p = quick_screen_pvalue(both, a_only, b_only, neither)
    exact_p = fisher_pvalue_from_counts(both, a_only, b_only, neither)
    assert screen_p < 0.2
    assert exact_p < 0.05


def test_quick_screen_pvalue_screens_out_a_clearly_random_table():
    # Roughly independent presence (both ~ expected under independence) --
    # should screen out well above any reasonable screen_alpha.
    both, a_only, b_only, neither = 3, 3, 3, 3
    assert quick_screen_pvalue(both, a_only, b_only, neither) > 0.2


def test_quick_screen_pvalue_handles_degenerate_zero_margin():
    assert quick_screen_pvalue(0, 0, 0, 10) == 1.0
    assert quick_screen_pvalue(10, 0, 0, 0) == 1.0


def test_benjamini_hochberg_sparse_matches_dense_when_total_m_equals_len():
    pvals = [0.001, 0.01, 0.02, 0.5, 0.9]
    dense = benjamini_hochberg(pvals)
    sparse = benjamini_hochberg_sparse(pvals, total_m=len(pvals))
    for d, s in zip(dense, sparse):
        assert abs(d - s) < 1e-9


def test_benjamini_hochberg_sparse_is_more_conservative_with_larger_total_m():
    # Same small p-values, but told there were many more (screened-out, p=1)
    # hypotheses in the true test space -- q-values must not get SMALLER.
    pvals = [0.001, 0.002, 0.003]
    q_small_m = benjamini_hochberg_sparse(pvals, total_m=3)
    q_large_m = benjamini_hochberg_sparse(pvals, total_m=1_000_000)
    for small, large in zip(q_small_m, q_large_m):
        assert large >= small


def test_benjamini_hochberg_sparse_rejects_total_m_smaller_than_survivor_count():
    try:
        benjamini_hochberg_sparse([0.01, 0.02], total_m=1)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_find_cooccurring_pairs_scaling_rewrite_matches_expected_calls():
    """Same fixture as test_find_cooccurring_pairs_applies_frequency_floor_and_fdr,
    re-run with a tight screen_alpha to confirm the prefilter doesn't discard
    a real signal, and with screening effectively disabled (screen_alpha=1.0)
    to confirm both code paths agree."""
    strains = [f"s{i}" for i in range(10)]
    pm = PresenceMatrix(families=["famA", "famB", "famRare1", "famRare2"], strains=strains)
    for s in strains[:6]:
        pm.set_call("famA", s, PRESENT)
        pm.set_call("famB", s, PRESENT)
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

    for screen_alpha in (0.05, 1.0):
        pairs = find_cooccurring_pairs(
            pm, frequency_table, clade_of_strain, outgroup_presence,
            min_strain_count=5, fdr_alpha=0.05, n_perms=100, seed=0,
            screen_alpha=screen_alpha,
        )
        reported = {(p["family_a"], p["family_b"]) for p in pairs}
        assert ("famA", "famB") in reported or ("famB", "famA") in reported, (
            f"screen_alpha={screen_alpha} wrongly discarded the real signal"
        )


def test_find_cooccurring_pairs_screen_alpha_never_tighter_than_fdr_alpha():
    """Regression test: at tiny sample sizes (n=3) the normal-approximation
    screen is a poor approximation and can wrongly discard a pair a
    permissive fdr_alpha would otherwise keep -- caught by the pipeline
    integration test's fdr_alpha=1.0 ("keep everything") fixture. The screen
    must auto-widen to at least fdr_alpha regardless of its own default."""
    strains = ["s1", "s2", "s3"]
    pm = PresenceMatrix(families=["famA", "famB"], strains=strains)
    for s in ["s1", "s2"]:
        pm.set_call("famA", s, PRESENT)
        pm.set_call("famB", s, PRESENT)
    frequency_table = [
        {"family": "famA", "bin": "shell"}, {"family": "famB", "bin": "shell"},
    ]
    clade_of_strain = {s: "cladeX" for s in strains}
    outgroup_presence = {"famA": (0, 0), "famB": (0, 0)}

    pairs = find_cooccurring_pairs(
        pm, frequency_table, clade_of_strain, outgroup_presence,
        min_strain_count=2, fdr_alpha=1.0, n_perms=20, strains=strains,
        screen_alpha=0.2,  # the default -- would wrongly discard this pair
    )
    reported = {(p["family_a"], p["family_b"]) for p in pairs}
    assert ("famA", "famB") in reported or ("famB", "famA") in reported


def test_main_treats_outgroup_absent_from_matrix_as_ambiguous_not_gain(
    monkeypatch, capsys, tmp_path
):
    """Regression test for the final-review fix-wave bug: when a matrix is
    built ingroup-only (e.g. build_presence_matrix.py --groups IN), main()
    must not compute outgroup_total from config.csv's full OUT-group list --
    that made every family's polarize_direction(0, N) come back "loss" or
    "gain" with total confidence, when the correct call is "ambiguous"
    (no outgroup data is actually present to check)."""
    from cooccurrence import main

    ingroup = [f"s{i}" for i in range(1, 11)]
    config_lines = ["GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup"]
    for i, s in enumerate(ingroup, start=1):
        clade = "cladeA" if i % 2 else "cladeB"
        config_lines.append(f"IN,Aspergillus fumigatus,{s},{s}.fa,{s}.dna.fa,{s},{clade}")
    config_lines.append("OUT,Aspergillus lentulus,out1,out1.fa,out1.dna.fa,out1,")
    config_file = tmp_path / "config.csv"
    config_file.write_text("\n".join(config_lines) + "\n")

    # Matrix built ingroup-only: no "out1" column at all. famA/famB co-occur
    # in 8 of 10 strains (a large enough margin to clear Fisher/BH at N=10).
    pm = PresenceMatrix(families=["famA", "famB"], strains=ingroup)
    for s in ingroup[:8]:
        pm.set_call("famA", s, PRESENT)
        pm.set_call("famB", s, PRESENT)
    matrix_file = tmp_path / "matrix.tsv"
    pm.to_tsv(matrix_file)

    freq_file = tmp_path / "frequency_table.tsv"
    freq_file.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famA\t0.8000\t8\tshell\n"
        "famB\t0.8000\t8\tshell\n"
    )

    output_file = tmp_path / "cooccurring_pairs.tsv"

    monkeypatch.setattr(sys, "argv", [
        "cooccurrence.py",
        "--matrix", str(matrix_file),
        "--frequency_table", str(freq_file),
        "--config", str(config_file),
        "--min_strain_count", "5",
        "--n_perms", "50",
        "--output", str(output_file),
    ])
    main()

    err = capsys.readouterr().err
    assert "none of config.csv's OUT-group strains are columns" in err

    rows = output_file.read_text().splitlines()[1:]
    assert rows, "expected the strongly co-occurring famA/famB pair to be reported"
    direction = rows[0].split("\t")[6]
    assert direction == "ambiguous", (
        f"expected 'ambiguous' when the matrix has no outgroup columns, got {direction!r}"
    )
