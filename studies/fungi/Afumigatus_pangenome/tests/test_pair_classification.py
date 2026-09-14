import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pair_classification import (
    load_family_positions, load_captain_families, has_captain_evidence,
    classify_pair,
)


def test_load_family_positions_groups_multi_copy_per_strain(tmp_path):
    f = tmp_path / "family_positions.tsv"
    f.write_text(
        "Short\tfamily\tcontig\trank\n"
        "s1\tfamA\tc1\t0\n"
        "s1\tfamA\tc1\t50\n"
        "s1\tfamB\tc1\t2\n"
    )
    result = load_family_positions(f)
    assert result == {"s1": {"famA": [("c1", 0), ("c1", 50)], "famB": [("c1", 2)]}}


def test_load_captain_families_resolves_hit_protein_to_its_family(tmp_path):
    tblout = tmp_path / "captain.tblout"
    tblout.write_text(
        "# comment\n"
        "s1|P1.1 -          DUF3435  PF11917.14   1e-40  100.0   0.3   1e-40  100.0   0.2   1   1   0   0   1   1   0   0 -\n"
    )
    member_to_rep = {"s1|P1.1": "famCaptain"}
    result = load_captain_families(str(tblout), member_to_rep)
    assert result == {"s1": {"famCaptain"}}


def test_has_captain_evidence_true_when_within_k():
    gene_position = {"s1": {"famA": [("c1", 5)], "famCaptain": [("c1", 8)]}}
    captain_families = {"s1": {"famCaptain"}}
    assert has_captain_evidence("famA", "famB", gene_position, captain_families, k=10) is True


def test_has_captain_evidence_false_when_far():
    gene_position = {"s1": {"famA": [("c1", 5)], "famCaptain": [("c1", 900)]}}
    captain_families = {"s1": {"famCaptain"}}
    assert has_captain_evidence("famA", "famB", gene_position, captain_families, k=10) is False


def test_classify_pair_insufficient_data_below_min_co_carrying():
    gene_position = {
        "s1": {"famA": [("c1", 5)], "famB": [("c1", 7)]},
    }
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.01, clade_composition={"cladeA": 1},
        gene_position=gene_position, captain_families={}, min_co_carrying=5,
    )
    assert label == "insufficient_data"


def test_classify_pair_starship_explained_when_physical_and_captain_nearby():
    strains = [f"s{i}" for i in range(6)]
    gene_position = {
        s: {"famA": [("c1", 5)], "famB": [("c1", 7)], "famCaptain": [("c1", 6)]}
        for s in strains
    }
    captain_families = {s: {"famCaptain"} for s in strains}
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.5, clade_composition={"cladeA": 6},
        gene_position=gene_position, captain_families=captain_families,
        min_co_carrying=5,
    )
    assert label == "starship_explained"
    assert frac == 1.0


def test_classify_pair_unexplained_physical_when_linked_but_no_captain():
    strains = [f"s{i}" for i in range(6)]
    gene_position = {s: {"famA": [("c1", 5)], "famB": [("c1", 7)]} for s in strains}
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.5, clade_composition={"cladeA": 6},
        gene_position=gene_position, captain_families={}, min_co_carrying=5,
    )
    assert label == "unexplained_physical"


def test_classify_pair_trans_when_unlinked_and_permutation_and_clades_clear():
    strains = [f"s{i}" for i in range(6)]
    gene_position = {
        s: {"famA": [("c1", 5)], "famB": [("c1", 900)]} for s in strains
    }
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.01,
        clade_composition={"cladeA": 3, "cladeB": 3},
        gene_position=gene_position, captain_families={}, min_co_carrying=5,
    )
    assert label == "trans"
    assert frac == 0.0


def test_classify_pair_trans_unconfirmed_when_permutation_fails():
    strains = [f"s{i}" for i in range(6)]
    gene_position = {
        s: {"famA": [("c1", 5)], "famB": [("c1", 900)]} for s in strains
    }
    # permutation_p=0.5 -- fails the perm_alpha=0.05 gate despite clearing FDR
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.5,
        clade_composition={"cladeA": 3, "cladeB": 3},
        gene_position=gene_position, captain_families={}, min_co_carrying=5,
    )
    assert label == "trans_unconfirmed"


def test_classify_pair_trans_unconfirmed_when_single_clade_only():
    strains = [f"s{i}" for i in range(6)]
    gene_position = {
        s: {"famA": [("c1", 5)], "famB": [("c1", 900)]} for s in strains
    }
    # Only one clade despite low permutation_p -- fails the >=2-clade gate.
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.01,
        clade_composition={"cladeA": 6},
        gene_position=gene_position, captain_families={}, min_co_carrying=5,
    )
    assert label == "trans_unconfirmed"


def test_classify_pair_ambiguous_linkage_in_the_middle_band():
    # 6 co-carrying strains, 3 linked -> frac=0.5... use 10 strains, 4 linked
    # -> frac=0.4, strictly between trans_threshold(0.05) and
    # physical_threshold(0.5).
    strains = [f"s{i}" for i in range(10)]
    gene_position = {}
    for i, s in enumerate(strains):
        if i < 4:
            gene_position[s] = {"famA": [("c1", 5)], "famB": [("c1", 7)]}  # linked
        else:
            gene_position[s] = {"famA": [("c1", 5)], "famB": [("c1", 900)]}  # far
    label, frac = classify_pair(
        "famA", "famB", permutation_p=0.01, clade_composition={"cladeA": 10},
        gene_position=gene_position, captain_families={}, min_co_carrying=5,
    )
    assert label == "ambiguous_linkage"
    assert frac == 0.4
