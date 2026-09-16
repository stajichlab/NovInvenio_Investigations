import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from find_accessory_islands import (
    load_is_core, load_strain_gene_orders, load_significant_physical_pairs,
    load_hit_families, find_significant_islands, build_pair_index, PHYSICAL_CLASSIFICATIONS,
)


def test_load_is_core_reads_bin_column(tmp_path):
    path = tmp_path / "freq.tsv"
    path.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famA\t1.0\t10\tcore\n"
        "famB\t0.9\t9\tsoft_core\n"
        "famC\t0.1\t1\tcloud\n"
    )
    is_core = load_is_core(str(path))
    assert is_core == {"famA": True, "famB": True, "famC": False}


def test_load_strain_gene_orders_sorts_by_rank(tmp_path):
    path = tmp_path / "positions.tsv"
    path.write_text(
        "Short\tfamily\tcontig\trank\n"
        "s1\tfamB\tc1\t2\n"
        "s1\tfamA\tc1\t0\n"
        "s1\tfamC\tc1\t1\n"
        "s2\tfamX\tc1\t0\n"
    )
    orders = load_strain_gene_orders(str(path))
    assert [g[0] for g in orders["s1"]] == ["famA", "famC", "famB"]
    assert [g[0] for g in orders["s2"]] == ["famX"]


def test_load_significant_physical_pairs_filters_classification(tmp_path):
    path = tmp_path / "pc.tsv"
    path.write_text(
        "family_a\tfamily_b\tclassification\tlinkage_fraction\tjaccard\t"
        "fisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n"
        "famA\tfamB\tstarship_explained\t1.0\t1.0\t1e-5\t1e-4\t0.0\tloss\t{}\n"
        "famC\tfamD\ttrans\t0.0\t0.5\t1e-5\t1e-4\t0.0\tloss\t{}\n"
        "famE\tfamF\tunexplained_physical\t1.0\t1.0\t1e-5\t1e-4\t0.0\tloss\t{}\n"
    )
    pairs = load_significant_physical_pairs(str(path))
    assert frozenset({"famA", "famB"}) in pairs
    assert frozenset({"famE", "famF"}) in pairs
    assert frozenset({"famC", "famD"}) not in pairs
    assert PHYSICAL_CLASSIFICATIONS == {"starship_explained", "unexplained_physical", "ambiguous_linkage"}


def test_load_hit_families_resolves_via_cluster_membership(tmp_path):
    path = tmp_path / "hits.tblout"
    path.write_text("s1|P1.1 - queryname - 1e-10 100.0 0.0 ...\n# comment\n")
    member_to_rep = {"s1|P1.1": "famA"}
    hits = load_hit_families(str(path), member_to_rep)
    assert hits == {"famA"}


def test_find_significant_islands_only_keeps_islands_with_a_significant_pair():
    # Strain s1: famA-famB-famC all non-core and adjacent -> one island.
    # famA/famB is a significant pair -> island kept.
    is_core = {"famA": False, "famB": False, "famC": False, "famCore": True}
    strain_gene_orders = {
        "s1": [
            ("famCore", "c1", 0, 0),
            ("famA", "c1", 1, 1),
            ("famB", "c1", 2, 2),
            ("famC", "c1", 3, 3),
            ("famCore", "c1", 4, 4),
        ],
    }
    significant_pairs = {frozenset({"famA", "famB"}): "unexplained_physical"}
    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    assert len(islands) == 1
    assert islands[0]["members"] == ["famA", "famB", "famC"]
    assert islands[0]["supporting_pairs"] == [
        (frozenset({"famA", "famB"}), "unexplained_physical")
    ]


def test_find_significant_islands_drops_islands_with_no_significant_pair():
    is_core = {"famA": False, "famB": False}
    strain_gene_orders = {"s1": [("famA", "c1", 0, 0), ("famB", "c1", 1, 1)]}
    # famA/famB adjacent but NOT in significant_pairs
    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs={})
    assert islands == []


def test_build_pair_index_indexes_by_both_family_members():
    significant_pairs = {
        frozenset({"famA", "famB"}): "unexplained_physical",
        frozenset({"famB", "famC"}): "starship_explained",
    }
    index = build_pair_index(significant_pairs)
    assert set(index.keys()) == {"famA", "famB", "famC"}
    assert (frozenset({"famA", "famB"}), "unexplained_physical") in index["famA"]
    assert len(index["famB"]) == 2  # famB participates in both pairs


def test_find_significant_islands_index_path_matches_pairs_with_a_hub_family():
    # famB is a "hub" appearing in multiple significant pairs; only the one
    # fully contained in the island (famA-famB) should surface, not the
    # famB-famZ pair whose other member isn't in this island at all.
    is_core = {"famA": False, "famB": False, "famC": False}
    strain_gene_orders = {
        "s1": [("famA", "c1", 0, 0), ("famB", "c1", 1, 1), ("famC", "c1", 2, 2)],
    }
    significant_pairs = {
        frozenset({"famA", "famB"}): "unexplained_physical",
        frozenset({"famB", "famZ"}): "starship_explained",  # famZ not in this island
    }
    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    assert len(islands) == 1
    assert islands[0]["supporting_pairs"] == [
        (frozenset({"famA", "famB"}), "unexplained_physical")
    ]


def test_find_significant_islands_requires_pair_fully_within_one_island():
    # famA and famC are both non-core but separated by a core gene ->
    # two separate single-gene islands, neither containing the whole pair.
    is_core = {"famA": False, "famCore": True, "famC": False}
    strain_gene_orders = {
        "s1": [("famA", "c1", 0, 0), ("famCore", "c1", 1, 1), ("famC", "c1", 2, 2)],
    }
    significant_pairs = {frozenset({"famA", "famC"}): "ambiguous_linkage"}
    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    assert islands == []
