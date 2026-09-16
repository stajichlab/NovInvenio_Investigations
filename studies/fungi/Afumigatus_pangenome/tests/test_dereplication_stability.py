import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from dereplication_stability import (
    load_pairwise_mash,
    load_ground_truth_clusters,
    load_path_to_short,
    single_linkage_groups,
    hierarchical_groups,
    groups_to_labels,
    ari_against_ground_truth,
)


def test_load_pairwise_mash_parses_and_skips_warn_lines(tmp_path):
    path = tmp_path / "dist.tsv"
    path.write_text(
        " WARN some cache message\n"
        "/data/a.fa\t/data/b.fa\t0.001\t0\t500/1000\n"
        "/data/a.fa\t/data/a.fa\t0\t0\t1000/1000\n"
    )
    pairs = load_pairwise_mash(str(path))
    assert pairs == [("/data/a.fa", "/data/b.fa", 0.001)]


def test_load_ground_truth_clusters(tmp_path):
    path = tmp_path / "s21_matches.tsv"
    path.write_text(
        "Short\tStrain\tS21_groupID\tcurrent_TaxonGroup\n"
        "Asfu_A\tA\tcluster1\tMash_clade_0\n"
        "Asfu_B\tB\tcluster2\tMash_clade_0\n"
    )
    truth = load_ground_truth_clusters(str(path))
    assert truth == {"Asfu_A": "cluster1", "Asfu_B": "cluster2"}


def test_load_path_to_short(tmp_path):
    path = tmp_path / "config.csv"
    path.write_text(
        "GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
        "IN,Aspergillus fumigatus,X,x.pep.fa,x.dna.fa,x.gff3,Asfu_X,\n"
    )
    mapping = load_path_to_short(str(path), data_dir="/data_dir")
    assert mapping == {"/data_dir/dna/x.dna.fa": "Asfu_X"}


def test_single_linkage_groups_chains_transitively():
    # a-b close, b-c close, a-c far: single-linkage still merges all three
    # (the "chaining" property under test).
    pairs = [("a", "b", 0.0005), ("b", "c", 0.0005), ("a", "c", 0.01)]
    groups = single_linkage_groups(pairs, threshold=0.001, nodes={"a", "b", "c"})
    assert groups == [{"a", "b", "c"}]


def test_single_linkage_groups_respects_threshold():
    pairs = [("a", "b", 0.0005), ("b", "c", 0.01)]
    groups = single_linkage_groups(pairs, threshold=0.001, nodes={"a", "b", "c"})
    assert {"a", "b"} in groups
    assert {"c"} in groups


def test_hierarchical_groups_complete_linkage_avoids_chaining():
    # Same chaining setup as above, but complete-linkage requires ALL
    # pairwise distances within a cluster to be below threshold, so a-c's
    # large distance should keep {a,b} and {c} split even though a-b and
    # b-c are both close.
    pairs = [("a", "b", 0.0005), ("b", "c", 0.0005), ("a", "c", 0.01)]
    groups = hierarchical_groups(pairs, threshold=0.001, nodes={"a", "b", "c"}, method="complete")
    assert {"a", "b", "c"} not in groups
    assert any(g == {"a", "b"} for g in groups) or all(len(g) == 1 for g in groups)


def test_groups_to_labels_assigns_same_label_within_group():
    groups = [{"a", "b"}, {"c"}]
    labels = groups_to_labels(groups)
    assert labels["a"] == labels["b"]
    assert labels["c"] != labels["a"]


def test_ari_against_ground_truth_perfect_agreement():
    pred = {"a": 0, "b": 0, "c": 1, "d": 1}
    truth = {"a": "X", "b": "X", "c": "Y", "d": "Y"}
    ari = ari_against_ground_truth(pred, truth)
    assert ari == 1.0


def test_ari_against_ground_truth_restricts_to_common_strains():
    # "e" has no ground truth label and must be excluded, not crash.
    pred = {"a": 0, "b": 0, "c": 1, "e": 5}
    truth = {"a": "X", "b": "X", "c": "Y"}
    ari = ari_against_ground_truth(pred, truth)
    assert ari == 1.0


def test_ari_against_ground_truth_no_overlap_returns_none():
    pred = {"a": 0}
    truth = {"z": "X"}
    assert ari_against_ground_truth(pred, truth) is None
