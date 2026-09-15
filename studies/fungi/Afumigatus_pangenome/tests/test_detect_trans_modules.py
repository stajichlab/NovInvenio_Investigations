import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from detect_trans_modules import load_trans_edges, build_graph, detect_modules


def test_load_trans_edges_filters_to_trans_classification_only(tmp_path):
    f = tmp_path / "pc.tsv"
    f.write_text(
        "family_a\tfamily_b\tclassification\tlinkage_fraction\tjaccard\t"
        "fisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n"
        "famA\tfamB\ttrans\t0.0\t0.5\t1e-4\t1e-3\t0.01\tgain\t{}\n"
        "famC\tfamD\tstarship_explained\t1.0\t0.5\t1e-4\t1e-3\t0.01\tgain\t{}\n"
        "famE\tfamF\ttrans\t0.0\t0.5\t1e-4\t1e-2\t0.01\tgain\t{}\n"
    )
    edges = load_trans_edges(str(f))
    pairs = {(a, b) for a, b, _ in edges}
    assert pairs == {("famA", "famB"), ("famE", "famF")}


def test_load_trans_edges_weight_is_higher_for_smaller_fdr_q(tmp_path):
    f = tmp_path / "pc.tsv"
    f.write_text(
        "family_a\tfamily_b\tclassification\tlinkage_fraction\tjaccard\t"
        "fisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n"
        "famA\tfamB\ttrans\t0.0\t0.5\t1e-10\t1e-9\t0.01\tgain\t{}\n"
        "famC\tfamD\ttrans\t0.0\t0.5\t1e-2\t1e-1\t0.01\tgain\t{}\n"
    )
    edges = load_trans_edges(str(f))
    weight_by_pair = {(a, b): w for a, b, w in edges}
    assert weight_by_pair[("famA", "famB")] > weight_by_pair[("famC", "famD")]


def test_load_trans_edges_clips_weight_for_zero_fdr_q(tmp_path):
    f = tmp_path / "pc.tsv"
    f.write_text(
        "family_a\tfamily_b\tclassification\tlinkage_fraction\tjaccard\t"
        "fisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n"
        "famA\tfamB\ttrans\t0.0\t0.5\t0.0\t0.0\t0.01\tgain\t{}\n"
    )
    edges = load_trans_edges(str(f))
    assert edges[0][2] == 300.0  # MAX_WEIGHT, not inf/nan


def test_build_graph_creates_correct_nodes_and_edges():
    edges = [("famA", "famB", 1.0), ("famB", "famC", 2.0)]
    g = build_graph(edges)
    assert g.vcount() == 3
    assert g.ecount() == 2
    assert set(g.vs["name"]) == {"famA", "famB", "famC"}


def test_detect_modules_recovers_two_disjoint_dense_cliques():
    # Two 5-node cliques with NO edges between them -- Leiden must put them
    # in separate modules.
    clique_a = [f"a{i}" for i in range(5)]
    clique_b = [f"b{i}" for i in range(5)]
    edges = []
    for i in range(5):
        for j in range(i + 1, 5):
            edges.append((clique_a[i], clique_a[j], 10.0))
            edges.append((clique_b[i], clique_b[j], 10.0))
    g = build_graph(edges)
    family_module = detect_modules(g, resolution=1.0, seed=0)

    modules_a = {family_module[f] for f in clique_a}
    modules_b = {family_module[f] for f in clique_b}
    assert len(modules_a) == 1  # all of clique_a in the same module
    assert len(modules_b) == 1  # all of clique_b in the same module
    assert modules_a != modules_b  # the two cliques are DIFFERENT modules


def test_detect_modules_renumbers_by_descending_size():
    # A big clique (6 nodes) plus a small disconnected pair (2 nodes) --
    # module 0 must be the big one.
    big = [f"big{i}" for i in range(6)]
    edges = [(big[i], big[j], 5.0) for i in range(6) for j in range(i + 1, 6)]
    edges.append(("small0", "small1", 5.0))
    g = build_graph(edges)
    family_module = detect_modules(g, resolution=1.0, seed=0)

    big_modules = {family_module[f] for f in big}
    assert big_modules == {0}
    assert family_module["small0"] == family_module["small1"]
    assert family_module["small0"] != 0
