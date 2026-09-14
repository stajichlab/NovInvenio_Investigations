import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from synteny_windows import parse_gff3_gene_order, sliding_windows, accessory_islands, linkage_fraction


def test_parse_gff3_gene_order_sorts_by_contig_and_start(tmp_path):
    gff3 = tmp_path / "strain.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "contig1\tsrc\tgene\t500\t600\t.\t+\t.\tID=geneB\n"
        "contig1\tsrc\tgene\t100\t200\t.\t+\t.\tID=geneA\n"
        "contig2\tsrc\tgene\t1\t50\t.\t+\t.\tID=geneC\n"
        "contig1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=geneA.t1;Parent=geneA\n"
    )
    order = parse_gff3_gene_order(gff3)
    assert order == [
        ("geneA", "contig1", 100, 200),
        ("geneB", "contig1", 500, 600),
        ("geneC", "contig2", 1, 50),
    ]


def test_sliding_windows_excludes_cross_contig():
    gene_order = [
        ("g1", "c1", 1, 10), ("g2", "c1", 20, 30), ("g3", "c1", 40, 50),
        ("g4", "c2", 1, 10), ("g5", "c2", 20, 30),
    ]
    windows = sliding_windows(gene_order, n=3)
    gene_id_windows = [[g[0] for g in w] for w in windows]
    assert ["g1", "g2", "g3"] in gene_id_windows
    # no window of size 3 crosses from c1 into c2
    assert not any("g4" in w and "g3" in w for w in gene_id_windows)


def test_accessory_islands_finds_maximal_noncore_runs():
    gene_order = [
        ("g1", "c1", 1, 10), ("g2", "c1", 20, 30), ("g3", "c1", 40, 50),
        ("g4", "c1", 60, 70), ("g5", "c1", 80, 90),
    ]
    is_core = {"g1": True, "g2": False, "g3": False, "g4": False, "g5": True}
    islands = accessory_islands(gene_order, is_core)
    assert [g[0] for g in islands[0]] == ["g2", "g3", "g4"]


def test_accessory_islands_never_merges_across_contig_boundary():
    # contig1: two consecutive non-core genes (an island)
    # contig2: one non-core gene followed by a core gene (a one-gene island)
    # contig3: one non-core gene, still "open" at the end of the list
    gene_order = [
        ("g1", "c1", 1, 10), ("g2", "c1", 20, 30),
        ("g3", "c2", 1, 10), ("g4", "c2", 20, 30),
        ("g5", "c3", 1, 10),
    ]
    is_core = {"g1": False, "g2": False, "g3": False, "g4": True, "g5": False}
    islands = accessory_islands(gene_order, is_core)
    gene_id_islands = [[g[0] for g in island] for island in islands]
    assert gene_id_islands == [["g1", "g2"], ["g3"], ["g5"]]
    # In particular: g2 (end of c1) and g3 (start of c2) must never be
    # merged into one island, even though both are non-core and adjacent
    # in gene_order.
    assert not any("g2" in island and "g3" in island for island in gene_id_islands)


def test_parse_gff3_gene_order_anchors_id_attribute(tmp_path):
    # A decoy attribute ending in "ID=" before the real "ID=" must not be
    # picked up by an unanchored regex search.
    gff3 = tmp_path / "decoy.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "contig1\tsrc\tgene\t100\t200\t.\t+\t.\torig_protein_ID=XP_1;ID=geneA\n"
    )
    order = parse_gff3_gene_order(gff3)
    assert order == [("geneA", "contig1", 100, 200)]


def test_linkage_fraction_measures_physical_proximity():
    gene_position = {
        "s1": {"famA": [("c1", 5)], "famB": [("c1", 7)]},   # 2 genes apart, within k
        "s2": {"famA": [("c1", 5)], "famB": [("c1", 500)]},  # far apart, same contig
        "s3": {"famA": [("c1", 5)]},                         # famB absent in s3
    }
    frac = linkage_fraction("famA", "famB", gene_position, k=10)
    assert frac == 0.5  # only s1 of {s1, s2} (both-present strains) is within k


def test_linkage_fraction_multi_copy_family_checks_all_copy_pairs():
    # famA has 2 copies in s1: one far from famB, one close -- must count as
    # linked because SOME copy pair is within k, not just the first stored.
    gene_position = {
        "s1": {"famA": [("c1", 5), ("c1", 900)], "famB": [("c1", 7)]},
    }
    frac = linkage_fraction("famA", "famB", gene_position, k=10)
    assert frac == 1.0


def test_linkage_fraction_multi_copy_family_still_zero_if_no_pair_is_close():
    gene_position = {
        "s1": {"famA": [("c1", 5), ("c1", 900)], "famB": [("c1", 500)]},
    }
    frac = linkage_fraction("famA", "famB", gene_position, k=10)
    assert frac == 0.0
