import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from build_family_positions import build_family_positions


def test_build_family_positions_single_copy_per_strain():
    rows = [
        ("s1", "P1.1", "contig1", 100, 200),
        ("s1", "P2.1", "contig1", 300, 400),
    ]
    member_to_rep = {"s1|P1.1": "s1|P1.1", "s1|P2.1": "s1|P1.1"}  # same family
    result = build_family_positions(rows, member_to_rep)
    assert result == {"s1": {"s1|P1.1": [("contig1", 0), ("contig1", 1)]}}


def test_build_family_positions_ranks_are_contig_relative_and_start_sorted():
    rows = [
        ("s1", "P3.1", "contig1", 500, 600),  # later start -> higher rank
        ("s1", "P1.1", "contig1", 100, 200),
        ("s1", "P2.1", "contig2", 50, 80),    # different contig
    ]
    member_to_rep = {
        "s1|P1.1": "famA", "s1|P2.1": "famB", "s1|P3.1": "famC",
    }
    result = build_family_positions(rows, member_to_rep)
    # sorted by (contig, start): P1.1(contig1,100)=rank0, P3.1(contig1,500)=rank1, P2.1(contig2,50)=rank2
    assert result["s1"]["famA"] == [("contig1", 0)]
    assert result["s1"]["famC"] == [("contig1", 1)]
    assert result["s1"]["famB"] == [("contig2", 2)]


def test_build_family_positions_multi_copy_family_gets_list_of_positions():
    rows = [
        ("s1", "P1.1", "contig1", 100, 200),
        ("s1", "P2.1", "contig1", 900, 950),
    ]
    # Both proteins belong to the SAME family (a duplicated/multi-copy family)
    member_to_rep = {"s1|P1.1": "famHAC", "s1|P2.1": "famHAC"}
    result = build_family_positions(rows, member_to_rep)
    assert result["s1"]["famHAC"] == [("contig1", 0), ("contig1", 1)]


def test_build_family_positions_skips_proteins_with_no_cluster_membership():
    rows = [("s1", "P1.1", "contig1", 100, 200)]
    member_to_rep = {}  # P1.1 was never clustered
    result = build_family_positions(rows, member_to_rep)
    assert result == {"s1": {}}


def test_build_family_positions_separates_strains():
    rows = [
        ("s1", "P1.1", "contig1", 100, 200),
        ("s2", "P1.1", "contig1", 100, 200),  # same protein_id, different strain
    ]
    member_to_rep = {"s1|P1.1": "famA", "s2|P1.1": "famA"}
    result = build_family_positions(rows, member_to_rep)
    assert set(result.keys()) == {"s1", "s2"}


def test_build_family_positions_merges_rescue_positions_interleaved_by_coordinate():
    # 2026-09-15 fix: a rescue-pass (GENOME_ONLY) hit has no protein_id, so
    # its family is already resolved directly -- it must still be
    # interleaved into the SAME per-strain rank ordering by genomic
    # coordinate, not appended after all annotated genes, or the rank
    # wouldn't reflect true physical distance for linkage_fraction.
    gene_rows = [
        ("s1", "P1.1", "contig1", 100, 200),   # rank 0
        ("s1", "P2.1", "contig1", 500, 600),   # rank 2 (rescue entry sits between)
    ]
    member_to_rep = {"s1|P1.1": "famA", "s1|P2.1": "famB"}
    rescue_rows = [("s1", "famRescued", "contig1", 300)]  # rank 1, between famA and famB

    result = build_family_positions(gene_rows, member_to_rep, rescue_position_rows=rescue_rows)
    assert result["s1"]["famA"] == [("contig1", 0)]
    assert result["s1"]["famRescued"] == [("contig1", 1)]
    assert result["s1"]["famB"] == [("contig1", 2)]


def test_build_family_positions_rescue_positions_default_to_none_unchanged():
    # Backward compatibility: omitting rescue_position_rows entirely must
    # behave exactly as before (existing callers/tests pass no such arg).
    rows = [("s1", "P1.1", "contig1", 100, 200)]
    member_to_rep = {"s1|P1.1": "famA"}
    assert build_family_positions(rows, member_to_rep) == {"s1": {"famA": [("contig1", 0)]}}
