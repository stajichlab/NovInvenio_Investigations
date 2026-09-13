import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from dereplicate_strains import (
    compute_assembly_stats, parse_mash_dist, choose_representatives,
    groups_to_shorts,
)


def test_compute_assembly_stats_n50_and_contig_count(tmp_path):
    fa = tmp_path / "strain.dna.fa"
    # two contigs: 100bp and 300bp -> total 400, N50 is the 300bp contig
    fa.write_text(">contig1\n" + "A" * 100 + "\n>contig2\n" + "C" * 300 + "\n")
    stats = compute_assembly_stats(fa)
    assert stats["n_contigs"] == 2
    assert stats["total_length"] == 400
    assert stats["n50"] == 300


def test_parse_mash_dist_groups_below_threshold():
    # mash dist -t output: query, ref1_dist, ref2_dist, ... one row per query
    # here: s1~s2 near-identical (0.0005), s3 distinct (0.05)
    lines = [
        "#query\ts1\ts2\ts3",
        "s1\t0.0000\t0.0005\t0.0500",
        "s2\t0.0005\t0.0000\t0.0520",
        "s3\t0.0500\t0.0520\t0.0000",
    ]
    groups = parse_mash_dist(lines, threshold=0.001)
    assert {"s1", "s2"} in groups
    assert {"s3"} in groups


def test_choose_representatives_picks_highest_n50():
    groups = [{"s1", "s2"}, {"s3"}]
    assembly_stats = {
        "s1": {"n50": 500_000, "n_contigs": 20, "total_length": 29_000_000},
        "s2": {"n50": 2_000_000, "n_contigs": 8, "total_length": 29_100_000},
        "s3": {"n50": 1_000_000, "n_contigs": 15, "total_length": 28_900_000},
    }
    reps = choose_representatives(groups, assembly_stats)
    assert reps == {"s1": "s2", "s2": "s2", "s3": "s3"}


def test_parse_mash_dist_empty_input():
    # Empty input should return empty groups
    groups = parse_mash_dist([], threshold=0.001)
    assert groups == []


def test_choose_representatives_deterministic_tiebreak_on_equal_n50():
    # When N50 values are equal, lexically later Short should be chosen
    # (deterministic, independent of set iteration order)
    groups = [{"s1", "s2"}]
    assembly_stats = {
        "s1": {"n50": 2_000_000, "n_contigs": 10, "total_length": 29_000_000},
        "s2": {"n50": 2_000_000, "n_contigs": 10, "total_length": 29_000_000},
    }
    reps = choose_representatives(groups, assembly_stats)
    # With deterministic tie-break on (n50, short), s2 > s1, so s2 is rep
    assert reps == {"s1": "s2", "s2": "s2"}


def test_parse_mash_dist_transitivity():
    # Test union-find transitivity: A~B, B~C should put all in same group
    # even though A and C are not directly below threshold
    lines = [
        "#query\tA\tB\tC",
        "A\t0.0000\t0.0005\t0.0100",  # A close to B, far from C
        "B\t0.0005\t0.0000\t0.0005",  # B close to both A and C
        "C\t0.0100\t0.0005\t0.0000",  # C close to B, far from A
    ]
    groups = parse_mash_dist(lines, threshold=0.001)
    # With transitivity: A~B and B~C means all in same group
    assert len(groups) == 1
    assert {"A", "B", "C"} in groups


def test_groups_to_shorts_translates_path_strings_to_short_ids():
    # Simulate mash output: full paths in groups, map them to Short IDs
    # This tests the critical translation that happens in main()
    path_groups = [
        {"/path/to/data/dna/s1.fa", "/path/to/data/dna/s2.fa"},
        {"/path/to/data/dna/s3.fa"},
    ]
    path_to_short = {
        "/path/to/data/dna/s1.fa": "s1",
        "/path/to/data/dna/s2.fa": "s2",
        "/path/to/data/dna/s3.fa": "s3",
    }
    short_groups = groups_to_shorts(path_groups, path_to_short)
    # Should convert path strings back to Short IDs
    assert {"s1", "s2"} in short_groups
    assert {"s3"} in short_groups
    assert len(short_groups) == 2
