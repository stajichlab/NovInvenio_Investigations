import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from dereplicate_strains import (
    compute_assembly_stats, parse_mash_dist, choose_representatives,
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
