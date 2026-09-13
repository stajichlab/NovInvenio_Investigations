import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from cluster_backend import two_tier_families


def test_two_tier_families_attaches_superfamily_label(tmp_path):
    tier1 = tmp_path / "tier1_cluster.tsv"
    # two tier-1 families: repA (with member memberA2) and repB (singleton)
    tier1.write_text("repA\trepA\nrepA\tmemberA2\nrepB\trepB\n")
    tier2 = tmp_path / "tier2_cluster.tsv"
    # tier-2 re-clusters the tier-1 reps: repA and repB fall under superfamily repA
    tier2.write_text("repA\trepA\nrepA\trepB\n")

    result = two_tier_families(tier1, tier2)

    assert result["repA"]["members"] == ["memberA2", "repA"]
    assert result["repA"]["superfamily"] == "repA"
    assert result["repB"]["members"] == ["repB"]
    assert result["repB"]["superfamily"] == "repA"
