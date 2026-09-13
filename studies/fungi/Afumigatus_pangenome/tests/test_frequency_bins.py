import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT
from frequency_bins import assign_bin, compute_frequency_table


def test_assign_bin_boundaries():
    assert assign_bin(freq=1.0, strain_count=100) == "core"
    assert assign_bin(freq=0.95, strain_count=95) == "core"
    assert assign_bin(freq=0.92, strain_count=92) == "soft_core"
    assert assign_bin(freq=0.50, strain_count=50) == "shell"
    assert assign_bin(freq=0.05, strain_count=5) == "cloud"
    assert assign_bin(freq=0.01, strain_count=1) == "singleton"


def test_compute_frequency_table():
    pm = PresenceMatrix(families=["famCore", "famSingleton"], strains=["s1", "s2", "s3", "s4"])
    for s in pm.strains:
        pm.set_call("famCore", s, PRESENT)
    pm.set_call("famSingleton", "s1", PRESENT)

    table = compute_frequency_table(pm)
    by_family = {row["family"]: row for row in table}
    assert by_family["famCore"]["bin"] == "core"
    assert by_family["famCore"]["strain_count"] == 4
    assert by_family["famSingleton"]["bin"] == "singleton"
    assert by_family["famSingleton"]["strain_count"] == 1
