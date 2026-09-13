import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY
from hac_screen import screen_family


def test_screen_family_reports_state_copy_number_and_starship():
    pm = PresenceMatrix(families=["hacFamily"], strains=["s1", "s2", "s3"])
    pm.set_call("hacFamily", "s1", PRESENT, copies=2)
    pm.set_call("hacFamily", "s2", GENOME_ONLY)
    # s3 left ABSENT (default)
    strain_to_starship = {"s1": "Nebuchadnezzar-h1"}

    rows = screen_family(pm, "hacFamily", strain_to_starship)
    by_strain = {r["strain"]: r for r in rows}

    assert by_strain["s1"]["state"] == PRESENT
    assert by_strain["s1"]["copy_number"] == 2
    assert by_strain["s1"]["starship"] == "Nebuchadnezzar-h1"
    assert by_strain["s2"]["state"] == GENOME_ONLY
    assert by_strain["s2"]["starship"] == "unknown"
    assert by_strain["s3"]["copy_number"] == 0
