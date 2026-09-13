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


def _run_hac_screen(args: list[str]):
    import subprocess
    script = Path(__file__).parent.parent / "bin" / "hac_screen.py"
    return subprocess.run(
        [sys.executable, str(script)] + args, capture_output=True, text=True
    )


def _fixture_matrix(tmp_path):
    pm = PresenceMatrix(families=["realFamily"], strains=["s1"])
    pm.set_call("realFamily", "s1", PRESENT, copies=2)
    matrix = tmp_path / "matrix.tsv"
    pm.to_tsv(matrix)
    smap = tmp_path / "starship.tsv"
    smap.write_text("strain\tstarship\ns1\tOsiris-h3\n")
    return matrix, smap


def test_main_errors_on_an_unknown_family_id(tmp_path):
    matrix, smap = _fixture_matrix(tmp_path)
    proc = _run_hac_screen([
        "--matrix", str(matrix), "--hac_family_id", "typo_family",
        "--haca_family_id", "realFamily", "--strain_starship_map", str(smap),
        "--output", str(tmp_path / "out.tsv"),
    ])
    assert proc.returncode != 0
    assert "not a family" in proc.stderr


def test_main_reports_copy_number_loaded_from_the_matrix_file(tmp_path):
    """Copy numbers must survive to_tsv/from_tsv -- this is the whole point of
    the copy-number sidecar; before it, this column was always 0 on real runs."""
    matrix, smap = _fixture_matrix(tmp_path)
    out = tmp_path / "out.tsv"
    proc = _run_hac_screen([
        "--matrix", str(matrix), "--hac_family_id", "realFamily",
        "--haca_family_id", "realFamily", "--strain_starship_map", str(smap),
        "--output", str(out),
    ])
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    assert rows[0][:5] == ["HAC", "s1", PRESENT, "2", "Osiris-h3"]
