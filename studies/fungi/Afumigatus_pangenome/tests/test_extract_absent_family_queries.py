import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY
from extract_absent_family_queries import (
    read_fasta_records, absent_families_by_strain,
)


def test_read_fasta_records_keys_on_first_whitespace_token():
    fasta_text = ">famA some description here\nMKV\nLLP\n>famB\nQQQ\n"
    tmp = Path("/tmp/test_extract_absent_family_queries_fasta.fa")
    tmp.write_text(fasta_text)
    try:
        records = read_fasta_records(str(tmp))
        assert set(records) == {"famA", "famB"}
        assert "".join(records["famA"]) == ">famA some description here\nMKV\nLLP\n"
        assert "".join(records["famB"]) == ">famB\nQQQ\n"
    finally:
        tmp.unlink()


def test_absent_families_by_strain_only_includes_absent_calls():
    strains = ["s1", "s2"]
    pm = PresenceMatrix(families=["famA", "famB", "famC"], strains=strains)
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s1", GENOME_ONLY)
    # famC left at its default (absent) for both strains; famA/famB default
    # absent for s2 too (never called present/genome_only for s2)
    by_strain = absent_families_by_strain(pm)
    assert by_strain["s1"] == ["famC"]
    assert by_strain["s2"] == ["famA", "famB", "famC"]


def test_main_writes_one_fasta_per_strain_and_a_manifest(tmp_path):
    from extract_absent_family_queries import main
    import argparse

    strains = ["s1", "s2"]
    pm = PresenceMatrix(families=["famA", "famB"], strains=strains)
    pm.set_call("famA", "s1", PRESENT)
    # famA absent in s2, famB absent in both
    matrix_path = tmp_path / "matrix.tsv"
    pm.to_tsv(matrix_path)

    rep_fasta = tmp_path / "reps.fa"
    rep_fasta.write_text(">famA desc\nMKV\n>famB desc\nQQQ\n")

    out_dir = tmp_path / "per_strain_queries"

    import sys as _sys
    old_argv = _sys.argv
    _sys.argv = [
        "extract_absent_family_queries.py",
        "--matrix", str(matrix_path),
        "--rep_fasta", str(rep_fasta),
        "--out_dir", str(out_dir),
    ]
    try:
        main()
    finally:
        _sys.argv = old_argv

    # s1: only famB absent
    s1_fa = (out_dir / "s1.absent.fa").read_text()
    assert "famA" not in s1_fa and "famB" in s1_fa
    # s2: both absent
    s2_fa = (out_dir / "s2.absent.fa").read_text()
    assert "famA" in s2_fa and "famB" in s2_fa

    manifest = (out_dir / "manifest.tsv").read_text().splitlines()
    assert manifest[0] == "strain\tn_absent_families\tquery_fasta"
    rows = {line.split("\t")[0]: line.split("\t")[1] for line in manifest[1:]}
    assert rows["s1"] == "1"
    assert rows["s2"] == "2"


def test_main_warns_on_missing_representative_sequence(tmp_path, capsys):
    from extract_absent_family_queries import main
    import sys as _sys

    strains = ["s1"]
    pm = PresenceMatrix(families=["famA", "famMissing"], strains=strains)
    matrix_path = tmp_path / "matrix.tsv"
    pm.to_tsv(matrix_path)

    rep_fasta = tmp_path / "reps.fa"
    rep_fasta.write_text(">famA desc\nMKV\n")  # famMissing has no rep sequence

    out_dir = tmp_path / "per_strain_queries"
    old_argv = _sys.argv
    _sys.argv = [
        "extract_absent_family_queries.py",
        "--matrix", str(matrix_path),
        "--rep_fasta", str(rep_fasta),
        "--out_dir", str(out_dir),
    ]
    try:
        main()
    finally:
        _sys.argv = old_argv

    err = capsys.readouterr().err
    assert "WARNING" in err and "1 family IDs" in err
