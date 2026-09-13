import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT
from rescue_pass import parse_tblastn_hits, apply_rescue, main


def test_parse_tblastn_hits_filters_by_identity_and_coverage():
    # outfmt6 with qcovs appended: qseqid sseqid pident length ... qcovs
    # family "famA" query, subject header encodes strain as "s2|contig1"
    lines = [
        "famA\ts2|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "famA\ts3|contig1\t60.0\t100\t0\t0\t1\t100\t500\t600\t1e-10\t80\t95",  # too low identity
        "famA\ts4|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t50",  # too low coverage
    ]
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0)
    assert hits == {("famA", "s2")}


def test_apply_rescue_only_upgrades_absent_calls():
    pm = PresenceMatrix(families=["famA"], strains=["s1", "s2", "s3"])
    pm.set_call("famA", "s1", PRESENT)
    # s2, s3 default to ABSENT
    applied, skipped = apply_rescue(pm, {("famA", "s1"), ("famA", "s2")})
    assert pm.call("famA", "s1") == PRESENT       # unchanged, was already PRESENT
    assert pm.call("famA", "s2") == GENOME_ONLY   # upgraded from ABSENT
    assert pm.call("famA", "s3") == ABSENT         # no rescue hit, stays ABSENT
    assert applied == 1                           # only famA/s2 was upgraded
    assert skipped == 0                           # all hits recognized


def test_parse_tblastn_hits_skips_blank_and_comment_lines():
    lines = [
        "famA\ts2|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "",                                         # blank line
        "# This is a comment line",                 # comment
        "famB\ts3|contig2\t92.0\t100\t0\t0\t1\t100\t500\t600\t1e-40\t200\t85",
    ]
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0)
    assert hits == {("famA", "s2"), ("famB", "s3")}  # only 2 valid lines parsed


def test_apply_rescue_skips_unrecognized_strain_or_family():
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", ABSENT)
    # Hits include: one valid (famA/s1), one with unrecognized strain (famA/s999),
    # one with unrecognized family (famZ/s1)
    applied, skipped = apply_rescue(pm, {("famA", "s1"), ("famA", "s999"), ("famZ", "s1")})
    assert pm.call("famA", "s1") == GENOME_ONLY   # only this one upgraded
    assert applied == 1
    assert skipped == 2                           # two unrecognized pairs


def test_apply_rescue_idempotence_with_genome_only():
    pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", GENOME_ONLY)       # pre-existing GENOME_ONLY
    # Hit on an already-GENOME_ONLY cell
    applied, skipped = apply_rescue(pm, {("famA", "s1")})
    assert pm.call("famA", "s1") == GENOME_ONLY   # unchanged (not ABSENT, so not upgraded)
    assert applied == 0                           # not upgraded since not ABSENT
    assert skipped == 0


def test_main_exits_when_all_hits_are_skipped(monkeypatch, capsys):
    """main() should hard-error if ALL hits are unrecognized (skipped)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create a matrix with strain s1, s2
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        # Create tblastn output with all hits on unrecognized strain s999
        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text("famA\ts999|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n")

        output_file = tmpdir / "output.tsv"

        # Mock sys.argv
        monkeypatch.setattr(sys, "argv", [
            "rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        # Should exit with code 1
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1

        # Check that error message was printed to stderr
        captured = capsys.readouterr()
        assert "All 1 parsed tblastn hits were skipped" in captured.err


def test_main_succeeds_when_all_hits_recognized_but_nothing_applied(monkeypatch, capsys):
    """main() should succeed (not hard-error) when all hits are recognized but
    nothing is applied (e.g., hits land on already-PRESENT/GENOME_ONLY cells)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create a matrix with famA/s1 already PRESENT
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", PRESENT)  # Already PRESENT, not ABSENT
        pm.set_call("famA", "s2", PRESENT)  # Also PRESENT
        pm.to_tsv(matrix_file)

        # Create tblastn output on recognized strains, but cells are already PRESENT
        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
            "famA\ts2|contig2\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
        )

        output_file = tmpdir / "output.tsv"

        # Mock sys.argv
        monkeypatch.setattr(sys, "argv", [
            "rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        # Should succeed (no SystemExit)
        main()

        # Check that output file was written
        assert output_file.exists()

        # Check that status message was printed
        captured = capsys.readouterr()
        assert "Rescue: 0 ABSENT→GENOME_ONLY, 0 skipped" in captured.err
        # Should NOT have the error message
        assert "All" not in captured.err or "were skipped" not in captured.err
