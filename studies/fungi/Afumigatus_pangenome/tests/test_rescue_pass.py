import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT
from rescue_pass import parse_tblastn_hits, apply_rescue


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
