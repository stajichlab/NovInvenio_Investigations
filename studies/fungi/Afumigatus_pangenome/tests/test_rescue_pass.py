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
    apply_rescue(pm, {("famA", "s1"), ("famA", "s2")})
    assert pm.call("famA", "s1") == PRESENT       # unchanged, was already PRESENT
    assert pm.call("famA", "s2") == GENOME_ONLY   # upgraded from ABSENT
    assert pm.call("famA", "s3") == ABSENT         # no rescue hit, stays ABSENT
