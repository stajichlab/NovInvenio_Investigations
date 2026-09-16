import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from cluster_homology_test import (
    parse_blast_outfmt6,
    best_hit_per_query,
    summarize_cluster_homology,
)

OUTFMT6_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore qlen slen"
)


def _row(qseqid, sseqid, pident, length, evalue, qlen, slen):
    return "\t".join([
        qseqid, sseqid, str(pident), str(length), "0", "0", "1", str(length),
        "1", str(length), str(evalue), "100", str(qlen), str(slen),
    ])


def test_parse_blast_outfmt6_parses_numeric_fields():
    lines = [_row("q1", "s1", 45.6, 158, 3.6e-44, 482, 189)]
    hits = parse_blast_outfmt6(lines)
    assert hits[0]["qseqid"] == "q1"
    assert hits[0]["pident"] == 45.6
    assert hits[0]["evalue"] == 3.6e-44
    assert hits[0]["qlen"] == 482
    assert hits[0]["qcov"] == 158 / 482


def test_best_hit_per_query_picks_lowest_evalue():
    hits = parse_blast_outfmt6([
        _row("q1", "s1", 30.0, 50, 1e-5, 100, 100),
        _row("q1", "s2", 90.0, 90, 1e-40, 100, 100),
    ])
    best = best_hit_per_query(hits)
    assert best["q1"]["sseqid"] == "s2"


def test_summarize_cluster_homology_real_cluster_case():
    # A genuine cluster homolog: most query genes hit distinct reference
    # genes at good coverage/identity.
    best_hits = {
        "q1": {"sseqid": "ref_geneA", "pident": 60.0, "qcov": 0.8, "evalue": 1e-30},
        "q2": {"sseqid": "ref_geneB", "pident": 55.0, "qcov": 0.75, "evalue": 1e-25},
        "q3": {"sseqid": "ref_geneC", "pident": 50.0, "qcov": 0.6, "evalue": 1e-15},
    }
    summary = summarize_cluster_homology(
        ["q1", "q2", "q3", "q4"], best_hits, evalue_threshold=1e-5, qcov_threshold=0.5,
    )
    assert summary["n_total"] == 4
    assert summary["n_with_qualifying_hit"] == 3
    assert summary["n_distinct_targets"] == 3
    assert summary["verdict"] == "supported"


def test_summarize_cluster_homology_single_coincidental_hit_case():
    # The real Phomopsin-island shape: one partial (33% query coverage),
    # moderate-identity hit to one reference gene, nothing else -- the low
    # coverage itself rules it out (a shared promiscuous domain, not a
    # full-length homolog), so it doesn't even count as "qualifying".
    best_hits = {
        "q1": {"sseqid": "ref_geneC", "pident": 45.6, "qcov": 158 / 482, "evalue": 3.6e-44},
    }
    summary = summarize_cluster_homology(
        ["q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8"], best_hits,
        evalue_threshold=1e-5, qcov_threshold=0.5,
    )
    assert summary["n_with_qualifying_hit"] == 0
    assert summary["n_distinct_targets"] == 0
    assert summary["verdict"] == "not_supported"


def test_summarize_cluster_homology_filters_low_coverage_hits():
    # A hit that passes the e-value bar but covers too little of the query
    # (e.g. one shared promiscuous domain) should not count as qualifying.
    best_hits = {"q1": {"sseqid": "ref_geneC", "pident": 45.6, "qcov": 0.1, "evalue": 1e-20}}
    summary = summarize_cluster_homology(
        ["q1"], best_hits, evalue_threshold=1e-5, qcov_threshold=0.5,
    )
    assert summary["n_with_qualifying_hit"] == 0
    assert summary["verdict"] == "not_supported"


def test_summarize_cluster_homology_one_qualifying_hit_is_still_not_supported():
    # Even a hit that DOES pass both e-value and coverage bars is not
    # enough on its own -- a real cluster needs >=2 genes agreeing.
    best_hits = {"q1": {"sseqid": "ref_geneC", "pident": 45.6, "qcov": 0.9, "evalue": 1e-40}}
    summary = summarize_cluster_homology(
        ["q1", "q2", "q3"], best_hits, evalue_threshold=1e-5, qcov_threshold=0.5,
    )
    assert summary["n_with_qualifying_hit"] == 1
    assert summary["verdict"] == "not_supported"


def test_summarize_cluster_homology_no_hits_at_all():
    summary = summarize_cluster_homology(["q1", "q2"], {}, evalue_threshold=1e-5, qcov_threshold=0.5)
    assert summary["n_with_qualifying_hit"] == 0
    assert summary["verdict"] == "not_supported"
