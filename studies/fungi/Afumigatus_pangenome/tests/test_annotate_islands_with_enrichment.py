import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from annotate_islands_with_enrichment import annotate_islands


def test_annotate_islands_flags_significant_domains_sorted_by_q():
    island_rows = [
        {"island_size": "5", "pfam_domains": "DomainA,DomainB,DomainC"},
    ]
    domain_qvalues = {"DomainA": 0.20, "DomainB": 0.001, "DomainC": 0.03}
    annotated = annotate_islands(island_rows, domain_qvalues, fdr_alpha=0.05)
    assert annotated[0]["n_significant_domains"] == "2"
    assert annotated[0]["significant_domains"] == "DomainB:1.00e-03,DomainC:3.00e-02"
    assert annotated[0]["min_fdr_q"] == "1.00e-03"


def test_annotate_islands_handles_no_domains():
    island_rows = [{"island_size": "2", "pfam_domains": "-"}]
    annotated = annotate_islands(island_rows, {"DomainA": 0.01}, fdr_alpha=0.05)
    assert annotated[0]["n_significant_domains"] == "0"
    assert annotated[0]["significant_domains"] == ""
    assert annotated[0]["min_fdr_q"] == ""


def test_annotate_islands_ignores_domains_not_in_enrichment_table():
    # A domain that appears on the island but was never tested (e.g. it
    # only occurs in families outside the eligible background) should be
    # silently skipped, not raise a KeyError.
    island_rows = [{"island_size": "3", "pfam_domains": "DomainX,DomainUntested"}]
    annotated = annotate_islands(island_rows, {"DomainX": 0.001}, fdr_alpha=0.05)
    assert annotated[0]["significant_domains"] == "DomainX:1.00e-03"


def test_annotate_islands_preserves_input_row_order():
    island_rows = [
        {"island_size": "10", "pfam_domains": "-"},
        {"island_size": "5", "pfam_domains": "-"},
    ]
    annotated = annotate_islands(island_rows, {}, fdr_alpha=0.05)
    assert [r["island_size"] for r in annotated] == ["10", "5"]
