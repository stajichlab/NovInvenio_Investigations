import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from map_pfam_to_go import load_pfam_to_go, add_go_terms


def test_load_pfam_to_go_parses_real_format(tmp_path):
    path = tmp_path / "pfam2go.txt"
    path.write_text(
        "!version date: 2026/07/06\n"
        "!comment line, skipped\n"
        "Pfam:PF03184 DDE_1 > GO:nucleic acid binding ; GO:0003676\n"
        "Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930\n"
    )
    mapping = load_pfam_to_go(str(path))
    assert mapping["PF03184"] == [("GO:0003676", "nucleic acid binding")]
    assert mapping["PF00001"] == [("GO:0004930", "G protein-coupled receptor activity")]


def test_load_pfam_to_go_collects_multiple_terms_per_accession(tmp_path):
    path = tmp_path / "pfam2go.txt"
    path.write_text(
        "Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930\n"
        "Pfam:PF00001 7tm_1 > GO:membrane ; GO:0016020\n"
    )
    mapping = load_pfam_to_go(str(path))
    assert len(mapping["PF00001"]) == 2


def test_add_go_terms_joins_by_pfam_accession():
    rows = [{"domain": "DDE_1", "pfam_accession": "PF03184"}]
    pfam_to_go = {"PF03184": [("GO:0003676", "nucleic acid binding")]}
    annotated = add_go_terms(rows, pfam_to_go)
    assert annotated[0]["n_go_terms"] == "1"
    assert annotated[0]["go_terms"] == "GO:0003676 (nucleic acid binding)"


def test_add_go_terms_handles_domain_with_no_go_mapping():
    rows = [{"domain": "DUF3435", "pfam_accession": "PF11917"}]
    annotated = add_go_terms(rows, pfam_to_go={})
    assert annotated[0]["n_go_terms"] == "0"
    assert annotated[0]["go_terms"] == ""


def test_add_go_terms_handles_missing_pfam_accession_column():
    # A row with no pfam_accession key at all (e.g. an unresolved domain)
    # must not crash -- treated the same as "no mapping found".
    rows = [{"domain": "Unresolved"}]
    annotated = add_go_terms(rows, pfam_to_go={"PF00001": [("GO:0004930", "x")]})
    assert annotated[0]["n_go_terms"] == "0"
