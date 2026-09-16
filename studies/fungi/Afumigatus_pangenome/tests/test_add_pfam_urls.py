import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from add_pfam_urls import load_name_to_accession, add_urls, PFAM_URL_TEMPLATE


def _domtblout_line(name, accession):
    cols = ["-"] * 22
    cols[0] = name
    cols[1] = accession
    return " ".join(cols)


def test_load_name_to_accession_strips_version_suffix(tmp_path):
    path = tmp_path / "hits.domtblout"
    path.write_text(_domtblout_line("DUF3435", "PF11917.14") + "\n# comment\n")
    mapping = load_name_to_accession([str(path)])
    assert mapping == {"DUF3435": "PF11917"}


def test_load_name_to_accession_merges_multiple_files(tmp_path):
    path1 = tmp_path / "a.domtblout"
    path1.write_text(_domtblout_line("DomainA", "PF00001.1") + "\n")
    path2 = tmp_path / "b.domtblout"
    path2.write_text(_domtblout_line("DomainB", "PF00002.1") + "\n")
    mapping = load_name_to_accession([str(path1), str(path2)])
    assert mapping == {"DomainA": "PF00001", "DomainB": "PF00002"}


def test_add_urls_builds_canonical_interpro_pfam_link():
    rows = [{"domain": "DUF3435", "fdr_q": "1e-10"}]
    annotated = add_urls(rows, {"DUF3435": "PF11917"})
    assert annotated[0]["pfam_accession"] == "PF11917"
    assert annotated[0]["pfam_url"] == "https://www.ebi.ac.uk/interpro/entry/pfam/PF11917/"


def test_add_urls_leaves_url_blank_when_accession_unknown():
    rows = [{"domain": "UnknownDomain", "fdr_q": "1e-2"}]
    annotated = add_urls(rows, {})
    assert annotated[0]["pfam_accession"] == ""
    assert annotated[0]["pfam_url"] == ""


def test_pfam_url_template_has_no_version_suffix():
    assert PFAM_URL_TEMPLATE.format(accession="PF00001") == "https://www.ebi.ac.uk/interpro/entry/pfam/PF00001/"
