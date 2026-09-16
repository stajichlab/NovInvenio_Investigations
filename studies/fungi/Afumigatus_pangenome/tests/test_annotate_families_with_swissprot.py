import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from annotate_families_with_swissprot import (
    parse_self_uniprot_id, parse_diamond_stitle, load_diamond_besthits, annotate_families,
)


def test_parse_self_uniprot_id_recognizes_reviewed_swissprot():
    result = parse_self_uniprot_id("Asfu_Af293|sp|Q4WKX2|FGND_ASPFU")
    assert result == ("sp", "Q4WKX2", "FGND_ASPFU")


def test_parse_self_uniprot_id_recognizes_unreviewed_trembl():
    result = parse_self_uniprot_id("Aslen_ref|tr|A0ABQ1AYS6|A0ABQ1AYS6_ASPLE")
    assert result == ("tr", "A0ABQ1AYS6", "A0ABQ1AYS6_ASPLE")


def test_parse_self_uniprot_id_returns_none_for_ordinary_ncbi_id():
    assert parse_self_uniprot_id("Asfu_HMR_AF_706|OXN03990.1") is None


def test_parse_diamond_stitle_splits_accession_name_and_description():
    stitle = "sp|D4AXL1|BLML_ARTBC Beta-lactamase-like protein ARB_00930 OS=Arthroderma benhamiae OX=663331 GN=ARB_00930 PE=1 SV=1"
    accession, name, description = parse_diamond_stitle(stitle)
    assert accession == "D4AXL1"
    assert name == "BLML_ARTBC"
    assert description == "Beta-lactamase-like protein ARB_00930"


def test_load_diamond_besthits_parses_real_column_layout(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(
        "famA\tsp|D4AXL1|BLML_ARTBC\tsp|D4AXL1|BLML_ARTBC Beta-lactamase-like protein OS=X OX=1\t47.5\t101\t1.83e-24\t99.4\n"
    )
    hits = load_diamond_besthits(str(path))
    assert hits["famA"]["accession"] == "D4AXL1"
    assert hits["famA"]["name"] == "BLML_ARTBC"
    assert hits["famA"]["source"] == "blast"
    assert hits["famA"]["pident"] == "47.5"


def test_annotate_families_prefers_self_annotation_over_blast_hit():
    # famA has BOTH a self-embedded UniProt ID and a (hypothetical) blast
    # hit -- the self-annotation (exact) must win over the inferred one.
    family_ids = ["Asfu_Af293|sp|Q4WKX2|FGND_ASPFU"]
    diamond_hits = {"Asfu_Af293|sp|Q4WKX2|FGND_ASPFU": {
        "source": "blast", "accession": "OTHER", "name": "OTHER_NAME",
        "description": "", "pident": "80.0", "evalue": "1e-10",
    }}
    annotated = annotate_families(family_ids, diamond_hits)
    assert annotated[family_ids[0]]["source"] == "self_sp"
    assert annotated[family_ids[0]]["accession"] == "Q4WKX2"


def test_annotate_families_falls_back_to_blast_hit_for_ncbi_style_id():
    family_ids = ["Asfu_HMR_AF_706|OXN03990.1"]
    diamond_hits = {"Asfu_HMR_AF_706|OXN03990.1": {
        "source": "blast", "accession": "D4AXL1", "name": "BLML_ARTBC",
        "description": "Beta-lactamase-like protein", "pident": "47.5", "evalue": "1.83e-24",
    }}
    annotated = annotate_families(family_ids, diamond_hits)
    assert annotated[family_ids[0]]["accession"] == "D4AXL1"


def test_annotate_families_leaves_unannotated_families_out():
    family_ids = ["Asfu_X|KAH000001.1"]
    annotated = annotate_families(family_ids, diamond_hits={})
    assert annotated == {}
