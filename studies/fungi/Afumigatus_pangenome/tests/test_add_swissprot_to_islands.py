import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from add_swissprot_to_islands import load_family_names, add_swissprot_names


def test_load_family_names_includes_description_when_present(tmp_path):
    path = tmp_path / "annot.tsv"
    path.write_text(
        "family\tsource\taccession\tname\tdescription\tpident\tevalue\n"
        "famA\tblast\tD4AXL1\tBLML_ARTBC\tBeta-lactamase-like protein\t47.5\t1e-24\n"
        "famB\tself_sp\tQ4WKX2\tFGND_ASPFU\t\t100.0\t0.0\n"
    )
    names = load_family_names(str(path))
    assert names["famA"] == "BLML_ARTBC (Beta-lactamase-like protein)"
    assert names["famB"] == "FGND_ASPFU"


def test_add_swissprot_names_joins_by_member_families():
    island_rows = [{"member_families": "famA,famB,famC"}]
    family_names = {"famA": "NAME_A (desc A)", "famC": "NAME_C"}
    annotated = add_swissprot_names(island_rows, family_names)
    assert annotated[0]["n_swissprot_hits"] == "2"
    assert annotated[0]["swissprot_names"] == "NAME_A (desc A); NAME_C"


def test_add_swissprot_names_handles_no_hits():
    island_rows = [{"member_families": "famX,famY"}]
    annotated = add_swissprot_names(island_rows, family_names={})
    assert annotated[0]["n_swissprot_hits"] == "0"
    assert annotated[0]["swissprot_names"] == "-"
