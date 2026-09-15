import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from build_summary_report import (
    band_counts, band_table_markdown, busco_table_markdown,
    hac_table_markdown, pair_classification_table_markdown, build_report,
)


def test_band_counts_tallies_by_bin():
    freq_table = [
        {"family": "f1", "bin": "core"}, {"family": "f2", "bin": "core"},
        {"family": "f3", "bin": "cloud"},
    ]
    assert band_counts(freq_table) == {"core": 2, "cloud": 1}


def test_band_table_markdown_includes_total_row():
    md = band_table_markdown({"core": 6, "cloud": 4})
    assert "core" in md and "cloud" in md
    assert "**Total**" in md and "10" in md


def test_band_table_markdown_uses_frequency_bins_py_actual_soft_core_spelling():
    # Regression: frequency_bins.py's real output uses "soft_core"
    # (underscore), not "soft-core" (hyphen) -- an earlier version of this
    # table hardcoded the wrong spelling, which silently split soft_core
    # into its own untracked extra row instead of merging into the
    # expected-order table (caught by a dry run against real study data).
    md = band_table_markdown({"core": 10, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 1})
    lines = md.splitlines()
    soft_core_lines = [l for l in lines if "soft_core" in l]
    assert len(soft_core_lines) == 1, f"expected exactly one soft_core row, got {soft_core_lines}"
    assert "| soft_core | 5 |" in soft_core_lines[0]


def test_busco_table_markdown_summarizes_completeness_and_contig_mix():
    rows = [
        {"Complete_pct": "99.0", "Contigs": "534"},
        {"Complete_pct": "97.3", "Contigs": "12"},
        {"Complete_pct": "98.5", "Contigs": "200"},
    ]
    md = busco_table_markdown(rows)
    assert "3" in md  # strain count
    assert "97.3" in md and "99.0" in md  # min/max complete
    # 2 strains > 100 contigs, 1 strain <= 20
    assert "97.3 / 98.5 / 99.0" in md


def test_hac_table_markdown_computes_presence_fractions():
    rows = [
        {"hacA_present": "Y", "hrmA_present": "Y"},
        {"hacA_present": "Y", "hrmA_present": "N"},
        {"hacA_present": "Y", "hrmA_present": "N"},
    ]
    md = hac_table_markdown(rows)
    assert "3/3" in md  # hacA present in all
    assert "1/3" in md  # hrmA present in 1


def test_pair_classification_table_markdown_sorts_by_count_descending():
    rows = [
        {"classification": "trans"}, {"classification": "trans"},
        {"classification": "starship_explained"},
    ]
    md = pair_classification_table_markdown(rows)
    trans_line = [l for l in md.splitlines() if l.startswith("| trans ")][0]
    starship_line = [l for l in md.splitlines() if l.startswith("| starship_explained ")][0]
    assert md.index(trans_line) < md.index(starship_line)
    assert "**Total FDR-significant pairs**" in md and "3" in md


def test_build_report_omits_optional_sections_when_not_given():
    freq_table = [{"family": "f1", "bin": "core"}]
    report = build_report(freq_table, None, None, None, "test label")
    assert "test label" in report
    assert "Pangenome composition" in report
    assert "BUSCO" not in report
    assert "HAC" not in report
    assert "Co-occurring pair classification" not in report


def test_build_report_includes_all_sections_when_given():
    freq_table = [{"family": "f1", "bin": "core"}]
    busco_rows = [{"Complete_pct": "99.0", "Contigs": "50"}]
    hac_rows = [{"hacA_present": "Y", "hrmA_present": "N"}]
    pc_rows = [{"classification": "trans"}]
    report = build_report(freq_table, busco_rows, hac_rows, pc_rows, "full test")
    assert "BUSCO" in report
    assert "HAC" in report
    assert "Co-occurring pair classification" in report
