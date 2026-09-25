import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from filter_config_by_taxon import filter_rows


def test_filter_rows_by_taxon_group():
    rows = [
        {"TaxonGroup": "Coccidioides immitis", "Short": "1M0"},
        {"TaxonGroup": "Coccidioides posadasii", "Short": "B3224"},
    ]
    result = filter_rows(rows, taxon_group="Coccidioides immitis")
    assert [r["Short"] for r in result] == ["1M0"]
