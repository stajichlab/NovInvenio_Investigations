import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from parse_starship_supplement import high_confidence_starships, cargo_gene_sets


def test_high_confidence_starships_combines_frequency_and_presence():
    table_s6 = pd.DataFrame({"nameID": ["Gnosis-h1", "Osiris-h4"], "freq": [0.599, 0.20]})
    table_s21 = pd.DataFrame({
        "isolateID": ["iso1", "iso2", "iso1", "iso2"],
        "nameID": ["Gnosis-h1", "Gnosis-h1", "Osiris-h4", "Osiris-h4"],
        "presence/absence": [1, 0, 0, 1],
    })
    result = high_confidence_starships(table_s6, table_s21)
    assert result["Gnosis-h1"]["population_freq"] == 0.599
    assert result["Gnosis-h1"]["presence"] == {"iso1": True, "iso2": False}
    assert result["Osiris-h4"]["presence"] == {"iso1": False, "iso2": True}


def test_cargo_gene_sets_groups_by_starship_id():
    table_s13 = pd.DataFrame({
        "starshipID": ["47-10_s00011", "47-10_s00011", "1F1SW-F4_e00006"],
        "geneID": ["47-10_000766", "47-10_000767", "1F1SW-F4_g001"],
    })
    result = cargo_gene_sets(table_s13)
    assert result["47-10_s00011"] == {"47-10_000766", "47-10_000767"}
    assert result["1F1SW-F4_e00006"] == {"1F1SW-F4_g001"}
