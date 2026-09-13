import sys
from pathlib import Path

import pandas as pd
import numpy as np

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


def test_high_confidence_starships_handles_nan_as_false():
    """NaN in presence/absence must not silently become True."""
    table_s6 = pd.DataFrame({"nameID": ["Gnosis-h1"], "freq": [0.599]})
    table_s21 = pd.DataFrame({
        "isolateID": ["iso1", "iso2", "iso3"],
        "nameID": ["Gnosis-h1", "Gnosis-h1", "Gnosis-h1"],
        "presence/absence": [1, 0, np.nan],  # iso3 has missing data
    })
    result = high_confidence_starships(table_s6, table_s21)
    # NaN should become False, not True
    assert result["Gnosis-h1"]["presence"] == {"iso1": True, "iso2": False, "iso3": False}


def test_cargo_gene_sets_deduplicates_duplicate_pairs():
    """Confirm set deduplication: same (starshipID, geneID) pair appearing multiple times."""
    table_s13 = pd.DataFrame({
        "starshipID": ["47-10_s00011", "47-10_s00011", "47-10_s00011", "47-10_s00011"],
        "geneID": ["47-10_000766", "47-10_000767", "47-10_000766", "47-10_000768"],
    })
    result = cargo_gene_sets(table_s13)
    # 4 input rows, but 47-10_000766 appears twice → 3 unique genes in the set
    assert result["47-10_s00011"] == {"47-10_000766", "47-10_000767", "47-10_000768"}
    assert len(result["47-10_s00011"]) == 3  # Confirm dedup by checking set size
