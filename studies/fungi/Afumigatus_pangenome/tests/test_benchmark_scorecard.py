import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from benchmark_scorecard import score_presence_recovery, score_cargo_grouping


def test_score_presence_recovery_computes_confusion_counts():
    predicted = {"iso1": True, "iso2": False, "iso3": True, "iso4": False}
    truth =     {"iso1": True, "iso2": True,  "iso3": True, "iso4": False}
    result = score_presence_recovery(predicted, truth)
    assert result["tp"] == 2
    assert result["fn"] == 1
    assert result["fp"] == 0
    assert result["tn"] == 1
    assert result["accuracy"] == 0.75


def test_score_cargo_grouping_purity_and_completeness():
    # our clustering put gene1, gene2 in family "famX", and gene3 alone in "famY"
    predicted_family_of_gene = {"gene1": "famX", "gene2": "famX", "gene3": "famY"}
    # ground truth: gene1, gene2, gene3 all belong to the same Starship's cargo
    truth_cargo_sets = {"starship1": {"gene1", "gene2", "gene3"}}
    result = score_cargo_grouping(predicted_family_of_gene, truth_cargo_sets)
    # famX correctly grouped 2 of the 3 true cargo genes together (some completeness),
    # and every gene it contains is a true cargo member (perfect purity)
    assert result["starship1"]["purity"] == 1.0
    assert round(result["starship1"]["completeness"], 4) == round(2 / 3, 4)
