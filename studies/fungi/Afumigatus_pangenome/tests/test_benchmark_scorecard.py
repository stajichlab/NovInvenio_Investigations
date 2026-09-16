import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from benchmark_scorecard import (
    score_presence_recovery,
    score_cargo_grouping,
    load_crosswalk,
    load_cargo_by_name,
    load_starships_by_short,
    build_protein_to_family,
    score_backend,
    score_negative_control,
)
from pangenome_matrix import PresenceMatrix


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


def test_load_crosswalk_reads_paper_id_to_study_protein(tmp_path):
    path = tmp_path / "crosswalk.tsv"
    path.write_text("paper_id\tstudy_protein_id\nXP_1.1\tAsfu_Af293|Afu1g00100-T-p1\n")
    result = load_crosswalk(path)
    assert result == {"XP_1.1": "Asfu_Af293|Afu1g00100-T-p1"}


def test_load_cargo_by_name_groups_geneids_under_nameid(tmp_path):
    path = tmp_path / "cargo.tsv"
    path.write_text("nameID\tgeneID\nGnosis-h1\tgeneA\nGnosis-h1\tgeneB\nOsiris-h4\tgeneC\n")
    result = load_cargo_by_name(path)
    assert result == {"Gnosis-h1": {"geneA", "geneB"}, "Osiris-h4": {"geneC"}}


def test_load_starships_by_short_parses_presence_json(tmp_path):
    path = tmp_path / "starships.tsv"
    path.write_text(
        'nameID\tpopulation_freq\tpresence_by_short_json\n'
        'Gnosis-h1\t0.599\t{"Asfu_Af293": true, "Asfu_W72310": false}\n'
    )
    result = load_starships_by_short(path)
    assert result["Gnosis-h1"]["population_freq"] == 0.599
    assert result["Gnosis-h1"]["presence"] == {"Asfu_Af293": True, "Asfu_W72310": False}


def test_build_protein_to_family_reads_cluster_tsv(tmp_path):
    path = tmp_path / "cluster.tsv"
    path.write_text("repA\trepA\nrepA\tmemberA2\nrepB\trepB\n")
    result = build_protein_to_family(path)
    assert result == {"repA": "repA", "memberA2": "repA", "repB": "repB"}


def test_score_backend_end_to_end():
    # paper gene "geneA" crosswalks to study protein "Asfu_strainA|p1", which
    # is a member of family "famX" -- family famX is present in strainA and
    # strainB but absent in strainC.
    crosswalk = {"geneA": "Asfu_strainA|p1", "geneB": "Asfu_strainA|p2"}
    protein_to_family = {"Asfu_strainA|p1": "famX", "Asfu_strainA|p2": "famX"}
    ground_truth_cargo = {"Gnosis-h1": {"geneA", "geneB"}}
    ground_truth_starships = {
        "Gnosis-h1": {
            "population_freq": 0.5,
            "presence": {"Asfu_strainA": True, "Asfu_strainB": True, "Asfu_strainC": False},
        }
    }
    matrix = PresenceMatrix(families=["famX"], strains=["Asfu_strainA", "Asfu_strainB", "Asfu_strainC"])
    matrix.set_call("famX", "Asfu_strainA", "present")
    matrix.set_call("famX", "Asfu_strainB", "present")
    matrix.set_call("famX", "Asfu_strainC", "absent")

    rows = score_backend("mmseqs", crosswalk, protein_to_family, ground_truth_cargo, ground_truth_starships, matrix)
    cargo_row = next(r for r in rows if r["control_type"] == "positive_cargo_grouping")
    assert cargo_row["purity"] == 1.0
    assert cargo_row["completeness"] == 1.0
    presence_row = next(r for r in rows if r["control_type"] == "positive_presence_recovery")
    assert presence_row["diagnostic_family"] == "famX"
    assert presence_row["tp"] == 2 and presence_row["fn"] == 0 and presence_row["fp"] == 0


def test_score_negative_control_reports_frequency_and_bin(tmp_path):
    negative_path = tmp_path / "negative.tsv"
    negative_path.write_text(
        "gene\tlocus\taccession\tcluster\tpredictedFunction\n"
        "afumA\tAFUB_071550\tXP_1.1\tCluster_33\tSqualene hopane cyclase\n"
        "unmatchedGene\tAFUB_999999\tXP_9.9\tNaN\tsomething\n"
    )
    crosswalk = {"XP_1.1": "Asfu_strainA|p1"}
    protein_to_family = {"Asfu_strainA|p1": "famCore"}
    matrix = PresenceMatrix(families=["famCore"], strains=["Asfu_strainA", "Asfu_strainB"])
    matrix.set_call("famCore", "Asfu_strainA", "present")
    matrix.set_call("famCore", "Asfu_strainB", "present")

    rows = score_negative_control("mmseqs", crosswalk, protein_to_family, negative_path, matrix)
    matched = next(r for r in rows if r["entity_id"] == "XP_1.1")
    assert matched["frequency"] == 1.0
    assert matched["freq_bin"] == "core"
    unmatched = next(r for r in rows if r["entity_id"] == "XP_9.9")
    assert unmatched["notes"] == "not crosswalked"
