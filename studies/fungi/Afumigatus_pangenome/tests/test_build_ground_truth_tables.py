import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from build_ground_truth_tables import (
    normalize_id,
    normalize_geneid_for_crosswalk,
    build_short_lookup,
    match_isolates_to_short,
    cargo_gene_sets_by_nameid,
    conserved_negative_control_genes,
)


def test_normalize_geneid_for_crosswalk_converts_af293_xp_style():
    assert normalize_geneid_for_crosswalk("AF293_XP-001481524.1") == "XP_001481524.1"
    assert normalize_geneid_for_crosswalk("AF293_XP-753152.2") == "XP_753152.2"


def test_normalize_geneid_for_crosswalk_leaves_internal_locus_tags_untouched():
    # non-AF293 strains have no public accession -- nothing to normalize to
    assert normalize_geneid_for_crosswalk("47-10_000766") == "47-10_000766"


def test_normalize_id_strips_non_alnum_and_lowercases():
    assert normalize_id("W72310-lr") == "w72310lr"
    assert normalize_id("ATCC 42202") == "atcc42202"


def test_build_short_lookup_indexes_both_short_and_strain():
    config = pd.DataFrame({"Short": ["Asfu_W72310"], "Strain": ["W72310-lr"]})
    lookup = build_short_lookup(config)
    assert lookup["w72310lr"] == "Asfu_W72310"
    assert lookup["asfuw72310"] == "Asfu_W72310"


def test_match_isolates_to_short_falls_back_to_original_id():
    short_lookup = {"w72310lr": "Asfu_W72310", "af293": "Asfu_Af293"}
    isolate_ids = ["isoX", "isoY", "isoZ"]
    original_ids = ["W72310-lr", "AF293", "unmatched-strain"]
    result = match_isolates_to_short(isolate_ids, original_ids, short_lookup)
    assert result == {"isoX": "Asfu_W72310", "isoY": "Asfu_Af293"}
    assert "isoZ" not in result


def test_cargo_gene_sets_by_nameid_unions_across_starship_instances():
    table_s13 = pd.DataFrame({
        "nameID": ["Gnosis-h1", "Gnosis-h1", "Osiris-h4"],
        "starshipID": ["strainA_s1", "strainB_s1", "strainC_s2"],
        "geneID": ["geneA1", "geneB1", "geneC1"],
    })
    result = cargo_gene_sets_by_nameid(table_s13)
    assert result["Gnosis-h1"] == {"geneA1", "geneB1"}
    assert result["Osiris-h4"] == {"geneC1"}


def test_conserved_negative_control_genes_keeps_only_real_accessions():
    table_s19 = pd.DataFrame({
        "gene": ["afumA", "hacA", "noAccessionGene"],
        "locus": ["AFUB_071550", "Afu3g04070", "Afu9g99999"],
        "accession": ["XP_748727.1", "AspGD only", None],
        "cluster": ["Cluster_33", None, None],
        "predictedFunction": ["x", "y", "z"],
    })
    result = conserved_negative_control_genes(table_s19)
    assert result["gene"].tolist() == ["afumA"]
