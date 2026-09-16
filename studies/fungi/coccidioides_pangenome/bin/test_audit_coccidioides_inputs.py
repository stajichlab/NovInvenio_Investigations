import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from audit_coccidioides_inputs import list_freeze_strains


def test_list_freeze_strains_parses_species_and_strain(tmp_path):
    pep_dir = tmp_path / "pep"
    pep_dir.mkdir()
    (pep_dir / "Coccidioides_immitis_1M0.proteins.fa").touch()
    (pep_dir / "Coccidioides_posadasii_B3224.proteins.fa").touch()
    result = list_freeze_strains(pep_dir)
    assert result == {
        "1M0": "Coccidioides immitis",
        "B3224": "Coccidioides posadasii",
    }


def test_missing_annotation_strains_reported(tmp_path):
    from audit_coccidioides_inputs import missing_annotation_strains
    freeze_strains = {"1M0", "B3224"}
    samples_csv = tmp_path / "samples.csv"
    samples_csv.write_text(
        "RunAcc,Strain,BioSample,Center,Experiment,Project,Organism,FileBase,Notes,LocusTag\n"
        "r1,1M0,,,,,,,,\n"
        "r2,B3224,,,,,,,,\n"
        "r3,NOTANNOT1,,,,,,,,\n"
    )
    result = missing_annotation_strains(samples_csv, freeze_strains)
    assert result == {"NOTANNOT1"}


def test_busco_join_handles_aaftf_suffix(tmp_path):
    from audit_coccidioides_inputs import load_busco_complete, low_busco_strains
    asm_stats = tmp_path / "asm_stats.tsv"
    asm_stats.write_text(
        "SampleID\tBUSCO_Complete\n"
        "1M0.AAFTF\t97.2\n"
        "NM_9861.AAFTF\t0.2\n"
    )
    busco = load_busco_complete(asm_stats)
    assert busco == {"1M0": 97.2, "NM_9861": 0.2}
    assert low_busco_strains(busco, min_complete=90.0) == {"NM_9861"}


def test_busco_join_skips_na_and_blank_values(tmp_path):
    from audit_coccidioides_inputs import load_busco_complete
    asm_stats = tmp_path / "asm_stats.tsv"
    asm_stats.write_text(
        "SampleID\tBUSCO_Complete\n"
        "1M0.AAFTF\t97.2\n"
        "NOVAL.AAFTF\tNA\n"
        "BLANK.AAFTF\t\n"
    )
    busco = load_busco_complete(asm_stats)
    assert busco == {"1M0": 97.2}
