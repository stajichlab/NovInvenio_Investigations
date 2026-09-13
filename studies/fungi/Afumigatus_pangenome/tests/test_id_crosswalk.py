import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from id_crosswalk import (
    parse_diamond_blastp_besthits,
    parse_diamond_blastp_hits_per_strain,
    parse_hmmsearch_tblout_hits,
)


def test_parse_diamond_blastp_besthits_keeps_highest_bitscore():
    # outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    lines = [
        "XP_748727.1\tstudy_protein_A\t99.0\t500\t0\t0\t1\t500\t1\t500\t0.0\t950",
        "XP_748727.1\tstudy_protein_B\t60.0\t300\t0\t0\t1\t300\t1\t300\t1e-20\t150",
        "AFUB_079030\tstudy_protein_C\t95.0\t400\t0\t0\t1\t400\t1\t400\t0.0\t800",
    ]
    result = parse_diamond_blastp_besthits(lines)
    assert result == {"XP_748727.1": "study_protein_A", "AFUB_079030": "study_protein_C"}


def test_parse_diamond_blastp_besthits_handles_non_monotonic_order():
    """Confirm highest bitscore wins regardless of file position (non-monotonic order)."""
    # outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    lines = [
        "XP_748727.1\tstudy_protein_A\t60.0\t300\t0\t0\t1\t300\t1\t300\t1e-20\t150",
        "XP_748727.1\tstudy_protein_B\t99.0\t500\t0\t0\t1\t500\t1\t500\t0.0\t950",
        "XP_748727.1\tstudy_protein_C\t90.0\t400\t0\t0\t1\t400\t1\t400\t1e-10\t300",
    ]
    result = parse_diamond_blastp_besthits(lines)
    # Bitscores in order: 150, 950, 300 — highest (950) should win despite not being first
    assert result == {"XP_748727.1": "study_protein_B"}


def test_parse_diamond_blastp_hits_per_strain_keeps_one_hit_per_query_strain():
    lines = [
        "hacA\tAsfu_strainA|KAH1.1\t100.0\t433\t0\t0\t1\t433\t1\t433\t1e-250\t731",
        # a second, worse hit for the SAME (query, strain) -- should be dropped
        "hacA\tAsfu_strainA|KAH2.1\t60.0\t200\t0\t0\t1\t200\t1\t200\t1e-20\t150",
        "hacA\tAsfu_strainB|KAH3.1\t98.0\t430\t0\t0\t1\t430\t1\t430\t1e-240\t700",
    ]
    result = parse_diamond_blastp_hits_per_strain(lines)
    assert set(result.keys()) == {("hacA", "Asfu_strainA"), ("hacA", "Asfu_strainB")}
    assert result[("hacA", "Asfu_strainA")]["subject"] == "Asfu_strainA|KAH1.1"
    assert result[("hacA", "Asfu_strainA")]["bitscore"] == 731.0
    assert result[("hacA", "Asfu_strainB")]["pident"] == 98.0


def test_parse_diamond_blastp_hits_per_strain_raises_when_all_ids_unprefixed():
    lines = ["hacA\tbare_protein_id\t100.0\t433\t0\t0\t1\t433\t1\t433\t1e-250\t731"]
    try:
        parse_diamond_blastp_hits_per_strain(lines)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_parse_hmmsearch_tblout_hits_groups_by_strain():
    # hmmsearch --tblout: target name, accession, query name, accession,
    # E-value(full seq), score(full seq), bias, ... (whitespace-separated)
    lines = [
        "# comment line, ignored",
        "Asfu_strainA|KAH1.1 -          PF11001  PF11001.15   1.2e-40  135.0   0.3   1.3e-40  134.6   0.2   1   1   0   0   1   1   0   0 -",
        "Asfu_strainA|KAH2.1 -          PF11001  PF11001.15   3.0e-10   30.0   0.1   4.0e-10   29.6   0.1   1   1   0   0   1   1   0   0 -",
        "Asfu_strainB|KAH3.1 -          PF11001  PF11001.15   5.5e-35  120.0   0.2   6.0e-35  119.8   0.2   1   1   0   0   1   1   0   0 -",
    ]
    result = parse_hmmsearch_tblout_hits(lines)
    assert set(result.keys()) == {"Asfu_strainA", "Asfu_strainB"}
    assert len(result["Asfu_strainA"]) == 2
    assert len(result["Asfu_strainB"]) == 1
    assert result["Asfu_strainB"][0]["evalue"] == 5.5e-35
