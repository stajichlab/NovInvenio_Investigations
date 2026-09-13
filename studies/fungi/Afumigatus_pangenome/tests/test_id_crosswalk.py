import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from id_crosswalk import parse_diamond_blastp_besthits


def test_parse_diamond_blastp_besthits_keeps_highest_bitscore():
    # outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    lines = [
        "XP_748727.1\tstudy_protein_A\t99.0\t500\t0\t0\t1\t500\t1\t500\t0.0\t950",
        "XP_748727.1\tstudy_protein_B\t60.0\t300\t0\t0\t1\t300\t1\t300\t1e-20\t150",
        "AFUB_079030\tstudy_protein_C\t95.0\t400\t0\t0\t1\t400\t1\t400\t0.0\t800",
    ]
    result = parse_diamond_blastp_besthits(lines)
    assert result == {"XP_748727.1": "study_protein_A", "AFUB_079030": "study_protein_C"}
