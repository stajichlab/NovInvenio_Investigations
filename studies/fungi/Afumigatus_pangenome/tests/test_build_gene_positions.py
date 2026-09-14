import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from build_gene_positions import parse_gff3_protein_positions


def test_parse_gff3_protein_positions_single_exon(tmp_path):
    gff3 = tmp_path / "strain.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "contig1\tGenbank\tgene\t100\t500\t.\t+\t.\tID=gene-A\n"
        "contig1\tGenbank\tCDS\t100\t500\t.\t+\t0\tID=cds-P1;protein_id=P1.1\n"
    )
    positions = parse_gff3_protein_positions(gff3)
    assert positions == {"P1.1": ("contig1", 100, 500)}


def test_parse_gff3_protein_positions_collapses_multi_exon_cds():
    # Same protein_id across two CDS rows (an intron between them) --
    # must collapse to (min start, max end), not two separate entries.
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
        f.write(
            "contig1\tGenbank\tCDS\t3203\t3229\t.\t+\t0\tprotein_id=P2.1\n"
            "contig1\tGenbank\tCDS\t3290\t3375\t.\t+\t0\tprotein_id=P2.1\n"
        )
        path = f.name
    positions = parse_gff3_protein_positions(path)
    assert positions == {"P2.1": ("contig1", 3203, 3375)}


def test_parse_gff3_protein_positions_ignores_non_cds_features():
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
        f.write(
            "contig1\tGenbank\tgene\t1\t1000\t.\t+\t.\tID=gene-A;protein_id=SHOULD_NOT_APPEAR\n"
            "contig1\tGenbank\tmRNA\t1\t1000\t.\t+\t.\tID=rna-A;protein_id=SHOULD_NOT_APPEAR\n"
            "contig1\tGenbank\tCDS\t10\t50\t.\t+\t0\tprotein_id=REAL.1\n"
        )
        path = f.name
    positions = parse_gff3_protein_positions(path)
    assert positions == {"REAL.1": ("contig1", 10, 50)}


def test_parse_gff3_protein_positions_skips_comments_and_blank_lines():
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
        f.write(
            "##gff-version 3\n"
            "\n"
            "# a comment\n"
            "contig1\tGenbank\tCDS\t10\t50\t.\t+\t0\tprotein_id=REAL.1\n"
        )
        path = f.name
    positions = parse_gff3_protein_positions(path)
    assert positions == {"REAL.1": ("contig1", 10, 50)}
