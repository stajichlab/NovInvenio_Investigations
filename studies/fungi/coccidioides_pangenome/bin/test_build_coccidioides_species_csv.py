import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_coccidioides_species_csv import row_for_strain, TAXON_ID_BY_SPECIES


def test_row_for_strain_immitis():
    row = row_for_strain(
        strain="1M0",
        species="Coccidioides immitis",
        pep_path=Path("/x/pep/Coccidioides_immitis_1M0.proteins.fa"),
        dna_path=Path("/x/DNA/Coccidioides_immitis_1M0.scaffolds.fa"),
        gff_path=Path("/x/GFF/Coccidioides_immitis_1M0.gff3"),
    )
    assert row["Short"] == "1M0"
    assert row["Species"] == "Coccidioides immitis"
    assert row["Strain"] == "1M0"
    assert row["Group"] == "IN"
    assert row["TaxonGroup"] == "Coccidioides immitis"
    assert row["Taxon_ID"] == "5501"
    assert row["Protein_Source"] == "local_faa"
    assert row["Protein_Accession"] == "/x/pep/Coccidioides_immitis_1M0.proteins.fa"
    assert row["Genome_Source"] == "local_genome"
    assert row["GFF3_Source"] == "local_gff3"


def test_taxon_id_table_matches_verified_values():
    # 5501/199306 verified via `taxonkit lineage` against the installed NCBI
    # taxdump on 2026-09-15 -- see spec Open Question 8.
    assert TAXON_ID_BY_SPECIES == {
        "Coccidioides immitis": "5501",
        "Coccidioides posadasii": "199306",
    }
