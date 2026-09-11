import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
import build_study_config as bsc  # noqa: E402


def _write_species_csv(path, rows, header=None):
    header = header or [
        "Short", "Species", "Strain", "Group", "TaxonGroup",
        "Protein_Source", "Protein_Accession", "Taxon_ID",
        "Genome_Source", "Genome_Accession",
        "GFF3_Source", "GFF3_Accession",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


def test_local_faa_and_local_genome_row(tmp_path, monkeypatch):
    # Arrange: a study dir whose only row sources everything from local files.
    study_dir = tmp_path / "studies" / "bacteria" / "toy_study"
    study_dir.mkdir(parents=True)
    local_faa = tmp_path / "src" / "Sp1.faa"
    local_faa.parent.mkdir(parents=True)
    local_faa.write_text(">seq1\nMAAA\n")
    local_genome = tmp_path / "src" / "Sp1.fna"
    local_genome.write_text(">contig1\nACGT\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp1", "Species": "Test species", "Strain": "T1",
        "Group": "IN", "TaxonGroup": "TestGroup",
        "Protein_Source": "local_faa", "Protein_Accession": str(local_faa),
        "Taxon_ID": "",
        "Genome_Source": "local_genome", "Genome_Accession": str(local_genome),
        "GFF3_Source": "", "GFF3_Accession": "",
    }])

    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(tmp_path / "data/ncbi"),
    ])

    # Act
    rc = bsc.main()

    # Assert
    assert rc == 0
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    row = rows[0]
    assert row["Short"] == "Sp1"
    assert row["Protein"] == "Sp1.pep.fa"
    assert row["DNA"] == "Sp1.dna.fa"
    assert row["GFF3"] == ""
    assert (study_dir / "data_dir" / "pep" / "Sp1.pep.fa").read_text() == ">seq1\nMAAA\n"
    assert (study_dir / "data_dir" / "dna" / "Sp1.dna.fa").read_text() == ">contig1\nACGT\n"
    manifest = (study_dir / "DATA_MANIFEST.yaml").read_text()
    assert "local_path:" in manifest
    assert "source_url:" in manifest
    assert str(local_faa) in manifest  # source path recorded in source_url field


def test_local_faa_with_ncbi_genome_row(tmp_path, monkeypatch):
    # Arrange: protein already local, genome must be "fetched" -- mock the
    # subprocess call and pre-seed the ncbi cache the way fetch_genome_assembly.py
    # would, so we test build_study_config.py's own cache-reading logic, not the
    # network fetch itself.
    study_dir = tmp_path / "studies" / "bacteria" / "toy_study2"
    study_dir.mkdir(parents=True)
    local_faa = tmp_path / "src" / "Sp2.faa"
    local_faa.parent.mkdir(parents=True)
    local_faa.write_text(">seq2\nMBBB\n")

    ncbi_cache = tmp_path / "data/ncbi"
    genome_dir = ncbi_cache / "GCF_000000001.1" / "extracted" / "ncbi_dataset" / "data" / "GCF_000000001.1"
    genome_dir.mkdir(parents=True)
    (genome_dir / "genomic.fna").write_text(">contig2\nTTTT\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp2", "Species": "Test species 2", "Strain": "T2",
        "Group": "OUT", "TaxonGroup": "TestGroup",
        "Protein_Source": "local_faa", "Protein_Accession": str(local_faa),
        "Taxon_ID": "",
        "Genome_Source": "ncbi", "Genome_Accession": "GCF_000000001.1",
        "GFF3_Source": "", "GFF3_Accession": "",
    }])

    calls = []

    def fake_run(cmd):
        calls.append(cmd)

    monkeypatch.setattr(bsc, "run", fake_run)
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(ncbi_cache),
    ])

    rc = bsc.main()

    assert rc == 0
    assert any("fetch_genome_assembly.py" in str(c) for call in calls for c in call)
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        row = next(csv.DictReader(fh))
    assert row["DNA"] == "Sp2.dna.fa"
    assert (study_dir / "data_dir" / "dna" / "Sp2.dna.fa").read_text() == ">contig2\nTTTT\n"


def test_unknown_source_errors(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "bacteria" / "bad_study"
    study_dir.mkdir(parents=True)
    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp3", "Species": "x", "Strain": "", "Group": "IN",
        "TaxonGroup": "x",
        "Protein_Source": "bogus", "Protein_Accession": "whatever",
        "Taxon_ID": "", "Genome_Source": "", "Genome_Accession": "",
        "GFF3_Source": "", "GFF3_Accession": "",
    }])
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py", "--study-dir", str(study_dir),
    ])
    try:
        bsc.main()
        assert False, "expected SystemExit"
    except SystemExit as e:
        assert "bogus" in str(e)


def test_ncbi_genome_with_local_gff3(tmp_path, monkeypatch):
    # Arrange: NCBI genome but local GFF3 (local overrides NCBI's GFF3)
    study_dir = tmp_path / "studies" / "bacteria" / "toy_study3"
    study_dir.mkdir(parents=True)
    local_faa = tmp_path / "src" / "Sp3.faa"
    local_faa.parent.mkdir(parents=True)
    local_faa.write_text(">seq3\nMCCC\n")
    local_gff3 = tmp_path / "src" / "Sp3.gff3"
    local_gff3.write_text("##gff-version 3\ncontig1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

    ncbi_cache = tmp_path / "data/ncbi"
    genome_dir = ncbi_cache / "GCF_000000002.1" / "extracted" / "ncbi_dataset" / "data" / "GCF_000000002.1"
    genome_dir.mkdir(parents=True)
    (genome_dir / "genomic.fna").write_text(">contig3\nGGGG\n")
    (genome_dir / "genomic.gff").write_text("##gff-version 3\ncontig3\t.\tgene\t1\t50\t.\t+\t.\tID=ncbi_gene1\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp3", "Species": "Test species 3", "Strain": "T3",
        "Group": "IN", "TaxonGroup": "TestGroup",
        "Protein_Source": "local_faa", "Protein_Accession": str(local_faa),
        "Taxon_ID": "",
        "Genome_Source": "ncbi", "Genome_Accession": "GCF_000000002.1",
        "GFF3_Source": "local_gff3", "GFF3_Accession": str(local_gff3),
    }])

    calls = []

    def fake_run(cmd):
        calls.append(cmd)

    monkeypatch.setattr(bsc, "run", fake_run)
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(ncbi_cache),
    ])

    rc = bsc.main()

    assert rc == 0
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        row = next(csv.DictReader(fh))
    # Verify genome and GFF3 were used
    assert row["DNA"] == "Sp3.dna.fa"
    assert row["GFF3"] == "Sp3.gff3"
    assert (study_dir / "data_dir" / "dna" / "Sp3.dna.fa").read_text() == ">contig3\nGGGG\n"
    # Verify local GFF3 was copied (not NCBI's)
    assert (study_dir / "data_dir" / "gff3" / "Sp3.gff3").read_text() == "##gff-version 3\ncontig1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n"
    # Verify both local and NCBI provenance are in manifest
    manifest = (study_dir / "DATA_MANIFEST.yaml").read_text()
    assert str(local_gff3) in manifest  # local GFF3 source recorded
    assert "source_url:" in manifest


def test_ncbi_genome_with_gff3_none(tmp_path, monkeypatch):
    # Arrange: NCBI genome available but explicitly declined (GFF3_Source=none).
    # This tests the bug fix: provenance should NOT include the declined GFF3.
    study_dir = tmp_path / "studies" / "bacteria" / "toy_study4"
    study_dir.mkdir(parents=True)
    local_faa = tmp_path / "src" / "Sp4.faa"
    local_faa.parent.mkdir(parents=True)
    local_faa.write_text(">seq4\nMDDD\n")

    ncbi_cache = tmp_path / "data/ncbi"
    genome_dir = ncbi_cache / "GCF_000000003.1" / "extracted" / "ncbi_dataset" / "data" / "GCF_000000003.1"
    genome_dir.mkdir(parents=True)
    (genome_dir / "genomic.fna").write_text(">contig4\nAAAA\n")
    (genome_dir / "genomic.gff").write_text("##gff-version 3\ncontig4\t.\tgene\t1\t30\t.\t+\t.\tID=ncbi_gene\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp4", "Species": "Test species 4", "Strain": "T4",
        "Group": "OUT", "TaxonGroup": "TestGroup",
        "Protein_Source": "local_faa", "Protein_Accession": str(local_faa),
        "Taxon_ID": "",
        "Genome_Source": "ncbi", "Genome_Accession": "GCF_000000003.1",
        "GFF3_Source": "none", "GFF3_Accession": "",
    }])

    calls = []

    def fake_run(cmd):
        calls.append(cmd)

    monkeypatch.setattr(bsc, "run", fake_run)
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(ncbi_cache),
    ])

    rc = bsc.main()

    assert rc == 0
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        row = next(csv.DictReader(fh))
    # Verify genome was used but GFF3 was declined
    assert row["DNA"] == "Sp4.dna.fa"
    assert row["GFF3"] == ""
    assert (study_dir / "data_dir" / "dna" / "Sp4.dna.fa").read_text() == ">contig4\nAAAA\n"
    # Verify GFF3 directory was NOT created (file was not copied)
    gff3_file = study_dir / "data_dir" / "gff3" / "Sp4.gff3"
    assert not gff3_file.exists(), "GFF3 file should not exist when GFF3_Source=none"
    # BUG FIX REGRESSION TEST: Verify manifest does NOT contain provenance for the declined NCBI GFF3
    manifest = (study_dir / "DATA_MANIFEST.yaml").read_text()
    # Count lines with "genomic.gff" -- should be 0 (not used, so no provenance)
    gff_lines = [line for line in manifest.split("\n") if "genomic.gff" in line]
    assert len(gff_lines) == 0, f"Manifest should not contain provenance for declined GFF3, but found: {gff_lines}"


def test_ncbi_sourced_protein(tmp_path, monkeypatch):
    # Arrange: Protein_Source=ncbi (protein from same NCBI Datasets genome package)
    study_dir = tmp_path / "studies" / "fungi" / "toy_study5"
    study_dir.mkdir(parents=True)

    ncbi_cache = tmp_path / "data/ncbi"
    genome_dir = ncbi_cache / "GCF_000005845.2" / "extracted" / "ncbi_dataset" / "data" / "GCF_000005845.2"
    genome_dir.mkdir(parents=True)
    (genome_dir / "genomic.fna").write_text(">chr1\nCCCC\n")
    (genome_dir / "genomic.gff").write_text("##gff-version 3\nchr1\t.\tgene\t1\t20\t.\t+\t.\tID=gene1\n")
    (genome_dir / "protein.faa").write_text(">prot1\nMEEE\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp5", "Species": "Test species 5", "Strain": "T5",
        "Group": "IN", "TaxonGroup": "TestGroup",
        "Protein_Source": "ncbi", "Protein_Accession": "GCF_000005845.2",
        "Taxon_ID": "",
        "Genome_Source": "ncbi", "Genome_Accession": "GCF_000005845.2",
        "GFF3_Source": "", "GFF3_Accession": "",
    }])

    calls = []

    def fake_run(cmd):
        calls.append(cmd)

    monkeypatch.setattr(bsc, "run", fake_run)
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(ncbi_cache),
    ])

    rc = bsc.main()

    assert rc == 0
    # Verify fetch_genome_assembly.py was called with --include-protein (resolve_protein calls it for NCBI protein fetch)
    assert any("--include-protein" in call for call in calls), \
        f"Expected fetch_genome_assembly.py --include-protein in calls, got: {calls}"
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        row = next(csv.DictReader(fh))
    # Verify protein, genome, and GFF3 are all present from NCBI
    assert row["Short"] == "Sp5"
    assert row["Protein"] == "Sp5.pep.fa"
    assert row["DNA"] == "Sp5.dna.fa"
    assert row["GFF3"] == "Sp5.gff3"
    # Verify files were copied correctly
    assert (study_dir / "data_dir" / "pep" / "Sp5.pep.fa").read_text() == ">prot1\nMEEE\n"
    assert (study_dir / "data_dir" / "dna" / "Sp5.dna.fa").read_text() == ">chr1\nCCCC\n"
    assert (study_dir / "data_dir" / "gff3" / "Sp5.gff3").read_text() == "##gff-version 3\nchr1\t.\tgene\t1\t20\t.\t+\t.\tID=gene1\n"
