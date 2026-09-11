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
