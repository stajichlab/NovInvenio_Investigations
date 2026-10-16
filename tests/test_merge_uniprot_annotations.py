import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin" / "merge_uniprot_annotations.py"


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def test_merge_passes_through_xrefs_column(tmp_path):
    annot = tmp_path / "annot.tsv"
    _write_tsv(
        annot,
        [{"accession": "A7UWL5", "gene_name": "", "description": "", "pfam_ids": "",
          "interpro_ids": "", "go_ids": "", "xrefs": "GeneID:5847462|KEGG:ncr:NCU10683"}],
        ["accession", "gene_name", "description", "pfam_ids", "interpro_ids", "go_ids", "xrefs"],
    )
    matrix = tmp_path / "matrix.tsv"
    _write_tsv(matrix, [{"protein_id": "A7UWL5"}], ["protein_id"])
    out = tmp_path / "out.tsv"

    subprocess.run(
        [sys.executable, str(BIN), "--matrix", str(matrix), "--annotations", str(annot), "--output", str(out)],
        check=True,
    )

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["uniprot_xrefs"] == "GeneID:5847462|KEGG:ncr:NCU10683"


def test_merge_lenient_missing_xrefs_column_is_empty_string(tmp_path):
    # An older annotations TSV (predates Task 1) has no xrefs column at all --
    # merge_uniprot_annotations.py must not KeyError on it.
    annot = tmp_path / "annot.tsv"
    _write_tsv(
        annot,
        [{"accession": "A7UWL5", "gene_name": "", "description": "", "pfam_ids": "",
          "interpro_ids": "", "go_ids": ""}],
        ["accession", "gene_name", "description", "pfam_ids", "interpro_ids", "go_ids"],
    )
    matrix = tmp_path / "matrix.tsv"
    _write_tsv(matrix, [{"protein_id": "A7UWL5"}], ["protein_id"])
    out = tmp_path / "out.tsv"

    subprocess.run(
        [sys.executable, str(BIN), "--matrix", str(matrix), "--annotations", str(annot), "--output", str(out)],
        check=True,
    )

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["uniprot_xrefs"] == ""


# --- annotation_source column (2026-10-14) -------------------------------------
import gzip  # noqa: E402

UNIPROT_COLS = ["accession", "taxon_id", "gene_name", "description", "go_ids", "pfam_ids",
                "pfam_names", "interpro_ids", "ec_numbers", "alphafold_id", "xrefs"]
HMMER_COLS = UNIPROT_COLS + ["pfam_best_ievalue", "pfam_best_domain_score"]
NOTE_BIN = REPO / "bin" / "annotation_source_note.py"


def _uniprot_file(path: Path) -> Path:
    _write_tsv(path, [
        {"accession": "P00498", "taxon_id": "559292", "gene_name": "HIS1", "description": "ATP PRT",
         "go_ids": "GO:0005737:IEA", "pfam_ids": "PF01634", "pfam_names": "HisG",
         "interpro_ids": "IPR001348", "ec_numbers": "2.4.2.17", "alphafold_id": "P00498", "xrefs": ""},
    ], UNIPROT_COLS)
    return path


def _hmmer_file(path: Path, marker_column: bool = True) -> Path:
    cols = HMMER_COLS if marker_column else UNIPROT_COLS
    row = {c: "" for c in cols}
    row.update({"accession": "OUM56101.1", "taxon_id": "73868", "description": "hypothetical protein",
                "pfam_ids": "PF16116", "pfam_names": "DUF4832"})
    if marker_column:
        row.update({"pfam_best_ievalue": "7.9e-24", "pfam_best_domain_score": "85.3"})
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerow(row)
    return path


def _merge(tmp_path, annots, matrix_ids, extra_args=()):
    matrix = tmp_path / "matrix.tsv"
    _write_tsv(matrix, [{"protein_id": p, "S1": "1"} for p in matrix_ids], ["protein_id", "S1"])
    out = tmp_path / "out.tsv"
    subprocess.run([sys.executable, str(BIN), "--matrix", str(matrix), "--annotations",
                    *[str(a) for a in annots], "--output", str(out), *extra_args], check=True)
    with open(out, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames), list(reader), out


def test_uniprot_only_adds_source_as_last_column(tmp_path):
    header, rows, out = _merge(tmp_path, [_uniprot_file(tmp_path / "UP1.tsv")], ["sp|P00498|HIS1_YEAST"])
    assert header == ["protein_id", "S1", "uniprot_gene_name", "uniprot_description", "uniprot_pfam_ids",
                      "uniprot_pfam_names", "uniprot_interpro_ids", "uniprot_go_ids", "uniprot_ec_numbers",
                      "uniprot_alphafold_id", "uniprot_xrefs", "annotation_source"]
    assert rows[0]["annotation_source"] == "uniprot"
    assert rows[0]["uniprot_go_ids"] == "GO:0005737:IEA"
    # Every other byte is as before: drop the last field of each line and compare
    # with the pre-2026-10-14 layout (DictWriter's \r\n line ends).
    lines = out.read_bytes().split(b"\r\n")
    assert lines[1].rsplit(b"\t", 1) == [
        b"sp|P00498|HIS1_YEAST\t1\tHIS1\tATP PRT\tPF01634\tHisG\tIPR001348\tGO:0005737:IEA\t2.4.2.17\tP00498\t",
        b"uniprot"]


def test_mixed_merge_labels_each_row_and_keeps_missing_empty(tmp_path):
    annots = [_uniprot_file(tmp_path / "UP1.tsv"), _hmmer_file(tmp_path / "GCA_1.tsv.gz")]
    header, rows, _ = _merge(tmp_path, annots, ["sp|P00498|HIS1_YEAST", "OUM56101.1", "XYZ.1"],
                             extra_args=["--lenient"])
    assert header[-1] == "annotation_source"
    by_id = {r["protein_id"]: r for r in rows}
    assert by_id["sp|P00498|HIS1_YEAST"]["annotation_source"] == "uniprot"
    assert by_id["OUM56101.1"]["annotation_source"] == "hmmer-pfam"
    assert by_id["OUM56101.1"]["uniprot_pfam_ids"] == "PF16116"
    assert by_id["OUM56101.1"]["uniprot_go_ids"] == ""
    assert by_id["XYZ.1"]["annotation_source"] == ""


def test_pfam_only_species_detected_from_sidecar_without_marker_column(tmp_path):
    f = _hmmer_file(tmp_path / "GCA_2.tsv.gz", marker_column=False)
    (tmp_path / "GCA_2.tsv.gz.provenance.yaml").write_text(
        "derived_by: 'bin/hmm_domtbl_to_annotations.py: python3 hmm_domtbl_to_annotations.py --domtbl x'\n")
    _, rows, _ = _merge(tmp_path, [f], ["OUM56101.1"])
    assert rows[0]["annotation_source"] == "hmmer-pfam"


def test_pfam_only_species_from_marker_column(tmp_path):
    _, rows, _ = _merge(tmp_path, [_hmmer_file(tmp_path / "GCA_3.tsv.gz")], ["OUM56101.1"])
    assert rows[0]["annotation_source"] == "hmmer-pfam"


PAGE = "<html><body>\n<div class=\"wrap\">\n  <header class=\"top\">\n    <h1 id=\"title\"></h1>\n  </header>\n  <main></main>\n</body></html>\n"


def _note(tmp_path, annots, pages):
    species = tmp_path / "species.csv"
    species.write_text("Short,Species,Protein_Source,Protein_Accession\n"
                       "Sacccere,Saccharomyces cerevisiae,uniprot,UP1\n"
                       "PiroE2sp,Piromyces sp. E2,ncbi,GCA_1\n")
    subprocess.run([sys.executable, str(NOTE_BIN), "--annotations", *[str(a) for a in annots],
                    "--species", str(species), "--html", *[str(p) for p in pages]], check=True)


def test_note_names_hmmer_species_and_is_idempotent(tmp_path):
    page = tmp_path / "novelties.html"
    page.write_text(PAGE)
    annots = [_uniprot_file(tmp_path / "UP1.tsv"), _hmmer_file(tmp_path / "GCA_1.tsv.gz")]
    _note(tmp_path, annots, [page, tmp_path / "losses.html"])  # missing page is skipped
    once = page.read_text()
    assert "Pfam from a HMMER search (no GO/InterPro) for: Piromyces sp. E2 (PiroE2sp)." in once
    assert "Saccharomyces" not in once
    assert once.index("</header>") < once.index("nii-annotation-source:start") < once.index("<main>")
    _note(tmp_path, annots, [page])
    assert page.read_text() == once


def test_note_absent_and_removed_for_uniprot_only(tmp_path):
    page = tmp_path / "core.html"
    page.write_text(PAGE)
    _note(tmp_path, [_uniprot_file(tmp_path / "UP1.tsv")], [page])
    assert page.read_text() == PAGE
    # A note left by an earlier mixed sync is removed.
    _note(tmp_path, [_uniprot_file(tmp_path / "UP1.tsv"), _hmmer_file(tmp_path / "GCA_1.tsv.gz")], [page])
    assert "nii-annotation-source" in page.read_text()
    _note(tmp_path, [_uniprot_file(tmp_path / "UP1.tsv")], [page])
    assert page.read_text() == PAGE


def test_note_from_matrix_names_only_species_with_hmmer_rows(tmp_path):
    # Merged matrix: one UniProt row (Sacccere) and one HMMER row (PiroE2sp).
    matrix = tmp_path / "presence_matrix.uniprot.tsv"
    _write_tsv(matrix, [
        {"protein_id": "sp|P00498|HIS1_YEAST", "source_proteome": "Sacccere", "annotation_source": "uniprot"},
        {"protein_id": "OUM56101.1", "source_proteome": "PiroE2sp", "annotation_source": "hmmer-pfam"},
    ], ["protein_id", "source_proteome", "annotation_source"])
    loss = tmp_path / "loss_presence_matrix.uniprot.tsv"
    _write_tsv(loss, [{"protein_id": "sp|P00498|HIS1_YEAST", "source_proteome": "Sacccere",
                       "annotation_source": "uniprot"}], ["protein_id", "source_proteome", "annotation_source"])
    species = tmp_path / "species.csv"
    species.write_text("Short,Species,Protein_Source,Protein_Accession\n"
                       "Sacccere,Saccharomyces cerevisiae,uniprot,UP1\n"
                       "PiroE2sp,Piromyces sp. E2,ncbi,GCA_1\n")
    nov, losses = tmp_path / "novelties.html", tmp_path / "losses.html"
    nov.write_text(PAGE)
    losses.write_text(PAGE)
    for m, page in ((matrix, nov), (loss, losses)):
        subprocess.run([sys.executable, str(NOTE_BIN), "--matrix", str(m), "--species", str(species),
                        "--html", str(page)], check=True)
    assert "for: Piromyces sp. E2 (PiroE2sp)." in nov.read_text()
    assert losses.read_text() == PAGE


# --- --fill-empty-uniprot files and local-FASTA names (2026-10-16) --------------
import hashlib  # noqa: E402

sys.path.insert(0, str(REPO / "lib"))
from annotation_sources import file_accession, file_source, species_labels  # noqa: E402

FILL_SIDECAR = """derived_by: 'bin/annotate_proteomes_pfam.py --fill-empty-uniprot: python3 bin/annotate_proteomes_pfam.py
  --study fungi/pezizo_set1 --species Amega --finalize --fill-empty-uniprot'
annotation_source: UniProt reference proteome (no GO/Pfam/InterPro in the UniProt
  entries) + Pfam-A HMM search (hmmsearch) for pfam_ids/pfam_names
checksum_sha256: {sha}
vocabularies_available:
- pfam
- pfam_clan (derived from Pfam-A.hmm.dat by the reader)
vocabularies_not_available:
- go
- interpro
"""


def _filled_uniprot_file(path: Path) -> Path:
    """UniProt columns only (no pfam_best_ievalue); pfam filled, GO/InterPro empty."""
    _hmmer_file(path, marker_column=False)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name(path.name + ".provenance.yaml").write_text(FILL_SIDECAR.format(sha=sha))
    return path


def test_fill_empty_uniprot_file_is_hmmer_pfam(tmp_path):
    f = _filled_uniprot_file(tmp_path / "UP001658139.tsv.gz")
    assert file_source(f) == "hmmer-pfam"
    _, rows, _ = _merge(tmp_path, [_uniprot_file(tmp_path / "UP1.tsv"), f],
                        ["sp|P00498|HIS1_YEAST", "OUM56101.1"])
    by_id = {r["protein_id"]: r["annotation_source"] for r in rows}
    assert by_id == {"sp|P00498|HIS1_YEAST": "uniprot", "OUM56101.1": "hmmer-pfam"}


def test_fill_empty_sidecar_is_ignored_after_file_is_rewritten(tmp_path):
    # build_study_config.py rewrites the file from the .dat.gz; the old sidecar
    # stays. Checksum mismatch -> the file is UniProt again.
    f = _filled_uniprot_file(tmp_path / "UP001658139.tsv.gz")
    plain = _uniprot_file(tmp_path / "rewritten.tsv")
    f.write_bytes(gzip.compress(plain.read_bytes()))
    assert file_source(f) == "uniprot"


def test_fill_empty_sidecar_needs_go_not_available(tmp_path):
    f = _filled_uniprot_file(tmp_path / "UP9.tsv.gz")
    side = tmp_path / "UP9.tsv.gz.provenance.yaml"
    side.write_text(side.read_text().replace("- go\n", ""))
    assert file_source(f) == "uniprot"


LOCAL_SPECIES = ("Short,Species,Protein_Source,Protein_Accession\n"
                 "KoxO_050916275,Klebsiella oxytoca,local_faa,/x/outgroup_faa/GCF_050916275.1.faa\n"
                 "UHM22bin8,Klebsiella oxytoca,local_faa,/x/bins/UHM22.bin.8.faa.gz\n"
                 "Sacccere,Saccharomyces cerevisiae,uniprot,UP1\n")


def test_species_labels_map_local_fasta_annotation_names(tmp_path):
    species = tmp_path / "species.csv"
    species.write_text(LOCAL_SPECIES)
    labels = species_labels(species)
    assert file_accession(tmp_path / "GCF_050916275.1.tsv.gz") == "GCF_050916275.1"
    assert labels["GCF_050916275.1"] == "Klebsiella oxytoca (KoxO_050916275)"
    assert labels["UHM22.bin.8"] == "Klebsiella oxytoca (UHM22bin8)"
    assert labels["UP1"] == "Saccharomyces cerevisiae (Sacccere)"
    assert labels["/x/outgroup_faa/GCF_050916275.1.faa"] == "Klebsiella oxytoca (KoxO_050916275)"


def test_note_names_local_fasta_species_from_annotation_files(tmp_path):
    species = tmp_path / "species.csv"
    species.write_text(LOCAL_SPECIES)
    page = tmp_path / "novelties.html"
    page.write_text(PAGE)
    annots = [_uniprot_file(tmp_path / "UP1.tsv"), _hmmer_file(tmp_path / "UHM22.bin.8.tsv.gz")]
    subprocess.run([sys.executable, str(NOTE_BIN), "--annotations", *[str(a) for a in annots],
                    "--species", str(species), "--html", str(page)], check=True)
    assert "for: Klebsiella oxytoca (UHM22bin8)." in page.read_text()
