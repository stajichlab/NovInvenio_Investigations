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
