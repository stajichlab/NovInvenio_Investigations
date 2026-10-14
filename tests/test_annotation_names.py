"""Tests for lib/annotation_names.py (annotation file name per species + write guards)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from annotation_names import (  # noqa: E402
    annotation_names, annotation_path, ensure_within, is_safe_name, raw_annotation_name,
)


def row(short, acc, src="local_faa"):
    return {"Short": short, "Protein_Accession": acc, "Protein_Source": src}


@pytest.mark.parametrize("acc,expect", [
    ("UP000193719", "UP000193719"),                       # UniProt: unchanged
    ("GCA_002104975.1", "GCA_002104975.1"),               # NCBI: unchanged (the ".1" stays)
    ("/bigdata/x/GCA_032206225.1_FOAM13_A1/GCA_032206225.1_FOAM13_A1_protein.faa.gz",
     "GCA_032206225.1_FOAM13_A1_protein"),
    ("/bigdata/stajichlab/jpere468/k/GCF_050916275.1.faa", "GCF_050916275.1"),
    ("/x/UHM1088.41098__UHM1088.41098_R.bin.56.faa", "UHM1088.41098__UHM1088.41098_R.bin.56"),
    ("/x/Sp1.pep.fa", "Sp1"),
    ("/x/Sp1.fasta.zst", "Sp1"),
    ("/x/Sp1.fa", "Sp1"),
    ("rel/dir/Sp2.pep", "Sp2"),
    ("/x/.faa", ".faa"),                                   # nothing but a suffix: kept, then unsafe
    ("", ""),
])
def test_raw_annotation_name(acc, expect):
    assert raw_annotation_name(acc) == expect


def test_names_uniprot_ncbi_unchanged():
    rows = [row("A", "UP000193719", "uniprot"), row("B", "GCA_002104975.1", "ncbi")]
    assert annotation_names(rows) == {"A": "UP000193719", "B": "GCA_002104975.1"}


def test_names_paths_and_fallbacks():
    rows = [
        row("K1", "/a/GCF_1.faa"),
        row("K2", "/b/same.faa"),
        row("K3", "/c/same.faa.gz"),      # different path, same stem -> both fall back to Short
        row("K4", ""),                     # empty -> Short
        row("K5", "/d/.faa"),              # unsafe stem -> Short
        row("K6", "/a/GCF_1.faa"),         # same accession as K1: shares the file
    ]
    names = annotation_names(rows)
    assert names == {"K1": "GCF_1", "K2": "K2", "K3": "K3", "K4": "K4", "K5": "K5", "K6": "GCF_1"}


def test_uniprot_never_falls_back():
    rows = [row("U", "UP000000001", "uniprot"), row("L", "/x/UP000000001.faa")]
    names = annotation_names(rows)
    assert names["U"] == "UP000000001" and names["L"] == "L"


def test_fallback_collision_raises():
    rows = [row("GCF_1", "/x/a.faa"), row("Z", "/y/GCF_1.faa"), row("Q", "/z/a.faa")]
    # 'a' collides (two paths) -> GCF_1 falls back to its Short 'GCF_1', which Z already uses.
    with pytest.raises(ValueError, match="two different proteomes"):
        annotation_names(rows)


def test_is_safe_name():
    assert is_safe_name("GCA_1.1") and is_safe_name("UHM1088.41098__x_R.bin.56")
    for bad in ("", ".", "..", "a/b", "/abs", "a b", "-x", ".hidden"):
        assert not is_safe_name(bad)


def test_annotation_path_guard(tmp_path):
    ann = tmp_path / "annotations"
    ann.mkdir()
    assert annotation_path(ann, "GCF_1") == (ann / "GCF_1.tsv.gz").resolve()
    for bad in ("/etc/passwd", "../x", "a/b", ""):
        with pytest.raises(ValueError):
            annotation_path(ann, bad)


def test_ensure_within(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    assert ensure_within(a / "x" / "y.txt", [a, b]) == (a / "x" / "y.txt").resolve()
    assert ensure_within(b, [a, b]) == b.resolve()
    with pytest.raises(ValueError, match="outside"):
        ensure_within(tmp_path / "c" / "z", [a, b])
    with pytest.raises(ValueError, match="outside"):
        ensure_within(a / ".." / "c", [a])
    with pytest.raises(ValueError, match="outside"):
        ensure_within(Path("/bigdata/stajichlab/shared/x.tsv.gz"), [a])
