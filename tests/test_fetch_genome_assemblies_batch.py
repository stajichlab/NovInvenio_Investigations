import subprocess
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
import fetch_genome_assemblies_batch as fgab  # noqa: E402


def test_chunked_splits_by_size():
    assert fgab.chunked(["a", "b", "c", "d", "e"], 2) == [["a", "b"], ["c", "d"], ["e"]]
    assert fgab.chunked(["a", "b"], 75) == [["a", "b"]]
    assert fgab.chunked([], 75) == []


def _write_fake_batch_zip(zip_path: Path, accessions: list[str], *, include_protein: bool) -> None:
    """Mimic a real `datasets download genome accession <acc...>` zip: one
    ncbi_dataset/data/<accession>/ subfolder per accession."""
    with zipfile.ZipFile(zip_path, "w") as zf:
        for acc in accessions:
            zf.writestr(f"ncbi_dataset/data/{acc}/{acc}_genomic.fna", f">contig\nACGT{acc}\n")
            zf.writestr(f"ncbi_dataset/data/{acc}/genomic.gff", "##gff-version 3\n")
            if include_protein:
                zf.writestr(f"ncbi_dataset/data/{acc}/protein.faa", f">p1\nMMMM{acc}\n")


def test_fetch_batch_chunks_and_reshapes_per_accession(tmp_path, monkeypatch):
    outdir = tmp_path / "ncbi_cache"
    accessions = ["GCA_1", "GCA_2", "GCA_3"]

    monkeypatch.setattr(fgab.shutil, "which", lambda name: "/usr/bin/fake_datasets")
    monkeypatch.setattr(fgab, "datasets_version", lambda exe: "fake-version 1.0")

    calls = []

    def fake_run(cmd, check=True, **kwargs):
        calls.append(cmd)
        if cmd[1:3] == ["download", "genome"]:
            # cmd shape: [exe, "download", "genome", "accession", *batch, "--include", ..., "--filename", path]
            batch = cmd[4 : cmd.index("--include")]
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            _write_fake_batch_zip(zip_path, batch, include_protein=True)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(fgab.subprocess, "run", fake_run)

    # batch_size=2 forces 2 chunks for 3 accessions: [GCA_1, GCA_2], [GCA_3]
    fgab.fetch_batch(accessions, outdir=outdir, include_protein=True, batch_size=2)

    download_calls = [c for c in calls if c[1:3] == ["download", "genome"]]
    assert len(download_calls) == 2, "expected 2 chunked datasets calls for batch_size=2 over 3 accessions"

    for acc in accessions:
        data_dir = outdir / acc / "extracted" / "ncbi_dataset" / "data" / acc
        fna = list(data_dir.glob("*.fna"))
        gff = list(data_dir.glob("*.gff"))
        faa = list(data_dir.glob("*.faa"))
        assert len(fna) == 1, f"{acc}: expected exactly one .fna"
        assert len(gff) == 1, f"{acc}: expected exactly one .gff"
        assert len(faa) == 1, f"{acc}: expected exactly one .faa"
        assert fna[0].with_suffix(fna[0].suffix + ".provenance.yaml").exists()

    # Per-batch scratch artifacts are cleaned up, not left behind.
    assert not list(outdir.glob("_batch_*"))


def test_fetch_batch_dedups_on_missing_accession(tmp_path, monkeypatch):
    """An accession NCBI doesn't return (invalid/superseded) is warned about,
    not a hard failure -- the other accessions in the same batch still resolve."""
    outdir = tmp_path / "ncbi_cache"
    accessions = ["GCA_OK", "GCA_MISSING"]

    monkeypatch.setattr(fgab.shutil, "which", lambda name: "/usr/bin/fake_datasets")
    monkeypatch.setattr(fgab, "datasets_version", lambda exe: "fake-version 1.0")

    def fake_run(cmd, check=True, **kwargs):
        if cmd[1:3] == ["download", "genome"]:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            _write_fake_batch_zip(zip_path, ["GCA_OK"], include_protein=False)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(fgab.subprocess, "run", fake_run)

    fgab.fetch_batch(accessions, outdir=outdir, include_protein=False, batch_size=75)

    assert (outdir / "GCA_OK" / "extracted" / "ncbi_dataset" / "data" / "GCA_OK").exists()
    assert not (outdir / "GCA_MISSING").exists()
