"""Tests for bin/annotate_proteomes_pfam.py (selection, job sizing, script text, merge, finalize)."""
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import annotate_proteomes_pfam as app  # noqa: E402
from test_hmm_domtbl_to_annotations import HEAD, HITS, search_line, trailer  # noqa: E402

SPECIES = ("Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID\n"
           "Uni1,U one,,IN,X,uniprot,UP000000001,1\n"
           "Ncbi1,N one,,IN,X,ncbi,GCA_000000001.1,2\n"
           "Ncbi2,N two,,OUT,X,ncbi,GCA_000000002.1,3\n"
           "Jgi1,J one,,OUT,X,local,JGI_Abc1,4\n")


def rows():
    import csv
    import io
    return list(csv.DictReader(io.StringIO(SPECIES)))


def test_select_species(tmp_path):
    ann = tmp_path / "annotations"
    ann.mkdir()
    (ann / "GCA_000000002.1.tsv.gz").write_bytes(b"")          # already annotated
    sel, skip = app.select_species(rows(), ann, {})
    assert [s.short for s in sel] == ["Ncbi1", "Jgi1"]
    assert dict(skip)["Uni1"].startswith("UniProt") and "exists" in dict(skip)["Ncbi2"]
    sel, _ = app.select_species(rows(), ann, {}, force=True)
    assert [s.short for s in sel] == ["Ncbi1", "Ncbi2", "Jgi1"]   # --force; UniProt never selected
    sel, _ = app.select_species(rows(), ann, {}, only=["Jgi1"])
    assert [s.short for s in sel] == ["Jgi1"] and sel[0].annotation == ann / "JGI_Abc1.tsv.gz"
    assert sel[0].taxon_id == "4"
    with pytest.raises(SystemExit, match="not in species.csv"):
        app.select_species(rows(), ann, {}, only=["Nope"])


@pytest.mark.parametrize("n,maxn,expect", [
    (0, 35000, []),
    (14612, 35000, [(0, 14612)]),
    (35000, 35000, [(0, 35000)]),
    (35001, 35000, [(0, 17501), (17501, 35001)]),
    (80000, 35000, [(0, 26667), (26667, 53334), (53334, 80000)]),
])
def test_plan_chunks(n, maxn, expect):
    chunks = app.plan_chunks(n, maxn)
    assert chunks == expect
    assert sum(e - s for s, e in chunks) == n
    assert all(e - s <= maxn for s, e in chunks)


def test_time_limit():
    assert app.time_limit(14612) == "02:00:00"      # 24 min estimate -> minimum 2 h
    assert app.time_limit(35000) == "03:00:00"      # 58 min estimate -> 2x + 1 h


def test_chunk_script_text(tmp_path):
    sp = app.Species("Ncbi1", "GCA_1", "ncbi", "2", Path("/abs/pep/Ncbi1.pep.fa"), tmp_path / "a.tsv.gz")
    txt = app.chunk_script(sp, 1, 100, 250, Path("/abs/work"), "/db/Pfam-A.hmm")
    assert "--constraint=milan" in txt and "-t 02:00:00" in txt
    assert "${BASH_SOURCE" not in txt
    assert "${SCRATCH:?" in txt
    assert "awk -v s=100 -v e=250" in txt and "expected 150" in txt
    assert "hmmsearch --cpu" in txt and "--cut_ga" in txt and "/db/Pfam-A.hmm" in txt
    assert "/abs/work/Ncbi1/chunk_001.domtbl.zst" in txt
    assert "#SBATCH -o /abs/work/Ncbi1/logs/pfam_Ncbi1_c001.%j.out" in txt
    fin = app.finalize_script(sp, "fungi/s", Path("/abs/work"), Path("/abs/ann"), "/db/Pfam-A.hmm", force=False)
    assert "--finalize" in fin and "--species Ncbi1" in fin and "--force" not in fin
    assert "--constraint=milan" in fin and "${BASH_SOURCE" not in fin


def test_chunk_awk_selects_records(tmp_path):
    """Run the awk filter of the chunk script on a small FASTA."""
    fa = tmp_path / "p.fa"
    fa.write_text("".join(f">p{i}\nMK\nVV\n" for i in range(5)))
    out = subprocess.run(["awk", "-v", "s=1", "-v", "e=3", "/^>/{n++} n>s && n<=e", str(fa)],
                         capture_output=True, text=True, check=True).stdout
    assert [ln for ln in out.splitlines() if ln.startswith(">")] == [">p1", ">p2"]
    assert out.count("MK") == 2


def test_merge_domtbls():
    t0 = HEAD + "".join(search_line(*h) for h in HITS[:3]) + trailer("hmmsearch")
    t1 = HEAD + search_line(*HITS[3]) + trailer("hmmsearch")
    empty = HEAD + trailer("hmmsearch")
    merged = app.merge_domtbls([t0, empty, t1])
    data = [ln for ln in merged.splitlines() if ln and not ln.startswith("#")]
    assert len(data) == 4 and data[3].startswith("P2.1")
    assert merged.rstrip().endswith("# [ok]") and "# Program:         hmmsearch" in merged
    assert app.merge_domtbls([empty]).count("# Program:") == 1
    with pytest.raises(SystemExit, match="ok"):
        app.merge_domtbls([t0, t1.replace("# [ok]\n", "")])


def _study(tmp_path):
    st = tmp_path / "studies" / "fungi" / "s"
    (st / "data_dir" / "pep").mkdir(parents=True)
    (st / "annotations").mkdir()
    (st / "species.csv").write_text(SPECIES)
    (st / "data_dir" / "pep" / "Ncbi1.pep.fa").write_text(">P1.1 kinase [Org]\nM\n>P0.1 x\nM\n>P2.1 y\nM\n")
    for s in ("Ncbi2", "Jgi1"):
        (st / "data_dir" / "pep" / f"{s}.pep.fa").write_text(">Q1\nM\n")
    return st


def test_write_jobs_and_dry_run(tmp_path):
    st = _study(tmp_path)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--dry-run"]) == 0
    assert not (st / "analysis").exists()                       # dry run writes nothing
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1",
                     "--max-proteins-per-job", "2"]) == 0
    d = st / "analysis" / "pfam_hmmsearch" / "Ncbi1"
    plan = json.loads((d / "plan.json").read_text())
    assert plan["n_proteins"] == 3 and plan["chunks"] == [[0, 2], [2, 3]]
    assert (d / "chunk_000.sbatch").is_file() and (d / "chunk_001.sbatch").is_file()
    assert (d / "finalize.sbatch").is_file() and plan["job_ids"] == []


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not in PATH")
def test_finalize_end_to_end(tmp_path):
    st = _study(tmp_path)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1",
                     "--max-proteins-per-job", "2"]) == 0
    d = st / "analysis" / "pfam_hmmsearch" / "Ncbi1"
    # Fake the two chunk outputs (what the sbatch jobs write).
    for i, text in enumerate([HEAD + "".join(search_line(*h) for h in HITS[:3]) + trailer("hmmsearch"),
                              HEAD + search_line(*HITS[3]) + trailer("hmmsearch")]):
        app.zstd_write(text, d / f"chunk_{i:03d}.domtbl.zst")
    args = ["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1", "--finalize"]
    assert app.main(args) == 0
    ann = st / "annotations" / "GCA_000000001.1.tsv.gz"
    with gzip.open(ann, "rt") as fh:
        lines = fh.read().splitlines()
    assert len(lines) == 4 and lines[1].split("\t")[5] == "PF00017|PF00069"
    assert Path(str(ann) + ".provenance.yaml").is_file()
    assert (st / "analysis" / "pfam_hmmsearch" / "Ncbi1.pfam.hmmsearch.domtbl.zst").is_file()
    before = ann.read_bytes()
    assert app.main(args) == 0                                  # exists: skipped, not rewritten
    assert ann.read_bytes() == before
    # Normal run now skips Ncbi1 (idempotent).
    sel, skip = app.select_species(app.read_csv(st / "species.csv"), st / "annotations", {})
    assert "Ncbi1" not in [s.short for s in sel]


def test_finalize_missing_chunk_stops(tmp_path):
    st = _study(tmp_path)
    app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1"])
    with pytest.raises(SystemExit, match="missing chunk"):
        app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1", "--finalize"])
    assert not (st / "annotations" / "GCA_000000001.1.tsv.gz").exists()
