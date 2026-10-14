"""Tests for bin/annotate_proteomes_pfam.py (selection, names, guard, job sizing, batches,
script text, merge, finalize, --fill-empty-uniprot)."""
import csv
import gzip
import io
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import annotate_proteomes_pfam as app  # noqa: E402
from test_hmm_domtbl_to_annotations import HEAD, HITS, search_line, trailer  # noqa: E402

SPECIES = ("Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID\n"
           "Uni1,U one,,IN,X,uniprot,UP000000001,1\n"
           "Ncbi1,N one,,IN,X,ncbi,GCA_000000001.1,2\n"
           "Ncbi2,N two,,OUT,X,ncbi,GCA_000000002.1,3\n"
           "Jgi1,J one,,OUT,X,local,JGI_Abc1,4\n")
# Local FASTA paths as Protein_Accession (fusarium_FOXY, UHM_* studies).
SPECIES_PATHS = ("Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID\n"
                 "Loc1,L one,,IN,X,local_faa,/bigdata/stajichlab/shared/x/GCA_032206225.1_A_protein.faa.gz,5\n"
                 "Loc2,L two,,OUT,X,local_faa,/bigdata/stajichlab/jpere468/k/GCF_050916275.1.faa,\n")
ANN_COLS = ["accession", "taxon_id", "gene_name", "description", "go_ids", "pfam_ids", "pfam_names",
            "interpro_ids", "ec_numbers", "alphafold_id", "xrefs"]


def rows(text=SPECIES):
    return list(csv.DictReader(io.StringIO(text)))


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
    assert [s.short for s in sel] == ["Jgi1"] and sel[0].annotation == (ann / "JGI_Abc1.tsv.gz").resolve()
    assert sel[0].taxon_id == "4"
    with pytest.raises(SystemExit, match="not in species.csv"):
        app.select_species(rows(), ann, {}, only=["Nope"])


def test_select_species_path_accession_stays_in_annotations(tmp_path):
    """Regression: an absolute Protein_Accession must not put the output beside the source file."""
    ann = tmp_path / "annotations"
    ann.mkdir()
    sel, _ = app.select_species(rows(SPECIES_PATHS), ann, {})
    assert [s.annotation for s in sel] == [(ann / "GCA_032206225.1_A_protein.tsv.gz").resolve(),
                                           (ann / "GCF_050916275.1.tsv.gz").resolve()]
    assert all(s.annotation.parent == ann.resolve() for s in sel)


def test_select_fill_empty_uniprot(tmp_path):
    ann = tmp_path / "annotations"
    ann.mkdir()
    sel, skip = app.select_species(rows(), ann, {}, fill_empty_uniprot=True, only=["Uni1"])
    assert sel == [] and "without" in dict(skip)["Uni1"]
    _write_ann(ann / "UP000000001.tsv.gz", [{"accession": "P1"}])        # no terms
    sel, _ = app.select_species(rows(), ann, {}, fill_empty_uniprot=True, only=["Uni1"])
    assert [(s.short, s.mode) for s in sel] == [("Uni1", "fill_uniprot")]
    _write_ann(ann / "UP000000001.tsv.gz", [{"accession": "P1", "go_ids": "GO:0000001:IEA"}])
    sel, skip = app.select_species(rows(), ann, {}, fill_empty_uniprot=True, only=["Uni1"])
    assert sel == [] and "has GO/Pfam/InterPro" in dict(skip)["Uni1"]


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
    assert app.time_limit(45000) == "04:00:00"      # 75 min estimate (one full batch)


def test_batch_capacity():
    assert app.batch_capacity(75, 0.1) == 45000
    assert app.batch_capacity(60, 0.1) == 36000


def test_pack_batches_small_proteomes():
    tasks = [(f"S{i:02d}", 0, 5000 + 10 * i) for i in range(34)]       # 34 small proteomes
    batches = app.pack_batches(tasks, 45000)
    assert sorted(t for b in batches for t in b) == sorted(tasks)      # every task once
    assert all(sum(t[2] for t in b) <= 45000 for b in batches)
    # Items of 5,000-5,330 proteins: at most 8 fit under 45,000, so 34 items need ceil(34/8) = 5.
    assert len(batches) == 5 and [len(b) for b in batches] == [8, 8, 8, 8, 2]
    assert app.pack_batches(tasks, 45000) == batches                   # deterministic


def test_pack_batches_large_and_oversize():
    tasks = [("Big", 0, 30000), ("Big", 1, 30000), ("Small", 0, 14000), ("Huge", 0, 50000)]
    batches = app.pack_batches(tasks, 45000)
    assert batches[0] == [("Huge", 0, 50000)]                          # larger than cap: alone
    loads = sorted(sum(t[2] for t in b) for b in batches)
    assert loads == [30000, 44000, 50000]
    assert app.pack_batches([], 45000) == []


def test_chunk_and_batch_script_text(tmp_path):
    sp = app.Species("Ncbi1", "GCA_1", "ncbi", "2", Path("/abs/pep/Ncbi1.pep.fa"), tmp_path / "a.tsv.gz")
    txt = app.chunk_script(sp, 1, 100, 250, Path("/abs/work"), "/db/Pfam-A.hmm")
    assert "--constraint=milan" in txt and "-t 02:00:00" in txt
    assert "${BASH_SOURCE" not in txt
    assert "${SCRATCH:?" in txt
    assert "awk -v s=100 -v e=250" in txt and "expected 150" in txt
    assert "hmmsearch --cpu" in txt and "--cut_ga" in txt and "/db/Pfam-A.hmm" in txt
    assert "/abs/work/Ncbi1/chunk_001.domtbl.zst" in txt
    b = app.batch_script(3, [("Ncbi1", 0, 20000, Path("/abs/work/Ncbi1/chunk_000.sbatch")),
                             ("Jgi1", 0, 15000, Path("/abs/work/Jgi1/chunk_000.sbatch"))],
                         Path("/abs/work/batches/t"))
    assert "--constraint=milan" in b and "-t 03:00:00" in b and "#SBATCH -J pfam_b003" in b
    assert "#SBATCH -o /abs/work/batches/t/logs/pfam_b003.%j.out" in b
    assert b.index("/abs/work/Ncbi1/chunk_000.sbatch") < b.index("/abs/work/Jgi1/chunk_000.sbatch")
    assert 'bash "$s"' in b and "failed+=" in b and "exit 1" in b
    assert "${BASH_SOURCE" not in b and "${SCRATCH:?" in b
    fin = app.finalize_script(["Ncbi1", "Jgi1"], "fungi/s", Path("/abs/work"), Path("/abs/ann"),
                              "/db/Pfam-A.hmm", force=False, batch_dir=Path("/abs/work/batches/t"))
    assert "--finalize" in fin and "--species Ncbi1 Jgi1" in fin and "--force" not in fin
    assert "--fill-empty-uniprot" not in fin
    assert "--constraint=milan" in fin and "${BASH_SOURCE" not in fin


def test_chunk_awk_selects_records(tmp_path):
    """Run the awk filter of the chunk script on a small FASTA."""
    fa = tmp_path / "p.fa"
    fa.write_text("".join(f">p{i}\nMK\nVV\n" for i in range(5)))
    out = subprocess.run(["awk", "-v", "s=1", "-v", "e=3", "/^>/{n++} n>s && n<=e", str(fa)],
                         capture_output=True, text=True, check=True).stdout
    assert [ln for ln in out.splitlines() if ln.startswith(">")] == [">p1", ">p2"]
    assert out.count("MK") == 2


def test_batch_script_runs_all_and_reports_failure(tmp_path):
    """The batch runs every chunk script; one failure does not stop the next one."""
    ok1, bad, ok2 = tmp_path / "a.sh", tmp_path / "b.sh", tmp_path / "c.sh"
    for p, body in ((ok1, f"touch {tmp_path}/a.done"), (bad, "exit 3"), (ok2, f"touch {tmp_path}/c.done")):
        p.write_text(f"#!/usr/bin/bash\nset -euo pipefail\n{body}\n")
    txt = app.batch_script(0, [("A", 0, 1, ok1), ("B", 0, 1, bad), ("C", 0, 1, ok2)], tmp_path)
    (tmp_path / "batch.sh").write_text(txt)
    res = subprocess.run(["bash", str(tmp_path / "batch.sh")], capture_output=True, text=True,
                         env={"SCRATCH": str(tmp_path), "PATH": "/usr/bin:/bin"})
    assert res.returncode == 1 and "1 chunk(s) failed" in res.stderr
    assert (tmp_path / "a.done").exists() and (tmp_path / "c.done").exists()


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


def test_fill_rows():
    rs = [{"accession": "P1", "pfam_ids": "", "pfam_names": "", "x": "keep"},
          {"accession": "P2", "pfam_ids": "", "pfam_names": "", "x": "keep"}]
    c = app.fill_rows(rs, {"P1": ("PF1|PF2", "A|B"), "P2": ("", "")})
    assert c == {"n_rows": 2, "n_rows_filled": 1}
    assert rs[0]["pfam_ids"] == "PF1|PF2" and rs[0]["pfam_names"] == "A|B" and rs[0]["x"] == "keep"
    assert rs[1]["pfam_ids"] == ""


def _write_ann(path, rs):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=ANN_COLS, delimiter="\t", lineterminator="\n", restval="")
    w.writeheader()
    w.writerows(rs)
    with gzip.open(path, "wt") as fh:
        fh.write(buf.getvalue())


def _study(tmp_path, species=SPECIES):
    st = tmp_path / "studies" / "fungi" / "s"
    (st / "data_dir" / "pep").mkdir(parents=True)
    (st / "annotations").mkdir()
    (st / "species.csv").write_text(species)
    (st / "data_dir" / "pep" / "Ncbi1.pep.fa").write_text(">P1.1 kinase [Org]\nM\n>P0.1 x\nM\n>P2.1 y\nM\n")
    for s in ("Ncbi2", "Jgi1", "Loc1", "Loc2"):
        (st / "data_dir" / "pep" / f"{s}.pep.fa").write_text(">Q1\nM\n")
    return st


def _batch_dir(st):
    dirs = sorted((st / "analysis" / "pfam_hmmsearch" / "batches").iterdir())
    assert len(dirs) == 1
    return dirs[0]


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
    bd = _batch_dir(st)
    assert (bd / "finalize.sbatch").is_file() and plan["job_ids"] == []
    assert plan["finalize_script"] == str(bd / "finalize.sbatch")
    batches = json.loads((bd / "batches.json").read_text())
    assert len(batches["batches"]) == 1                          # 3 proteins: one batch job
    assert batches["batches"][0]["chunks"] == [["Ncbi1", 0, 2], ["Ncbi1", 1, 1]]


def test_batches_pack_several_species(tmp_path):
    st = _study(tmp_path)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s"]) == 0
    bd = _batch_dir(st)
    batches = json.loads((bd / "batches.json").read_text())
    assert len(batches["batches"]) == 1
    assert sorted(c[0] for c in batches["batches"][0]["chunks"]) == ["Jgi1", "Ncbi1", "Ncbi2"]
    fin = (bd / "finalize.sbatch").read_text()
    assert "--species Ncbi1 Ncbi2 Jgi1" in fin


def test_path_accessions_dry_run_and_write(tmp_path, capsys):
    st = _study(tmp_path, SPECIES_PATHS)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--dry-run"]) == 0
    err = capsys.readouterr().err
    assert f"output {(st / 'annotations' / 'GCF_050916275.1.tsv.gz').resolve()}" in err
    outs = [ln.split("; output ", 1)[1] for ln in err.splitlines() if "; output " in ln]
    assert len(outs) == 2
    assert all(Path(p).parent == (st / "annotations").resolve() for p in outs)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s"]) == 0
    plan = json.loads((st / "analysis" / "pfam_hmmsearch" / "Loc2" / "plan.json").read_text())
    assert plan["annotation"] == str((st / "annotations" / "GCF_050916275.1.tsv.gz").resolve())


def test_guard_refuses_outside(tmp_path):
    g = app.Guard([tmp_path / "annotations", tmp_path / "work"])
    g(tmp_path / "work" / "x" / "plan.json")
    with pytest.raises(ValueError, match="outside"):
        g(Path("/bigdata/stajichlab/shared/projects/x.tsv.gz"))
    with pytest.raises(ValueError, match="outside"):
        app.write_text_guarded(g, tmp_path / "elsewhere.txt", "x")
    assert not (tmp_path / "elsewhere.txt").exists()


def _fake_chunks(st, short, texts):
    d = st / "analysis" / "pfam_hmmsearch" / short
    for i, text in enumerate(texts):
        app.zstd_write(text, d / f"chunk_{i:03d}.domtbl.zst")


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not in PATH")
def test_finalize_end_to_end(tmp_path):
    st = _study(tmp_path)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1",
                     "--max-proteins-per-job", "2"]) == 0
    # Fake the two chunk outputs (what the batch job writes).
    _fake_chunks(st, "Ncbi1", [HEAD + "".join(search_line(*h) for h in HITS[:3]) + trailer("hmmsearch"),
                               HEAD + search_line(*HITS[3]) + trailer("hmmsearch")])
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


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not in PATH")
def test_finalize_path_accession_writes_in_annotations(tmp_path):
    st = _study(tmp_path, SPECIES_PATHS)
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Loc2"]) == 0
    _fake_chunks(st, "Loc2", [HEAD + trailer("hmmsearch")])
    assert app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Loc2", "--finalize"]) == 0
    files = sorted(p.name for p in (st / "annotations").iterdir())
    assert files == ["GCF_050916275.1.tsv.gz", "GCF_050916275.1.tsv.gz.provenance.yaml"]


def test_finalize_missing_chunk_stops(tmp_path):
    st = _study(tmp_path)
    app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1"])
    with pytest.raises(SystemExit, match="missing chunk"):
        app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1", "--finalize"])
    assert not (st / "annotations" / "GCA_000000001.1.tsv.gz").exists()


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not in PATH")
def test_finalize_continues_after_missing_species(tmp_path):
    """One species without chunk tables does not stop the others; the run still fails."""
    st = _study(tmp_path)
    app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1", "Jgi1"])
    _fake_chunks(st, "Jgi1", [HEAD + trailer("hmmsearch")])
    with pytest.raises(SystemExit, match="1 species"):
        app.main(["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Ncbi1", "Jgi1",
                  "--finalize"])
    assert (st / "annotations" / "JGI_Abc1.tsv.gz").is_file()
    assert not (st / "annotations" / "GCA_000000001.1.tsv.gz").exists()


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not in PATH")
def test_fill_empty_uniprot_end_to_end(tmp_path):
    st = _study(tmp_path)
    (st / "data_dir" / "pep" / "Uni1.pep.fa").write_text(
        ">tr|P1|P1_X kinase OS=x\nM\n>tr|P2|P2_X other OS=x\nM\n>sp|P4|P4_X no row OS=x\nM\n")
    ann = st / "annotations" / "UP000000001.tsv.gz"
    orig_rows = [{"accession": "P1", "taxon_id": "1", "gene_name": "g1", "description": "Kinase"},
                 {"accession": "P2", "taxon_id": "1", "gene_name": "g2", "description": "Other"},
                 {"accession": "P3", "taxon_id": "1", "gene_name": "g3", "description": "Not in FASTA"}]
    _write_ann(ann, orig_rows)
    orig_bytes = ann.read_bytes()
    base = ["--repo-root", str(tmp_path), "--study", "fungi/s", "--species", "Uni1"]
    assert app.main(base) == 0                                   # without the flag: nothing to do
    assert not (st / "analysis").exists()
    assert app.main(base + ["--fill-empty-uniprot"]) == 0
    fin = (_batch_dir(st) / "finalize.sbatch").read_text()
    assert "--fill-empty-uniprot" in fin
    hits = [("tr|P1|P1_X", "Kinase", "PF00069.32", "1e-50", "180.5", 40),
            ("sp|P4|P4_X", "SH2", "PF00017.30", "2e-10", "35.0", 5)]
    _fake_chunks(st, "Uni1", [HEAD + "".join(search_line(*h) for h in hits) + trailer("hmmsearch")])
    assert app.main(base + ["--finalize", "--fill-empty-uniprot"]) == 0
    with gzip.open(ann, "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == ANN_COLS                        # same columns
        got = list(rdr)
    assert [r["accession"] for r in got] == ["P1", "P2", "P3"]   # same rows, same order
    assert got[0]["pfam_ids"] == "PF00069" and got[0]["pfam_names"] == "Kinase"
    assert got[0]["gene_name"] == "g1" and got[0]["description"] == "Kinase"
    assert got[1]["pfam_ids"] == "" and got[2]["pfam_ids"] == ""
    orig = st / "annotations_orig" / "UP000000001.tsv.gz"
    assert orig.read_bytes() == orig_bytes                       # original kept, outside annotations/
    rec = yaml.safe_load(Path(str(ann) + ".provenance.yaml").read_text())
    assert rec["counts"]["n_rows_filled"] == 1
    assert rec["counts"]["n_proteins_with_pfam_but_no_row"] == 1          # P4
    assert rec["counts"]["n_rows_not_in_fasta"] == 1                      # P3
    assert rec["vocabularies_not_available"] == ["go", "interpro"]
    assert rec["original_file"]["path"] == str(orig.resolve())
    # Second run: the file now has terms -> not selected, not changed.
    filled = ann.read_bytes()
    assert app.main(base + ["--finalize", "--fill-empty-uniprot"]) == 0
    assert ann.read_bytes() == filled
