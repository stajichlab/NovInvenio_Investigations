import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_site import render_markdown, render_run_index  # noqa: E402
from sync_pangenome_report import main  # noqa: E402

REPO_ROOT = Path(__file__).parent.parent

REPORT_MD = """# Pangenome Island + Pfam Enrichment Report

## Pipeline diagnostics

- **rescue_redundancy** [OK]: 1/10 (10.0%) rescuable cells

## Pangenome composition

Total families: 3

![Composition](figures/core_shell_cloud_pie.png)

| Locus | Size | Pfam domains |
|---|---|---|
| s1:c1:1-100 | 5 | <script>alert(1)</script> |
"""


def _fake_study(root: Path, run: str = "run_a", with_optional: bool = True) -> Path:
    study = root / "studies" / "fungi" / "demo_pangenome"
    study.mkdir(parents=True)
    (study / "species.csv").write_text("Short,Group\ns1,IN\ns2,IN\nout,OUT\n")
    pg = study / "results" / run / "output" / "pangenome"
    (pg / "report" / "figures").mkdir(parents=True)
    (pg / "report" / "figures_pdf").mkdir()
    (pg / "report_tables").mkdir()
    (pg / "report" / "report.md").write_text(REPORT_MD)
    (pg / "report" / "figures" / "core_shell_cloud_pie.png").write_bytes(b"\x89PNG fake")
    (pg / "report" / "figures_pdf" / "core_shell_cloud_pie.pdf").write_bytes(b"%PDF fake")
    (pg / "report" / "pangenome_openness.tsv").write_text("metric\tvalue\ngamma\t0.33\n")
    (pg / "frequency_table.tsv").write_text(
        "family\tfrequency\tstrain_count\tbin\nf1\t1\t3\tcore\nf2\t0.6\t2\tshell\nf3\t0.3\t1\tsingleton\n")
    (pg / "report_tables" / "per_strain_summary.tsv").write_text("Short\tn_families\ns1\t3\ns2\t2\nout\t1\n")
    # A large per-pair table that must NOT be published.
    (pg / "pair_classification.tsv").write_text("a\tb\n")
    if with_optional:
        (pg / "island_synteny.html").write_text("<html>synteny</html>")
        (pg / "assembly_quality_report.md").write_text("# Assembly Quality vs Pangenome Content QC\n\nrho table\n")
    return study


def _run(root: Path, *extra: str) -> int:
    return main(["--study", "fungi/demo_pangenome", "--repo-root", str(root), *extra])


def test_stages_run_page_assets_and_study_index(tmp_path):
    _fake_study(tmp_path)
    assert _run(tmp_path, "--run", "run_a") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"

    html = (run_docs / "report.html").read_text()
    assert "Pangenome composition" in html
    assert 'src="figures/core_shell_cloud_pie.png"' in html
    assert "<table>" in html
    assert "archive/frequency_table.tsv.gz" in html
    assert 'href="island_synteny.html"' in html
    assert 'href="assembly_quality.html"' in html

    assert (run_docs / "figures" / "core_shell_cloud_pie.png").read_bytes() == b"\x89PNG fake"
    assert (run_docs / "figures_pdf" / "core_shell_cloud_pie.pdf").exists()
    assert (run_docs / "island_synteny.html").exists()
    assert "rho table" in (run_docs / "assembly_quality.html").read_text()
    with gzip.open(run_docs / "archive" / "pangenome_openness.tsv.gz", "rt") as fh:
        assert fh.read() == "metric\tvalue\ngamma\t0.33\n"
    assert not any("pair_classification" in p.name for p in (run_docs / "archive").iterdir())
    assert "report.html" in (run_docs / "index.html").read_text()

    meta = json.loads((run_docs / "run.json").read_text())
    assert meta["run"] == "run_a"
    assert meta["n_families"] == 3
    assert meta["n_strains"] == 3
    assert meta["source_dir"] == "studies/fungi/demo_pangenome/results/run_a/output/pangenome"
    assert len(meta["report_md_sha256"]) == 64

    index = (tmp_path / "docs" / "fungi" / "demo_pangenome" / "report.html").read_text()
    assert 'href="run_a/report.html"' in index


def test_raw_html_in_report_md_is_escaped(tmp_path):
    _fake_study(tmp_path)
    _run(tmp_path, "--run", "run_a")
    html = (tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "report.html").read_text()
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;" in html


def test_optional_outputs_absent_means_no_dead_links(tmp_path):
    _fake_study(tmp_path, with_optional=False)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    html = (run_docs / "report.html").read_text()
    assert "island_synteny.html" not in html
    assert "assembly_quality.html" not in html
    assert "diagnostics.tsv.gz" not in html


def test_resync_removes_stale_release_files(tmp_path):
    study = _fake_study(tmp_path)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert (run_docs / "island_synteny.html").exists()
    pg = study / "results" / "run_a" / "output" / "pangenome"
    (pg / "island_synteny.html").unlink()
    (pg / "report" / "figures" / "core_shell_cloud_pie.png").unlink()
    _run(tmp_path, "--run", "run_a")
    assert not (run_docs / "island_synteny.html").exists()
    assert not (run_docs / "figures" / "core_shell_cloud_pie.png").exists()


def test_two_runs_both_listed(tmp_path):
    study = _fake_study(tmp_path, run="run_a")
    pg_b = study / "results" / "run_b" / "output"
    import shutil
    shutil.copytree(study / "results" / "run_a" / "output", pg_b)
    _run(tmp_path, "--run", "run_a")
    _run(tmp_path, "--run", "run_b")
    index = (tmp_path / "docs" / "fungi" / "demo_pangenome" / "report.html").read_text()
    assert 'href="run_a/report.html"' in index and 'href="run_b/report.html"' in index


def test_archive_gzip_is_deterministic(tmp_path):
    _fake_study(tmp_path)
    _run(tmp_path, "--run", "run_a")
    gz = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "archive" / "frequency_table.tsv.gz"
    first = gz.read_bytes()
    _run(tmp_path, "--run", "run_a")
    assert gz.read_bytes() == first


def test_ambiguous_pangenome_dir_fails(tmp_path):
    study = _fake_study(tmp_path)
    import shutil
    shutil.copytree(study / "results" / "run_a" / "output", study / "results" / "run_a" / "output2")
    with pytest.raises(SystemExit, match="exactly one"):
        _run(tmp_path, "--run", "run_a")


def test_undefined_study_fails(tmp_path):
    with pytest.raises(SystemExit, match="species.csv"):
        _run(tmp_path, "--run", "run_a")


def test_render_markdown_disables_raw_html():
    assert "<b>" not in render_markdown("<b>x</b>")


def test_render_run_index_empty_and_counts():
    assert "No runs published yet" in render_run_index("fungi", "s", [])
    html = render_run_index("fungi", "s", [{"run": "r", "published": "2026-09-23",
                                            "n_families": 54421, "n_strains": 530}])
    assert "54,421" in html and "530" in html


@pytest.mark.parametrize("path,ignored", [
    ("docs/fungi/demo/run_a/figures/x.png", True),
    ("docs/fungi/demo/run_a/figures_pdf/x.pdf", True),
    ("docs/fungi/demo/run_a/archive/x.tsv.gz", True),
    ("docs/fungi/demo/run_a/island_synteny.html", True),
    ("docs/fungi/demo/run_a/assembly_quality.html", True),
    ("docs/fungi/demo/run_a/report.html", False),
    ("docs/fungi/demo/run_a/run.json", False),
    ("docs/fungi/demo/report.html", False),
])
def test_gitignore_matches_data_classes(path, ignored):
    r = subprocess.run(["git", "check-ignore", "-q", "--no-index", path], cwd=REPO_ROOT)
    assert (r.returncode == 0) == ignored, path
