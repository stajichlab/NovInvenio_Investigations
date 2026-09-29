import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_site import render_markdown  # noqa: E402
from study_runs import render_run_index  # noqa: E402
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


def _fake_study(root: Path, run: str = "run_a", with_optional: bool = True,
                publish: bool = True) -> Path:
    study = root / "studies" / "fungi" / "demo_pangenome"
    study.mkdir(parents=True)
    (study / "species.csv").write_text("Short,Group\ns1,IN\ns2,IN\nout,OUT\n")
    if publish:
        (study / "publish.yaml").write_text("study: demo_pangenome\n")
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
        (pg / "clinker").mkdir()
        (pg / "clinker" / "L001.html").write_text("<html>clinker</html>")
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
    assert "No runs published yet" in render_run_index("fungi", "s", [], None)
    html = render_run_index("fungi", "s", [{"run": "r", "published": "2026-09-23",
                                            "n_families": 54421, "n_strains": 530}], "r")
    assert "54,421" in html and "530" in html


def test_missing_publish_yaml_skips_without_error(tmp_path, capsys):
    _fake_study(tmp_path, publish=False)
    assert _run(tmp_path, "--run", "run_a") == 0
    assert not (tmp_path / "docs").exists()
    assert "publish.yaml not found" in capsys.readouterr().err


def test_publish_yaml_study_name_sets_docs_dir(tmp_path):
    study = _fake_study(tmp_path)
    (study / "publish.yaml").write_text("study: cocci\n")
    _run(tmp_path, "--run", "run_a")
    assert (tmp_path / "docs" / "fungi" / "cocci" / "run_a" / "report.html").exists()


def test_single_run_becomes_current_and_index_redirects(tmp_path):
    _fake_study(tmp_path)
    _run(tmp_path, "--run", "run_a")
    study_docs = tmp_path / "docs" / "fungi" / "demo_pangenome"
    assert json.loads((study_docs / "study.json").read_text()) == {"current": "run_a"}
    assert 'url=run_a/report.html' in (study_docs / "index.html").read_text()


def test_bad_run_name_rejected(tmp_path):
    _fake_study(tmp_path)
    with pytest.raises(SystemExit, match="not a valid name"):
        _run(tmp_path, "--run", "bad name")


@pytest.mark.parametrize("path,ignored", [
    ("docs/fungi/demo/run_a/figures/x.png", True),
    ("docs/fungi/demo/run_a/figures_pdf/x.pdf", True),
    ("docs/fungi/demo/run_a/archive/x.tsv.gz", True),
    ("docs/fungi/demo/run_a/island_synteny.html", True),
    ("docs/fungi/demo/run_a/clinker/L001.html", True),
    ("docs/fungi/demo/run_a/assembly_quality.html", True),
    ("docs/fungi/demo/run_a/report.html", False),
    ("docs/fungi/demo/run_a/run.json", False),
    ("docs/fungi/demo/report.html", False),
])
def test_gitignore_matches_data_classes(path, ignored):
    r = subprocess.run(["git", "check-ignore", "-q", "--no-index", path], cwd=REPO_ROOT)
    assert (r.returncode == 0) == ignored, path


def test_run_json_gets_pipeline_commit_from_ni_run_record(tmp_path):
    study = _fake_study(tmp_path)
    launch = study / ".nf_launch" / "run_a"
    launch.mkdir(parents=True)
    (launch / "ni_run.json").write_text(json.dumps(
        {"pipeline": "stajichlab/NovInvenio", "pipeline_commit": "a" * 40}))
    assert _run(tmp_path, "--run", "run_a") == 0
    meta = json.loads((tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "run.json").read_text())
    assert meta["pipeline"] == "pangenome.nf"
    assert meta["pipeline_repo"] == "stajichlab/NovInvenio"
    assert meta["pipeline_commit"] == "a" * 40


def test_run_json_without_ni_run_record_is_unchanged(tmp_path):
    _fake_study(tmp_path)
    assert _run(tmp_path, "--run", "run_a") == 0
    meta = json.loads((tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "run.json").read_text())
    assert "pipeline_commit" not in meta and "pipeline_repo" not in meta


def test_clinker_pages_are_staged_next_to_the_synteny_viewer(tmp_path):
    _fake_study(tmp_path)
    assert _run(tmp_path, "--run", "run_a") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert (run_docs / "clinker" / "L001.html").read_text() == "<html>clinker</html>"


def test_resync_replaces_the_clinker_dir(tmp_path):
    study = _fake_study(tmp_path)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    (run_docs / "clinker" / "L999.html").write_text("stale")
    pg = study / "results" / "run_a" / "output" / "pangenome"
    (pg / "clinker" / "L001.html").unlink()
    _run(tmp_path, "--run", "run_a")
    assert not (run_docs / "clinker" / "L999.html").exists()
    assert not (run_docs / "clinker" / "L001.html").exists()


def test_publish_script_and_pages_deploy_ship_clinker():
    publish = (REPO_ROOT / "bin" / "publish_report_release.sh").read_text()
    static = (REPO_ROOT / ".github" / "workflows" / "static.yml").read_text()
    assert "archive clinker island_synteny.html" in publish
    assert "for d in figures figures_pdf archive clinker; do" in static


# ---- N4: NII publishes clinker pages for the top 25 loci only ----

def _fake_study_with_n_clinker_pages(root: Path, n: int, run: str = "run_a") -> Path:
    study = _fake_study(root, run=run)
    pg = study / "results" / run / "output" / "pangenome"
    clinker = pg / "clinker"
    for f in clinker.iterdir():
        f.unlink()
    for i in range(1, n + 1):
        (clinker / f"L{i:03d}.html").write_text(f"<html>clinker {i}</html>")
    # a non-locus file must be ignored, not counted or copied.
    (clinker / "assets.css").write_text("body{}")
    return study


def test_default_publishes_only_the_top_25_of_30_clinker_pages(tmp_path):
    _fake_study_with_n_clinker_pages(tmp_path, 30)
    assert _run(tmp_path, "--run", "run_a") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    staged = sorted(p.name for p in (run_docs / "clinker").glob("L*.html"))
    assert staged == [f"L{i:03d}.html" for i in range(1, 26)]
    assert not (run_docs / "clinker" / "assets.css").exists()


def test_clinker_publish_top_minus_one_publishes_all_and_injects_nothing(tmp_path):
    study = _fake_study_with_n_clinker_pages(tmp_path, 30)
    assert _run(tmp_path, "--run", "run_a", "--clinker_publish_top", "-1") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    staged = sorted(p.name for p in (run_docs / "clinker").glob("L*.html"))
    assert staged == [f"L{i:03d}.html" for i in range(1, 31)]
    assert "CLINKER_PUBLISHED" not in (run_docs / "island_synteny.html").read_text()
    meta = json.loads((run_docs / "run.json").read_text())
    assert meta["clinker_published"] == 30
    assert meta["clinker_total"] == 30
    assert meta["clinker_full_dir"] == str((study / "results" / "run_a" / "output" /
                                            "pangenome" / "clinker").resolve())


def test_clinker_publish_top_zero_publishes_none(tmp_path):
    _fake_study_with_n_clinker_pages(tmp_path, 30)
    assert _run(tmp_path, "--run", "run_a", "--clinker_publish_top", "0") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert not any((run_docs / "clinker").glob("L*.html")) if (run_docs / "clinker").exists() else True
    meta = json.loads((run_docs / "run.json").read_text())
    assert meta["clinker_published"] == 0
    assert meta["clinker_total"] == 30


def test_run_json_records_clinker_published_and_total(tmp_path):
    _fake_study_with_n_clinker_pages(tmp_path, 30)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    meta = json.loads((run_docs / "run.json").read_text())
    assert meta["clinker_published"] == 25
    assert meta["clinker_total"] == 30
    assert meta["clinker_full_dir"].endswith("/clinker")


def test_island_synteny_gets_the_clinker_published_script_injected_once(tmp_path):
    study = _fake_study_with_n_clinker_pages(tmp_path, 30)
    pg = study / "results" / "run_a" / "output" / "pangenome"
    (pg / "island_synteny.html").write_text("<html><head></head><body>synteny</body></html>")
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    html = (run_docs / "island_synteny.html").read_text()
    assert html.count("window.CLINKER_PUBLISHED=") == 1
    assert html.index("<script>window.CLINKER_PUBLISHED=") < html.index("</head>")
    keys = json.loads(html.split("window.CLINKER_PUBLISHED=", 1)[1].split(";", 1)[0])
    assert keys == [f"L{i:03d}" for i in range(1, 26)]

    # a second sync must not duplicate the injected script.
    _run(tmp_path, "--run", "run_a")
    html2 = (run_docs / "island_synteny.html").read_text()
    assert html2.count("window.CLINKER_PUBLISHED=") == 1


def test_all_files_published_means_no_injection(tmp_path):
    _fake_study_with_n_clinker_pages(tmp_path, 10)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    html = (run_docs / "island_synteny.html").read_text()
    assert "CLINKER_PUBLISHED" not in html
    meta = json.loads((run_docs / "run.json").read_text())
    assert meta["clinker_published"] == 10
    assert meta["clinker_total"] == 10


def test_no_clinker_dir_means_no_clinker_fields_in_run_json(tmp_path):
    _fake_study(tmp_path, with_optional=False)
    _run(tmp_path, "--run", "run_a")
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    meta = json.loads((run_docs / "run.json").read_text())
    assert "clinker_published" not in meta
    assert "clinker_total" not in meta
    assert "clinker_full_dir" not in meta


def test_rerender_record_goes_into_run_json(tmp_path):
    # NovInvenio #212 rollout: report steps re-run offline at a newer commit
    # than the pipeline run; run.json must say so.
    study = _fake_study(tmp_path)
    rec = {"report_steps_commit": "8c55750", "steps": ["FREQUENCY_BINS", "REPORT_RENDER"],
           "rerendered": "2026-09-28"}
    (study / "results" / "run_a" / "output" / "pangenome" / "report" / "rerender.json").write_text(
        json.dumps(rec))
    assert _run(tmp_path, "--run", "run_a") == 0
    meta = json.loads((tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "run.json").read_text())
    assert meta["report_rerender"] == rec


def test_no_rerender_record_no_key(tmp_path):
    _fake_study(tmp_path)
    assert _run(tmp_path, "--run", "run_a") == 0
    meta = json.loads((tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "run.json").read_text())
    assert "report_rerender" not in meta


def test_group_class_overlap_table_archived(tmp_path):
    # NovInvenio #212 PR 2 table.
    study = _fake_study(tmp_path)
    pg = study / "results" / "run_a" / "output" / "pangenome"
    (pg / "report_tables" / "group_class_overlap.tsv").write_text(
        "ingroup_class\toutgroup_class\tn_families\ncore\tcore\t1\n")
    assert _run(tmp_path, "--run", "run_a") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert (run_docs / "archive" / "group_class_overlap.tsv.gz").exists()
    assert "archive/group_class_overlap.tsv.gz" in (run_docs / "report.html").read_text()


def test_empty_tables_are_not_listed_as_downloads(tmp_path):
    # A 0-byte table means "not computed" (e.g. outgroup not binned); no dead link.
    study = _fake_study(tmp_path)
    pg = study / "results" / "run_a" / "output" / "pangenome"
    (pg / "report_tables" / "group_class_overlap.tsv").write_text("")
    assert _run(tmp_path, "--run", "run_a") == 0
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert not (run_docs / "archive" / "group_class_overlap.tsv.gz").exists()
    assert "group_class_overlap" not in (run_docs / "report.html").read_text()
