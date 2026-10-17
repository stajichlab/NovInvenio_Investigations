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
    assert "archive clinker locus_summary island_synteny.html" in publish
    assert "locus_summary.html" in publish
    assert "for d in figures figures_pdf archive clinker locus_summary starship_overlap; do" in static
    assert "assembly_quality.html locus_summary.html starship_overlap.html; do" in static


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


def test_stage_locus_summary_publishes_subset(tmp_path):
    import sync_pangenome_report as sp
    fig = tmp_path / "fig"
    (fig / "locus_strips").mkdir(parents=True)
    (fig / "locus_table.tsv").write_text(
        "key\ttype\tsize\tn_strains\tcontent_n\tcontent_frac\tn_states\ttier\tstarship\tpair\n"
        "L001\tswap_candidate\t7\t293\t272\t0.928\t19\tfull\t\tA(1) vs B(2)\n"
        "L002\tunresolved\t5\t100\t3\t0.03\t4\tcontig_end\tLamia-h4\t\n")
    (fig / "locus_overview.png").write_bytes(b"png")
    (fig / "locus_strips" / "L001.png").write_bytes(b"png")
    run = tmp_path / "docs" / "set" / "run"
    run.mkdir(parents=True)
    assert sp.stage_locus_summary(fig, run)
    assert (run / "locus_summary" / "locus_overview.png").is_file()
    assert not (run / "locus_summary" / "locus_strips").exists()
    page = (run / "locus_summary.html").read_text()
    assert "swap_candidate 1" in page and "Lamia-h4" in page
    assert not sp.stage_locus_summary(tmp_path / "missing", run)


# ---- Starships, islands and SM clusters page (--starship_overlap_dir) ----

_SSO_SUMMARY = (
    "run\tmetric\tvalue\tdefinition\n"
    "r1\tn_families\t100\trows\n"
    "r1\tn_islands\t20\trows\n"
    "r1\tannotated_strains\tAf293\tstrains with called elements\n"
    "r1\tn_elements\t2\t\n"
    "r1\tn_elements_high\t1\t\n"
    "r1\tcaptain_hits\t7\ttblout targets with full-sequence E <= 1e-05\n"
    "r1\tcaptain_hits_unmapped\t0\t\n"
    "r1\tn_starship_all_families\t5\t\n"
    "r1\tn_starship_all_families_in_island\t4\t\n"
    "r1\tn_islands_with_starship_all_family\t3\t\n"
    "r1\tn_captain_families\t6\t\n"
    "r1\tn_captain_families_in_island\t2\t\n"
    "r1\tn_islands_with_captain_family\t2\t\n")
_SSO_ENR = (
    "annotation\tuniverse\tN_universe\tK_annotated\tn_island\tk_both\texpected_k\tfold\todds_ratio\t"
    "p_hypergeom_greater\tfrac_annotated_in_island\tfrac_universe_in_island\n"
    "starship_all\tnoncore\t80\t5\t30\t4\t1.88\t2.13\tinf\t1.23e-02\t0.8\t0.375\n"
    "captain\tnoncore\t80\t6\t30\t2\t2.25\t0.89\t0.8\t6.54e-01\t0.3333\t0.375\n")


def _sso_dir(root: Path, with_sm: bool) -> Path:
    d = root / "analysis" / "starship_island_overlap" / "r1"
    d.mkdir(parents=True)
    (d / "summary.tsv").write_text(_SSO_SUMMARY)
    (d / "enrichment.tsv").write_text(_SSO_ENR)
    (d / "per_starship.tsv").write_text("element_id\tlabel\tconfidence\tn_families\nE1\tLamia h4\thigh\t3\n")
    (d / "PROVENANCE.yaml").write_text("derived_by: bin/starship_island_overlap.py\n")
    (d / "per_island.tsv.gz").write_bytes(gzip.compress(b"island_idx\n0\n"))
    figs = d.parent / "figures"
    figs.mkdir()
    (figs / "r1_starship_islands.png").write_bytes(b"png")
    (figs / "r1_extra_starship_islands.png").write_bytes(b"other run")
    if with_sm:
        (d / "sm_by_type.tsv").write_text(
            "species\tsm_type\tn_clusters\tn_with_captain_gene\tpct_with_captain_gene\n"
            "all\tT1PKS\t12\t3\t25.00\n")
        (d / "sm_starship_overlap.tsv").write_text(
            "metric\tvalue\tdefinition\nsm_clusters\t12\tantiSMASH 8.0.1 regions\n"
            "sm_clusters_overlapping_starship_ge1gene\tNA\tno Starship element calls\n"
            "sm_clusters_overlapping_starship_ge50pct\tNA\tno Starship element calls\n")
        (d / "PROVENANCE_sm.yaml").write_text("source_annotation: own funannotate (not the freeze)\n")
        (figs / "r1_sm.png").write_bytes(b"png")
    return d


def test_stage_starship_overlap_with_sm_tables(tmp_path):
    import re as _re
    import sync_pangenome_report as sp
    d = _sso_dir(tmp_path, with_sm=True)
    run = tmp_path / "docs" / "set" / "run"
    run.mkdir(parents=True)
    assert sp.stage_starship_overlap(d, run)
    dest = run / "starship_overlap"
    for f in ("summary.tsv", "enrichment.tsv", "per_starship.tsv", "per_island.tsv.gz",
              "sm_by_type.tsv", "PROVENANCE.yaml", "PROVENANCE_sm.yaml", "starship_islands.png", "sm.png"):
        assert (dest / f).is_file(), f
    assert (dest / "starship_islands.png").read_bytes() == b"png"   # exact-name figure, not r1_extra_
    page = (run / "starship_overlap.html").read_text()
    assert "Starship families, all elements" in page and "1.23e-02" in page
    assert "4 of 5 (0.8)" in page and "30 of 80 (0.375)" in page and "3 of 20" in page
    assert "Lamia h4" in page and "Captain-gene (DUF3435) families" in page
    assert "antiSMASH 8.0.1 regions" in page and "T1PKS" in page
    assert "no Starship element calls" in page          # >= 1 gene / >= 50% rule lines
    assert "own funannotate (not the freeze)" in page
    assert "Data and limits" in page and "exist only for Af293" in page
    assert "<script src" not in page and "http" not in page.split("<footer>")[0]
    for href in _re.findall(r'(?:href|src)="(starship_overlap/[^"]+)"', page):
        assert (run / href).is_file(), href


def test_stage_starship_overlap_without_sm_or_optional_tables(tmp_path):
    import sync_pangenome_report as sp
    d = _sso_dir(tmp_path, with_sm=False)
    for f in ("enrichment.tsv", "per_starship.tsv"):
        (d / f).unlink()
    (d / "summary.tsv").write_text("\n".join(l for l in _SSO_SUMMARY.splitlines()
                                             if "starship" not in l) + "\n")
    run = tmp_path / "run"
    run.mkdir()
    assert sp.stage_starship_overlap(d, run)
    page = (run / "starship_overlap.html").read_text()
    assert "no Starship element coordinates" in page
    assert "Not available: no sm_by_type.tsv" in page
    assert "enrichment.tsv: not available" in page
    assert not (run / "starship_overlap" / "sm.png").exists()


def test_stage_starship_overlap_needs_summary(tmp_path):
    import sync_pangenome_report as sp
    run = tmp_path / "run"
    run.mkdir()
    assert not sp.stage_starship_overlap(tmp_path / "missing", run)
    assert not (run / "starship_overlap.html").exists()


def test_starship_overlap_page_is_linked_and_removed_on_resync(tmp_path):
    _fake_study(tmp_path)
    d = _sso_dir(tmp_path, with_sm=True)
    run_docs = tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a"
    assert _run(tmp_path, "--run", "run_a", "--starship_overlap_dir", str(d)) == 0
    assert 'href="starship_overlap.html"' in (run_docs / "report.html").read_text()
    assert (run_docs / "starship_overlap" / "summary.tsv").is_file()
    assert _run(tmp_path, "--run", "run_a") == 0     # no option: page and dir are gone
    assert "starship_overlap" not in (run_docs / "report.html").read_text()
    assert not (run_docs / "starship_overlap").exists()
    assert not (run_docs / "starship_overlap.html").exists()


def test_publish_script_ships_starship_overlap():
    publish = (REPO_ROOT / "bin" / "publish_report_release.sh").read_text()
    gitignore = (REPO_ROOT / ".gitignore").read_text()
    assert "starship_overlap starship_overlap.html; do" in publish
    assert "docs/*/*/*/starship_overlap/\n" in gitignore
    assert "docs/*/*/*/starship_overlap.html\n" in gitignore


# ---- separate SM dir (Af293 antiSMASH, full_v072_antismash/ schema of 2026-10-17) ----
# Rows copied from studies/fungi/Afumigatus_pangenome/analysis/starship_island_overlap/
# full_v072_antismash/ (summary.tsv, sm_by_type.tsv, sm_starship_overlap.tsv,
# sm_permutation.tsv, sm_family_class_direction.tsv).

_AS_SUMMARY = (
    "metric\tvalue\tsource\n"
    "antismash_version\t8.0.4\tAf293.json.gz\n"
    "af293_cds_in_antismash\t9630\tCDS features in antiSMASH JSON (from GFF3)\n"
    "af293_cds_mapped_to_member\t9627\tcoord:9623; fasta_GN:4; none:3\n"
    "sm_clusters\t44\tantiSMASH regions (sm_clusters.tsv.gz)\n"
    "sm_clusters_overlap_ge1gene\t3\trule A (sm_clusters.tsv.gz)\n"
    "sm_clusters_overlap_ge50pct\t1\trule B (sm_clusters.tsv.gz)\n"
    "sm_clusters_with_captain\t3\tDUF3435 hit E<=1e-5 in a cluster gene\n")
_AS_BY_TYPE = (
    "type\tn_clusters\tn_overlap_ge1gene\tn_overlap_ge50pct\tn_with_captain\n"
    "terpene\t12\t1\t0\t1\n"
    "T1PKS\t8\t1\t1\t1\n"
    "NRPS\t6\t0\t0\t0\n"
    "ALL\t44\t3\t1\t3\n")
_AS_PER_ELEMENT = (
    "element_id\tlabel\tconfidence\tcontig\tstart\tend\tn_genes_in_element\tn_sm_clusters_ge1gene\t"
    "n_sm_clusters_ge50pct\tsm_clusters\tn_sm_cluster_genes_in_element\tcaptain_genes_in_element\n"
    "AF293_e00029\tLamia h4\thigh\tCM000169.1\t6458\t394211\t71\t1\t1\t"
    "CM000169.1.region001(T1PKS;13/13 genes)\t13\tAFUA_1G00890\n"
    "AF293_e00039\tGnosis h1\thigh\tCM000174.1\t2214324\t2290180\t21\t0\t0\t\t0\tAFUA_6G09300\n"
    "AF293_e00041\tJanus (likely) h2 (likely)\tunknown\tCM000175.1\t1807654\t1860556\t6\t0\t0\t\t0\t\n")
_AS_PERM = (
    "statistic\tobserved\tn_perm\tseed\tnull_mean\tnull_q95\tnull_max\tp_perm_greater\tp_perm_less\n"
    "clusters_overlap_ge1gene\t3\t1000\t20261017\t1.884\t4\t6\t0.28272\t0.8951\n"
    "clusters_overlap_ge50pct\t1\t1000\t20261017\t1.206\t3\t5\t0.71528\t0.64036\n")
_AS_FAM = (
    "gene_set\tbin\tdirection\tin_island\tn_families\n"
    "core_gene_families\tcloud\tgain\tno\t3\n"
    "core_gene_families\tcloud\tgain\tyes\t3\n"
    "core_gene_families\tcore\tconserved\tno\t44\n")


def _antismash_dir(root: Path) -> Path:
    d = root / "analysis" / "starship_island_overlap" / "r1_antismash"
    (d / "figures").mkdir(parents=True)
    (d / "summary.tsv").write_text(_AS_SUMMARY)
    (d / "sm_by_type.tsv").write_text(_AS_BY_TYPE)
    (d / "sm_starship_overlap.tsv").write_text(_AS_PER_ELEMENT)
    (d / "sm_permutation.tsv").write_text(_AS_PERM)
    (d / "sm_family_class_direction.tsv").write_text(_AS_FAM)
    (d / "PROVENANCE.yaml").write_text("antismash:\n  version: antiSMASH 8.0.4\n")
    (d / "PROVENANCE_analysis.yaml").write_text("derived_by: antismash_sm_starship.py\n")
    (d / "figures" / "af293_antismash_sm_starship.png").write_bytes(b"png")
    return d


def test_stage_starship_overlap_separate_antismash_dir(tmp_path, capsys):
    import re as _re
    import sync_pangenome_report as sp
    d = _sso_dir(tmp_path, with_sm=False)
    a = _antismash_dir(tmp_path)
    run = tmp_path / "run"
    run.mkdir()
    assert sp.stage_starship_overlap(d, run, a)
    dest = run / "starship_overlap"
    # the overlap dir's own summary/provenance are not overwritten by the SM dir's
    assert (dest / "summary.tsv").read_text() == _SSO_SUMMARY
    assert (dest / "sm_summary.tsv").read_text() == _AS_SUMMARY
    assert (dest / "PROVENANCE_sm.yaml").read_text().startswith("antismash:")
    assert (dest / "PROVENANCE_sm_analysis.yaml").is_file()
    assert (dest / "af293_antismash_sm_starship.png").is_file()
    page = (run / "starship_overlap.html").read_text()
    assert "antiSMASH 8.0.4 regions" in page
    assert "Table S7 elements and SM clusters" in page and "CM000169.1.region001(T1PKS;13/13 genes)" in page
    assert "SM cluster position permutation" in page and "0.28272" in page
    assert "SM core-gene families by class, direction and island" in page
    assert "<td>sm_clusters_overlap_ge1gene</td><td class=\"n\">3</td>" in page
    assert "<td>sm_clusters</td><td class=\"n\">44</td>" in page
    assert "does not map SM cluster genes" not in page
    assert "mapped to pangenome families" in page and "9627 of 9630" in page
    assert "WARNING" not in capsys.readouterr().err
    for href in _re.findall(r'(?:href|src)="(starship_overlap/[^"]+)"', page):
        assert (run / href).is_file(), href


def test_stage_starship_overlap_metric_schema_still_works(tmp_path):
    import sync_pangenome_report as sp
    d = _sso_dir(tmp_path, with_sm=True)   # cocci-style metric/value/definition table
    run = tmp_path / "run"
    run.mkdir()
    assert sp.stage_starship_overlap(d, run)
    page = (run / "starship_overlap.html").read_text()
    assert "antiSMASH 8.0.1 regions" in page
    assert "Table S7 elements and SM clusters" not in page
    assert "does not map SM cluster genes" in page   # no sm_family_class_direction.tsv


def test_stage_starship_overlap_unexpected_schemas_do_not_crash(tmp_path, capsys):
    import sync_pangenome_report as sp
    d = _sso_dir(tmp_path, with_sm=False)
    (d / "enrichment.tsv").write_text("foo\tbar\n1\t2\n")
    a = _antismash_dir(tmp_path)
    (a / "sm_starship_overlap.tsv").write_text("foo\tbar\n1\t2\n")
    (a / "summary.tsv").write_text("name\tnumber\nx\t1\n")
    run = tmp_path / "run"
    run.mkdir()
    assert sp.stage_starship_overlap(d, run, a)
    page = (run / "starship_overlap.html").read_text()
    assert "Not available: sm_starship_overlap.tsv has unexpected columns foo, bar." in page
    assert "enrichment.tsv: not available" in page
    err = capsys.readouterr().err
    assert "sm_starship_overlap.tsv: unexpected columns" in err and "expected columns" in err
    (d / "summary.tsv").write_text("a\tb\n1\t2\n")          # overlap summary itself wrong: no page
    run2 = tmp_path / "run2"
    run2.mkdir()
    assert not sp.stage_starship_overlap(d, run2, a)
    assert not (run2 / "starship_overlap.html").exists()
