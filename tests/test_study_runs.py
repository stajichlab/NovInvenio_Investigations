import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

import generate_docs  # noqa: E402
import migrate_docs_to_runs  # noqa: E402
import study_pages  # noqa: E402
from study_runs import (  # noqa: E402
    REDIRECT_MARKER, REDIRECT_MAX_BYTES, load_publish_target, rebuild_study_pages,
    render_redirect, set_current,
)

REPO_ROOT = Path(__file__).parent.parent


def _run_json(study_docs: Path, run: str, published: str = "2026-09-01", **extra) -> None:
    (study_docs / run).mkdir(parents=True, exist_ok=True)
    (study_docs / run / "run.json").write_text(json.dumps({"run": run, "published": published, **extra}))


def _folder(root: Path, domain: str, folder: str, publish: str | None,
            params: str = "", staged: bool = True) -> Path:
    d = root / "studies" / domain / folder
    d.mkdir(parents=True)
    (d / "species.csv").write_text("Short,Group,TaxonGroup\na,IN,Pezizo\nb,IN,Pezizo\nc,OUT,Other\n")
    if publish is not None:
        (d / "publish.yaml").write_text(publish)
    if params:
        (d / "run_params.txt").write_text(params)
    if staged:
        (d / "config.csv").write_text("x\n")
        (d / "data_dir").mkdir()
    return d


# --- publish.yaml ----------------------------------------------------------


def test_publish_target_missing_is_none(tmp_path):
    assert load_publish_target(tmp_path) is None


def test_publish_target_study_and_run(tmp_path):
    (tmp_path / "publish.yaml").write_text("study: pezizo_set1\nrun: mmseqs-cov0.3\n")
    t = load_publish_target(tmp_path)
    assert (t.study, t.run) == ("pezizo_set1", "mmseqs-cov0.3")


def test_publish_target_study_only(tmp_path):
    (tmp_path / "publish.yaml").write_text("study: cocci\n")
    assert load_publish_target(tmp_path).run is None


@pytest.mark.parametrize("bad", ["study: a/b\n", "study: ..\n", "study: x\nrun: 'a b'\n", "run: x\n"])
def test_publish_target_rejects_bad_names(tmp_path, bad):
    (tmp_path / "publish.yaml").write_text(bad)
    with pytest.raises(ValueError):
        load_publish_target(tmp_path)


# --- study pages -----------------------------------------------------------


def test_single_run_becomes_current(tmp_path):
    _run_json(tmp_path, "r1")
    rebuild_study_pages(tmp_path, "fungi", "s")
    assert json.loads((tmp_path / "study.json").read_text()) == {"current": "r1"}
    assert "url=r1/report.html" in (tmp_path / "index.html").read_text()


def test_several_runs_without_current_marks_none(tmp_path):
    _run_json(tmp_path, "r1")
    _run_json(tmp_path, "r2")
    rebuild_study_pages(tmp_path, "fungi", "s")
    assert not (tmp_path / "study.json").exists()
    assert "url=report.html" in (tmp_path / "index.html").read_text()


def test_current_run_listed_first_and_marked(tmp_path):
    _run_json(tmp_path, "old", published="2026-01-01")
    _run_json(tmp_path, "new", published="2026-09-01")
    set_current(tmp_path, "old")
    rebuild_study_pages(tmp_path, "fungi", "s")
    html = (tmp_path / "report.html").read_text()
    assert html.index('href="old/report.html"') < html.index('href="new/report.html"')
    assert "Current result" in html and "archived" in html


def test_set_current_requires_existing_run(tmp_path):
    with pytest.raises(ValueError, match="run.json not found"):
        set_current(tmp_path, "nope")


def test_stale_current_is_replaced_when_one_run(tmp_path):
    _run_json(tmp_path, "r1")
    (tmp_path / "study.json").write_text('{"current": "gone"}')
    rebuild_study_pages(tmp_path, "fungi", "s")
    assert json.loads((tmp_path / "study.json").read_text()) == {"current": "r1"}


def test_redirect_stub_is_small_and_marked():
    stub = render_redirect("../x/novelties.html")
    assert REDIRECT_MARKER in stub and len(stub.encode()) <= REDIRECT_MAX_BYTES


# --- bin/study_pages.py ----------------------------------------------------


def test_target_exit_codes(tmp_path, capsys):
    _folder(tmp_path, "fungi", "none", None)
    _folder(tmp_path, "fungi", "studyonly", "study: s\n")
    _folder(tmp_path, "fungi", "full", "study: s\nrun: r\n")
    root = ["--repo-root", str(tmp_path)]
    assert study_pages.main(root + ["target", "fungi/none"]) == 3
    assert study_pages.main(root + ["target", "fungi/studyonly"]) == 4
    assert study_pages.main(root + ["target", "fungi/full"]) == 0
    assert capsys.readouterr().out.strip() == "s r"


def test_record_run_writes_metadata(tmp_path):
    _folder(tmp_path, "fungi", "pezizo_set1_cluster", "study: pezizo_set1\nrun: mmseqs\n",
            params="# comment\n--run_tool diamond\n\n--cluster_tool mmseqs  # inline\n")
    report = tmp_path / "docs" / "fungi" / "pezizo_set1" / "mmseqs" / "report.html"
    report.parent.mkdir(parents=True)
    report.write_text("<html></html>")
    study_pages.main(["--repo-root", str(tmp_path), "record-run", "fungi/pezizo_set1_cluster",
                      "--report", str(report), "--published", "2026-09-10", "--current"])
    meta = json.loads((report.parent / "run.json").read_text())
    assert meta["params"] == ["--run_tool diamond", "--cluster_tool mmseqs"]
    assert meta["n_species"] == 3
    assert meta["published"] == "2026-09-10"
    assert meta["source_dir"] == "results/pezizo_set1_cluster"
    assert meta["pipeline"] == "main.nf"
    assert json.loads((report.parent.parent / "study.json").read_text()) == {"current": "mmseqs"}


# --- bin/generate_docs.py --------------------------------------------------


@pytest.fixture
def gallery(tmp_path, monkeypatch):
    (tmp_path / "conf").mkdir()
    (tmp_path / "conf" / "domains.yaml").write_text("- {slug: fungi, name: Fungi, desc: d}\n")
    monkeypatch.setattr(generate_docs, "STUDIES_ROOT", tmp_path / "studies")
    monkeypatch.setattr(generate_docs, "DOCS_ROOT", tmp_path / "docs")
    monkeypatch.setattr(generate_docs, "DOMAINS_YAML", tmp_path / "conf" / "domains.yaml")
    monkeypatch.setattr(generate_docs, "PIXI_TOML", tmp_path / "pixi.toml")
    return tmp_path


def test_gallery_cards(gallery):
    _folder(gallery, "fungi", "pezizo_set1", "study: pezizo_set1\nrun: pw\n")
    _folder(gallery, "fungi", "pezizo_set1_cluster", "study: pezizo_set1\nrun: mm\n")
    _folder(gallery, "fungi", "newstudy", "study: newstudy\nrun: pw\n")
    _folder(gallery, "fungi", "pending_one", "study: pending_one\nrun: pw\n", staged=False)
    _folder(gallery, "fungi", "dmnd_test", None)
    study_docs = gallery / "docs" / "fungi" / "pezizo_set1"
    _run_json(study_docs, "pw")
    _run_json(study_docs, "mm")
    set_current(study_docs, "pw")

    assert generate_docs.main() == 0
    html = (gallery / "docs" / "fungi" / "index.html").read_text()
    assert 'href="pezizo_set1/pw/report.html"' in html
    assert 'href="pezizo_set1/report.html">2 runs</a>' in html
    assert "2 ingroup" in html
    assert "newstudy" in html and "data staged" in html
    assert "pending one" in html and "not yet started" in html
    assert "dmnd" not in html
    assert html.count("pezizo set1") == 1


# --- bin/migrate_docs_to_runs.py -------------------------------------------


def _git(root, *args):
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "config", "user.email", "t@t")
    _git(tmp_path, "config", "user.name", "t")
    (tmp_path / ".gitignore").write_text("docs/*/*/novelties.html\ndocs/*/*/*/novelties.html\n")
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "study_pages.py").symlink_to(REPO_ROOT / "bin" / "study_pages.py")
    _folder(tmp_path, "fungi", "pz", "study: pz\nrun: pw\n", params="--run_tool diamond\n")
    _folder(tmp_path, "fungi", "pz_cluster", "study: pz\nrun: mm\n")
    for f in ("pz", "pz_cluster"):
        d = tmp_path / "docs" / "fungi" / f
        (d / "alignments").mkdir(parents=True)
        (d / "report.html").write_text(f"<html>{f} report</html>")
        (d / "alignment.html").write_text("<html>viewer</html>")
        (d / "index.html").write_text("<html>old redirect</html>")
        (d / "novelties.html").write_text("real novelties " * 500)
        (d / "alignments" / "shard.json.gz").write_bytes(b"x")
    _git(tmp_path, "add", ".gitignore", "studies", "docs/fungi/pz/report.html",
         "docs/fungi/pz/alignment.html", "docs/fungi/pz/index.html",
         "docs/fungi/pz_cluster/report.html", "docs/fungi/pz_cluster/alignment.html",
         "docs/fungi/pz_cluster/index.html")
    _git(tmp_path, "commit", "-q", "-m", "init", "--date", "2026-09-10T00:00:00")
    return tmp_path


def test_migrate_dry_run_changes_nothing(repo):
    before = sorted(p.relative_to(repo) for p in (repo / "docs").rglob("*"))
    assert migrate_docs_to_runs.main(["--repo-root", str(repo)]) == 0
    assert sorted(p.relative_to(repo) for p in (repo / "docs").rglob("*")) == before


def test_migrate_apply(repo):
    assert migrate_docs_to_runs.main(["--repo-root", str(repo), "--apply"]) == 0
    study = repo / "docs" / "fungi" / "pz"
    pw, mm = study / "pw", study / "mm"
    for run in (pw, mm):
        assert (run / "report.html").exists() and (run / "alignment.html").exists()
        assert (run / "novelties.html").exists() and (run / "alignments" / "shard.json.gz").exists()
    meta = json.loads((pw / "run.json").read_text())
    assert meta["params"] == ["--run_tool diamond"]
    assert meta["published"] == subprocess.run(
        ["git", "-C", str(repo), "log", "-1", "--format=%cs"], capture_output=True, text=True).stdout.strip()
    assert "migrated from docs/fungi/pz/" in meta["note"]

    # same-dir case: report.html is now the run list; other old names are stubs
    assert "pw/report.html" in (study / "report.html").read_text()
    assert "url=pw/novelties.html" in (study / "novelties.html").read_text()
    assert "url=pw/alignment.html" in (study / "alignment.html").read_text()
    # other-dir case: every old name is a stub into ../pz/mm/
    old = repo / "docs" / "fungi" / "pz_cluster"
    for name in ("report.html", "novelties.html", "alignment.html"):
        assert f"url=../pz/mm/{name}" in (old / name).read_text()
    assert "url=../pz/mm/report.html" in (old / "index.html").read_text()

    # stubs are staged (force-added past the ignore rule); real reports are not
    staged = subprocess.run(["git", "-C", str(repo), "diff", "--cached", "--name-only"],
                            capture_output=True, text=True).stdout.split()
    assert "docs/fungi/pz_cluster/novelties.html" in staged
    assert "docs/fungi/pz/mm/novelties.html" not in staged
    assert "docs/fungi/pz/pw/run.json" in staged


def test_migrate_is_rerunnable_and_moves_leftovers(repo):
    migrate_docs_to_runs.main(["--repo-root", str(repo), "--apply"])
    # Simulate a second checkout that got the migration by merge: an untracked
    # real report is still at the old place.
    (repo / "docs" / "fungi" / "pz_cluster" / "core.html").write_text("real core " * 500)
    assert migrate_docs_to_runs.main(["--repo-root", str(repo), "--apply"]) == 0
    assert (repo / "docs" / "fungi" / "pz" / "mm" / "core.html").read_text().startswith("real core")
    core_stub = repo / "docs" / "fungi" / "pz_cluster" / "core.html"
    assert "url=../pz/mm/core.html" in core_stub.read_text()
    staged = subprocess.run(["git", "-C", str(repo), "diff", "--cached", "--name-only"],
                            capture_output=True, text=True).stdout.split()
    assert "docs/fungi/pz_cluster/core.html" in staged
    # stubs untouched
    assert REDIRECT_MARKER in (repo / "docs" / "fungi" / "pz_cluster" / "novelties.html").read_text()
    assert migrate_docs_to_runs.main(["--repo-root", str(repo), "--apply"]) == 0


# --- guard on the real repo ------------------------------------------------


def test_tracked_set_level_reports_are_only_redirect_stubs():
    """Set-level novelties/core/losses.html are gitignored; the only ones that
    may be committed are old-URL redirect stubs (force-added by the migration).
    A real report committed there is exactly the bloat CLAUDE.md forbids."""
    out = subprocess.run(["git", "-C", str(REPO_ROOT), "ls-files", "docs"],
                         capture_output=True, text=True, check=True).stdout.split()
    for rel in out:
        parts = rel.split("/")
        if len(parts) == 4 and parts[3] in ("novelties.html", "core.html", "losses.html"):
            p = REPO_ROOT / rel
            assert p.stat().st_size <= REDIRECT_MAX_BYTES, rel
            assert REDIRECT_MARKER in p.read_text(), rel
