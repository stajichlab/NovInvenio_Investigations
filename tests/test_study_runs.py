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
    REDIRECT_MARKER, REDIRECT_MAX_BYTES, inject_run_report_notes, load_publish_target,
    rebuild_study_pages, render_redirect, set_current,
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


def test_publish_target_notes(tmp_path):
    (tmp_path / "publish.yaml").write_text("study: s\nrun: r\nnotes:\n  - first\n  - ' second '\n")
    assert load_publish_target(tmp_path).notes == ("first", "second")
    (tmp_path / "publish.yaml").write_text("study: s\n")
    assert load_publish_target(tmp_path).notes == ()


@pytest.mark.parametrize("bad", ["notes: text", "notes: [1]", "notes: ['']"])
def test_publish_target_rejects_bad_notes(tmp_path, bad):
    (tmp_path / "publish.yaml").write_text(f"study: s\n{bad}\n")
    with pytest.raises(ValueError):
        load_publish_target(tmp_path)


_REPORT = ('<body>\n<header class="top">\n  <h1>x</h1>\n  </header>\n\n'
           '  <section class="card">\n    <div class="cards">\n'
           '      <a class="card-link" href="novelties.html">n</a>\n    </div>\n  </section>\n</body>\n')


def _report_dir(tmp_path):
    import hashlib
    d = tmp_path / "run"
    d.mkdir()
    (d / "report.html").write_text(_REPORT)
    (d / "run.json").write_text(json.dumps(
        {"run": "run", "report_sha256": hashlib.sha256(_REPORT.encode()).hexdigest()}))
    return d


def test_inject_notes_adds_replaces_and_removes(tmp_path):
    import hashlib
    d = _report_dir(tmp_path)
    inject_run_report_notes(d, ["a <b> note"])
    html = (d / "report.html").read_text()
    assert html.index("</header>") < html.index("a &lt;b&gt; note") < html.index('class="cards"')
    assert html.count("nii-run-notes:start") == 1
    sha = json.loads((d / "run.json").read_text())["report_sha256"]
    assert sha == hashlib.sha256((d / "report.html").read_bytes()).hexdigest()
    assert "already current" in inject_run_report_notes(d, ["a <b> note"])
    inject_run_report_notes(d, ["other"])
    html = (d / "report.html").read_text()
    assert "other" in html and "a &lt;b&gt;" not in html and html.count("nii-run-notes:start") == 1
    inject_run_report_notes(d, [])
    assert (d / "report.html").read_text() == _REPORT
    assert json.loads((d / "run.json").read_text())["report_sha256"] == hashlib.sha256(_REPORT.encode()).hexdigest()


def test_inject_notes_leaves_mismatched_sha_alone(tmp_path):
    d = _report_dir(tmp_path)
    (d / "run.json").write_text(json.dumps({"run": "run", "report_sha256": "0" * 64}))
    assert "WARNING" in inject_run_report_notes(d, ["n"])
    assert json.loads((d / "run.json").read_text())["report_sha256"] == "0" * 64


def test_inject_notes_without_header_warns(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    (d / "report.html").write_text("<body></body>")
    assert "WARNING" in inject_run_report_notes(d, ["n"])
    assert (d / "report.html").read_text() == "<body></body>"


def test_inject_notes_command_reads_publish_yaml(tmp_path):
    _folder(tmp_path, "fungi", "s1", "study: s\nrun: r\nnotes: [flagged]\n")
    run = tmp_path / "docs" / "fungi" / "s" / "r"
    run.mkdir(parents=True)
    (run / "report.html").write_text(_REPORT)
    assert study_pages.main(["--repo-root", str(tmp_path), "inject-notes", "fungi/s1"]) == 0
    assert "flagged" in (run / "report.html").read_text()


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


def test_nested_folder_target_and_record_run(tmp_path, capsys):
    """A folder one level deeper (studies/fungi/<group>/<contrast>/) publishes under
    the flat docs/<domain>/<study>/<run>/ layout; its results dir is named after
    the folder's last component, as bin/run_study.sh --project uses basename."""
    _folder(tmp_path, "fungi", "yeast_convergence/yc_x", "study: yc_x\nrun: pairwise-diamond\n")
    root = ["--repo-root", str(tmp_path)]
    assert study_pages.main(root + ["target", "fungi/yeast_convergence/yc_x"]) == 0
    assert capsys.readouterr().out.strip() == "yc_x pairwise-diamond"
    report = tmp_path / "docs" / "fungi" / "yc_x" / "pairwise-diamond" / "report.html"
    report.parent.mkdir(parents=True)
    report.write_text("<html></html>")
    study_pages.main(root + ["record-run", "fungi/yeast_convergence/yc_x", "--report", str(report),
                             "--published", "2026-10-22"])
    meta = json.loads((report.parent / "run.json").read_text())
    assert meta["folder"] == "studies/fungi/yeast_convergence/yc_x"
    assert meta["source_dir"] == "results/yc_x"


@pytest.mark.parametrize("bad", ["fungi", "fungi/", "fungi/../x", "fungi/./x", "/x"])
def test_folder_argument_rejects_bad_paths(tmp_path, bad):
    with pytest.raises(SystemExit):
        study_pages.main(["--repo-root", str(tmp_path), "target", bad])


def test_sync_reports_domain_is_first_component():
    """sync_reports.sh must not use dirname for DOMAIN (wrong for nested folders)."""
    text = (REPO_ROOT / "bin" / "sync_reports.sh").read_text()
    assert 'DOMAIN="${STUDY%%/*}"' in text
    assert 'DOMAIN="$(dirname "$STUDY")"' not in text


def test_publish_target_related(tmp_path):
    (tmp_path / "publish.yaml").write_text(
        "study: s\nrun: r\nrelated:\n  - {title: Comparison, href: yeast_convergence/index.html}\n")
    assert load_publish_target(tmp_path).related == (("Comparison", "yeast_convergence/index.html"),)
    (tmp_path / "publish.yaml").write_text("study: s\n")
    assert load_publish_target(tmp_path).related == ()


@pytest.mark.parametrize("entry", [
    "{title: X, href: 'https://example.org/'}", "{title: X, href: /abs/path.html}",
    "{title: X, href: ../other/index.html}", "{title: '', href: a/index.html}", "{href: a/index.html}",
    "just a string",
])
def test_publish_target_related_rejects_bad_entries(tmp_path, entry):
    (tmp_path / "publish.yaml").write_text(f"study: s\nrelated:\n  - {entry}\n")
    with pytest.raises(ValueError, match="related"):
        load_publish_target(tmp_path)


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
    assert "enrichment" not in html

    # An enrichment page in the current run gets a card link on the next run.
    (study_docs / "pw" / "enrichment.html").write_text("x")
    assert generate_docs.main() == 0
    html = (gallery / "docs" / "fungi" / "index.html").read_text()
    assert 'href="pezizo_set1/pw/enrichment.html">enrichment</a>' in html
    assert "enrichment_loss" not in html


def test_gallery_lists_nested_folders(gallery):
    """Unpublished contrasts inside a grouping folder (no species.csv of its own)
    get cards too; a grouping folder's own README is not mistaken for a study."""
    _folder(gallery, "fungi", "yeast_convergence/yc_alpha", "study: yc_alpha\nrun: pw\n", staged=False)
    _folder(gallery, "fungi", "yeast_convergence/yc_beta", "study: yc_beta\nrun: pw\n")
    _folder(gallery, "fungi", "yeast_convergence/yc_unpub", None)
    (gallery / "studies" / "fungi" / "yeast_convergence" / "batches").mkdir()
    assert generate_docs.main() == 0
    html = (gallery / "docs" / "fungi" / "index.html").read_text()
    assert "yc alpha" in html and "yc beta" in html
    assert "yc unpub" not in html


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


# --- enrichment link text from the page's vocabulary tag (2026-10-16) ----------
from study_runs import (  # noqa: E402
    RUN_EXTRA_PAGES, inject_run_report_links, page_vocabularies, run_extra_links, vocab_meta_tag,
)

_LINK_REPORT = ('<html><body>\n      <a class="card-link" href="losses.html"><span class="card-link-title">L'
           "</span></a>\n</body></html>\n")


def _page_with(vocabs):
    head = "" if vocabs is None else vocab_meta_tag(vocabs) + "\n"
    return f"<!doctype html>\n<html><head>\n{head}<title>x</title></head><body></body></html>\n"


def test_vocab_meta_tag_order_and_parse(tmp_path):
    assert vocab_meta_tag({"go", "pfam", "interpro"}) == \
        '<meta name="nii-enrichment-vocabularies" content="pfam,interpro,go">'
    (tmp_path / "p.html").write_text(_page_with(["pfam"]))
    assert page_vocabularies(tmp_path / "p.html") == ["pfam"]
    (tmp_path / "q.html").write_text(_page_with(None))
    assert page_vocabularies(tmp_path / "q.html") is None


@pytest.mark.parametrize("vocabs", [None, ["pfam", "interpro", "go"]])
def test_full_or_untagged_pages_keep_default_link_text(tmp_path, vocabs):
    (tmp_path / "enrichment.html").write_text(_page_with(vocabs))
    (tmp_path / "enrichment_loss.html").write_text(_page_with(vocabs))
    (tmp_path / "report.html").write_text(_LINK_REPORT)
    assert run_extra_links(tmp_path) == [("enrichment.html", "enrichment"),
                                         ("enrichment_loss.html", "loss enrichment")]
    inject_run_report_links(tmp_path)
    report = (tmp_path / "report.html").read_text()
    for fname, _, title, desc in RUN_EXTRA_PAGES:
        assert (f'<a class="card-link" href="{fname}"><span class="card-link-title">{title}</span>'
                f'<span class="card-link-desc">{desc}</span></a>') in report


def test_pfam_only_pages_name_the_vocabulary(tmp_path):
    (tmp_path / "enrichment.html").write_text(_page_with(["pfam"]))
    (tmp_path / "enrichment_loss.html").write_text(_page_with(["pfam"]))
    (tmp_path / "report.html").write_text(_LINK_REPORT)
    assert run_extra_links(tmp_path) == [("enrichment.html", "enrichment (Pfam only)"),
                                         ("enrichment_loss.html", "loss enrichment (Pfam only)")]
    inject_run_report_links(tmp_path)
    report = (tmp_path / "report.html").read_text()
    assert ("Pfam term enrichment of the novelty candidates against the proteome background. "
            "InterPro and GO: not available.") in report
    assert "Pfam, InterPro and GO" not in report
    # Run list (rebuild_study_pages) uses the same labels.
    (tmp_path / "run.json").write_text(json.dumps({"run": tmp_path.name, "published": "2026-10-16"}))
    pages = rebuild_study_pages(tmp_path.parent, "fungi", "s")
    assert pages and "enrichment (Pfam only)" in (tmp_path.parent / "report.html").read_text()


def test_gallery_domain_pages(gallery):
    """conf/domains.yaml `pages:` adds a card for a hand-made page that is not a study."""
    (gallery / "conf" / "domains.yaml").write_text(
        "- slug: fungi\n  name: Fungi\n  desc: d\n  pages:\n"
        "    - {title: yeast_convergence comparison, href: yeast_convergence/index.html, desc: 'Across 16 contrasts'}\n")
    _folder(gallery, "fungi", "newstudy", "study: newstudy\nrun: pw\n")
    assert generate_docs.main() == 0
    html = (gallery / "docs" / "fungi" / "index.html").read_text()
    assert 'href="yeast_convergence/index.html"' in html
    assert "yeast_convergence comparison" in html and "Across 16 contrasts" in html


def test_gallery_domain_with_only_pages_is_written(gallery):
    (gallery / "conf" / "domains.yaml").write_text(
        "- slug: fungi\n  name: Fungi\n  desc: d\n  pages:\n    - {title: T, href: x/index.html, desc: D}\n")
    assert generate_docs.main() == 0
    assert 'href="x/index.html"' in (gallery / "docs" / "fungi" / "index.html").read_text()


# --- gallery groups ---------------------------------------------------------


def _groups_yaml(gallery, extra=""):
    (gallery / "conf" / "domains.yaml").write_text(
        "- slug: fungi\n  name: Fungi\n  desc: d\n  groups:\n"
        "    - slug: yc\n      title: yeast convergence\n      member_label: contrasts\n"
        "      desc: grouped\n      page: {title: comparison, href: yc/comparison.html}\n" + extra)


def _nested_pair(gallery):
    _folder(gallery, "fungi", "yc/yc_alpha", "study: yc_alpha\nrun: pw\n")
    _folder(gallery, "fungi", "yc/yc_beta", "study: yc_beta\nrun: pw\n")
    _folder(gallery, "fungi", "solo", "study: solo\nrun: pw\n")


def test_publish_target_group_key(tmp_path):
    (tmp_path / "publish.yaml").write_text("study: s\n")
    t = load_publish_target(tmp_path)
    assert (t.group, t.group_set) == (None, False)
    (tmp_path / "publish.yaml").write_text("study: s\ngroup: g\n")
    assert load_publish_target(tmp_path).group == "g"
    (tmp_path / "publish.yaml").write_text("study: s\ngroup: null\n")
    t = load_publish_target(tmp_path)
    assert (t.group, t.group_set) == (None, True)
    (tmp_path / "publish.yaml").write_text("study: s\ngroup: ../x\n")
    with pytest.raises(ValueError):
        load_publish_target(tmp_path)


def test_gallery_group_card_and_page(gallery):
    _groups_yaml(gallery)
    _nested_pair(gallery)
    (gallery / "docs" / "fungi" / "yc").mkdir(parents=True)
    (gallery / "docs" / "fungi" / "yc" / "comparison.html").write_text("x")
    assert generate_docs.main() == 0
    domain = (gallery / "docs" / "fungi" / "index.html").read_text()
    assert 'href="yc/index.html"' in domain and "2 contrasts" in domain
    assert "yc alpha" not in domain and "yc beta" not in domain and "solo" in domain
    group = (gallery / "docs" / "fungi" / "yc" / "index.html").read_text()
    assert "ni-group-page" in group and "2 contrasts" in group
    assert "yc alpha" in group and "yc beta" in group and "solo" not in group
    assert 'href="../yc/comparison.html"' in group
    top = (gallery / "docs" / "index.html").read_text()
    assert "3 studies" in top


def test_gallery_group_member_links_are_relative_to_group_page(gallery):
    _groups_yaml(gallery)
    _nested_pair(gallery)
    for name in ("yc_alpha", "yc_beta"):
        d = gallery / "docs" / "fungi" / name
        _run_json(d, "pw")
        set_current(d, "pw")
        (d / "pw" / "enrichment.html").write_text("x")
    assert generate_docs.main() == 0
    group = (gallery / "docs" / "fungi" / "yc" / "index.html").read_text()
    assert 'href="../yc_alpha/pw/report.html"' in group
    assert 'href="../yc_alpha/report.html">1 run</a>' in group
    assert 'href="../yc_alpha/pw/enrichment.html"' in group


def test_gallery_unlisted_group_stays_ungrouped(gallery, capsys):
    _nested_pair(gallery)  # domains.yaml has no `groups:`
    assert generate_docs.main() == 0
    assert "yc alpha" in (gallery / "docs" / "fungi" / "index.html").read_text()
    assert not (gallery / "docs" / "fungi" / "yc").exists()
    assert "not listed under `groups:`" in capsys.readouterr().err


def test_gallery_group_opt_out_and_override(gallery):
    _groups_yaml(gallery)
    _folder(gallery, "fungi", "yc/yc_alpha", "study: yc_alpha\nrun: pw\ngroup: null\n")
    _folder(gallery, "fungi", "yc/yc_beta", "study: yc_beta\nrun: pw\n")
    _folder(gallery, "fungi", "other", "study: other\nrun: pw\ngroup: yc\n")
    assert generate_docs.main() == 0
    group = (gallery / "docs" / "fungi" / "yc" / "index.html").read_text()
    assert "yc beta" in group and "other" in group and "yc alpha" not in group


def test_gallery_one_member_group_shows_member_and_page(gallery, capsys):
    _groups_yaml(gallery)
    _folder(gallery, "fungi", "yc/yc_alpha", "study: yc_alpha\nrun: pw\n")
    assert generate_docs.main() == 0
    html = (gallery / "docs" / "fungi" / "index.html").read_text()
    assert "yc alpha" in html and 'href="yc/comparison.html"' in html
    assert not (gallery / "docs" / "fungi" / "yc").exists()
    assert "1 member" in capsys.readouterr().err


def test_gallery_group_conflict_between_folders_fails(gallery):
    _groups_yaml(gallery)
    _folder(gallery, "fungi", "yc/yc_alpha", "study: yc_alpha\nrun: pw\n")
    _folder(gallery, "fungi", "yc_alpha_flat", "study: yc_alpha\nrun: mm\n")
    with pytest.raises(SystemExit, match="different groups"):
        generate_docs.main()


def test_gallery_group_slug_equal_to_study_name_fails(gallery):
    _groups_yaml(gallery)
    _folder(gallery, "fungi", "yc", "study: yc\nrun: pw\n")
    with pytest.raises(SystemExit, match="same name as a study"):
        generate_docs.main()


def test_gallery_group_refuses_to_overwrite_foreign_index(gallery):
    _groups_yaml(gallery)
    _nested_pair(gallery)
    d = gallery / "docs" / "fungi" / "yc"
    d.mkdir(parents=True)
    (d / "index.html").write_text("<html>hand made</html>")
    with pytest.raises(SystemExit, match="not a generated group page"):
        generate_docs.main()
    # Once the page has moved to its new address, the old file may be replaced.
    (d / "comparison.html").write_text("moved")
    assert generate_docs.main() == 0
    assert "ni-group-page" in (d / "index.html").read_text()


def test_gallery_warns_about_stale_group_page(gallery, capsys):
    _groups_yaml(gallery)
    _nested_pair(gallery)
    (gallery / "docs" / "fungi" / "yc").mkdir(parents=True)
    (gallery / "docs" / "fungi" / "yc" / "comparison.html").write_text("x")
    assert generate_docs.main() == 0
    (gallery / "conf" / "domains.yaml").write_text("- {slug: fungi, name: Fungi, desc: d}\n")
    capsys.readouterr()
    assert generate_docs.main() == 0
    assert "no longer exists" in capsys.readouterr().err
