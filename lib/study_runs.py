"""Study/run layout for the docs/ gallery (spec:
notes/superpowers/specs/2026-09-24-study-run-site-layout-design.md).

  docs/<domain>/<study>/report.html     run list, current run highlighted
  docs/<domain>/<study>/index.html      redirect -> current run's report.html
  docs/<domain>/<study>/study.json      {"current": "<run>"}
  docs/<domain>/<study>/<run>/run.json  run metadata (written by a sync script)

A studies/<domain>/<folder>/ is published only when it has a publish.yaml
naming its target study (and, for main.nf folders, its run). Publishing is
opt-in: a folder without one is skipped, never guessed.
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from html import escape
from pathlib import Path

import yaml

from site_pages import _page

SLUG_RE = re.compile(r"^[A-Za-z0-9._-]+$")

# A redirect stub is tiny and carries this marker. Tests use both to check
# that no real report is ever committed at set level.
REDIRECT_MARKER = 'http-equiv="refresh"'
REDIRECT_MAX_BYTES = 2048


@dataclass(frozen=True)
class PublishTarget:
    study: str
    run: str | None
    notes: tuple[str, ...] = ()
    related: tuple[tuple[str, str], ...] = ()   # (title, href relative to docs/<domain>/)
    group: str | None = None     # gallery group slug from `group:`; None = no value given
    group_set: bool = False      # True when the file has a `group:` key (null opts out)


def validate_related(entries, where: str) -> tuple[tuple[str, str], ...]:
    """publish.yaml `related:` -> ((title, href), ...). href is a path under docs/<domain>/:
    no scheme, no leading slash, no `..` segment."""
    if entries is None:
        return ()
    if not isinstance(entries, list):
        raise ValueError(f"{where}: related must be a list of {{title, href}}")
    out = []
    for e in entries:
        title = e.get("title") if isinstance(e, dict) else None
        href = e.get("href") if isinstance(e, dict) else None
        if not (isinstance(title, str) and title.strip() and isinstance(href, str) and href.strip()):
            raise ValueError(f"{where}: each related entry needs a non-empty title and href, got {e!r}")
        href = href.strip()
        if "://" in href or href.startswith("/") or ".." in href.split("/"):
            raise ValueError(f"{where}: related href must be a path under docs/<domain>/, got {href!r}")
        out.append((title.strip(), href))
    return tuple(out)


def validate_slug(name: str, what: str) -> str:
    if not isinstance(name, str) or not SLUG_RE.match(name) or name in (".", ".."):
        raise ValueError(f"{what} {name!r} is not a valid name (allowed: letters, digits, . _ -)")
    return name


def load_publish_target(folder_dir: Path) -> PublishTarget | None:
    """Read studies/<domain>/<folder>/publish.yaml. None when the file is
    absent (the folder is not published)."""
    path = folder_dir / "publish.yaml"
    if not path.exists():
        return None
    with open(path) as fh:
        data = yaml.safe_load(fh) or {}
    study = validate_slug(data.get("study"), f"{path}: study")
    run = data.get("run")
    if run is not None:
        run = validate_slug(run, f"{path}: run")
    notes = data.get("notes") or []
    if not isinstance(notes, list) or not all(isinstance(n, str) and n.strip() for n in notes):
        raise ValueError(f"{path}: notes must be a list of non-empty strings")
    related = validate_related(data.get("related"), f"{path}: related")
    group_set = "group" in data
    group = data.get("group")
    if group is not None:
        group = validate_slug(group, f"{path}: group")
    return PublishTarget(study=study, run=run, notes=tuple(n.strip() for n in notes), related=related,
                         group=group, group_set=group_set)


def render_redirect(target: str) -> str:
    """Client-side redirect stub (GitHub Pages has no server-side redirects).
    `target` is a relative URL."""
    t = escape(target, quote=True)
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        f'<meta http-equiv="refresh" content="0; url={t}">\n'
        "<title>Redirecting…</title>\n</head>\n<body>\n"
        f'<p>Moved to <a href="{t}">{t}</a>.</p>\n'
        f"<script>location.replace({json.dumps(target)});</script>\n"
        "</body>\n</html>\n"
    )


def load_runs(study_docs: Path) -> list[dict]:
    runs = []
    for meta_path in sorted(study_docs.glob("*/run.json")):
        with open(meta_path) as fh:
            runs.append(json.load(fh))
    return runs


def read_current(study_docs: Path) -> str | None:
    path = study_docs / "study.json"
    if not path.exists():
        return None
    with open(path) as fh:
        return json.load(fh).get("current")


def set_current(study_docs: Path, run: str) -> None:
    validate_slug(run, "run")
    if not (study_docs / run / "run.json").exists():
        raise ValueError(f"{study_docs / run}/run.json not found -- stage the run first")
    (study_docs / "study.json").write_text(json.dumps({"current": run}, indent=2) + "\n")


# Optional pages in a run folder. The run list, the gallery card and the run
# report link to a page only when its file exists, so a rerun of
# bin/generate_docs.py keeps the links and adds none for other runs.
# bin/generate_enrichment_page.py writes these files.
RUN_EXTRA_PAGES = (
    ("enrichment.html", "enrichment", "Functional enrichment (novelty)",
     "Pfam, InterPro and GO term enrichment of the novelty candidates against the proteome background."),
    ("enrichment_loss.html", "loss enrichment", "Functional enrichment (loss)",
     "Pfam, InterPro and GO term enrichment of the loss candidates against the proteome background."),
)


# bin/generate_enrichment_page.py writes the vocabularies that have data into
# each page as <meta name="nii-enrichment-vocabularies" content="pfam,interpro,go">.
# The link text names them when one or more is missing. A page without the
# tag (written before 2026-10-16) keeps the RUN_EXTRA_PAGES text.
VOCAB_META_NAME = "nii-enrichment-vocabularies"
VOCAB_LABELS = (("pfam", "Pfam"), ("interpro", "InterPro"), ("go", "GO"))
_VOCAB_META_RE = re.compile(
    r'<meta name="' + re.escape(VOCAB_META_NAME) + r'" content="([^"]*)"\s*/?>')


def vocab_meta_tag(vocabs: list[str] | tuple[str, ...] | set[str]) -> str:
    """The <meta> tag for the available vocabularies, in VOCAB_LABELS order."""
    keys = [k for k, _ in VOCAB_LABELS if k in set(vocabs)]
    return f'<meta name="{VOCAB_META_NAME}" content="{",".join(keys)}">'


def page_vocabularies(page: Path) -> list[str] | None:
    """Vocabularies from the page's meta tag, or None when the page has no tag."""
    try:
        head = page.read_text(encoding="utf-8", errors="replace")[:4096]
    except OSError:
        return None
    m = _VOCAB_META_RE.search(head)
    if m is None:
        return None
    return [v.strip() for v in m.group(1).split(",") if v.strip()]


def _join_labels(labels: list[str]) -> str:
    return labels[0] if len(labels) == 1 else ", ".join(labels[:-1]) + " and " + labels[-1]


def extra_page_text(run_dir: Path, fname: str, label: str, title: str, desc: str) -> tuple[str, str, str]:
    """(short label, title, description) for one RUN_EXTRA_PAGES entry. Unchanged
    when the page has no vocabulary tag or has all of Pfam, InterPro and GO."""
    vocabs = page_vocabularies(run_dir / fname)
    if vocabs is None:
        return label, title, desc
    have = [lbl for key, lbl in VOCAB_LABELS if key in vocabs]
    missing = [lbl for key, lbl in VOCAB_LABELS if key not in vocabs]
    if not missing:
        return label, title, desc
    what = desc.split(" term enrichment ", 1)[1] if " term enrichment " in desc else ""
    if not have:
        return (f"{label} (no annotation)", title,
                "No Pfam, InterPro or GO annotation; no enrichment test was run.")
    names = _join_labels(have)
    new_desc = (f"{names} term enrichment {what}" if what else f"{names} term enrichment.")
    new_desc += f" {_join_labels(missing)}: not available."
    return f"{label} ({names} only)", title, new_desc


def run_extra_links(run_dir: Path) -> list[tuple[str, str]]:
    """[(filename, short label), ...] for the RUN_EXTRA_PAGES that exist in run_dir.
    The label names the vocabularies when the page has fewer than all three
    (extra_page_text)."""
    return [(fname, extra_page_text(run_dir, fname, label, title, desc)[0])
            for fname, label, title, desc in RUN_EXTRA_PAGES if (run_dir / fname).is_file()]


_LINKS_START = "<!-- nii-run-extra-links:start -->"
_LINKS_END = "<!-- nii-run-extra-links:end -->"
_LINKS_BLOCK_RE = re.compile(
    r"[ \t]*" + re.escape(_LINKS_START) + r".*?" + re.escape(_LINKS_END) + r"\n?", re.S)
# The pipeline's run report (nf_NovInvenio COLLATE_REPORTS) lists its pages
# as <a class="card-link" ...> elements. The extra links go after the last
# of these, so they use the same style.
_CARD_LINK_RE = re.compile(r'<a class="card-link" href="(?:novelties|core|losses)\.html">.*?</a>\n?', re.S)
_INDENT = "      "


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def inject_run_report_links(run_dir: Path) -> str:
    """Write the RUN_EXTRA_PAGES links that exist in run_dir into
    run_dir/report.html, between marker comments. Idempotent: an earlier
    block is replaced, and no block is written when no extra page exists.

    report.html is a pipeline page; a sync of a new run overwrites it, so
    bin/sync_reports.sh calls this again (through
    bin/generate_enrichment_page.py --links-only) before it records the run.

    When run.json's report_sha256 matched report.html before the edit, it is
    updated to the new content. When it did not match, run.json is not
    changed. Returns a one-line status."""
    report = run_dir / "report.html"
    if not report.is_file():
        return f"no {report}; nothing to link"
    old_bytes = report.read_bytes()
    old = old_bytes.decode("utf-8")
    text = _LINKS_BLOCK_RE.sub("", old)
    pages = [(p[0], *extra_page_text(run_dir, *p)) for p in RUN_EXTRA_PAGES if (run_dir / p[0]).is_file()]
    if pages:
        anchors = list(_CARD_LINK_RE.finditer(text))
        if not anchors:
            return f"WARNING: {report} has no card-link anchor; links not added"
        block = _INDENT + _LINKS_START + "\n" + "".join(
            f'{_INDENT}<a class="card-link" href="{escape(fname)}"><span class="card-link-title">'
            f'{escape(title)}</span><span class="card-link-desc">{escape(desc)}</span></a>\n'
            for fname, _, title, desc in pages
        ) + _INDENT + _LINKS_END + "\n"
        end = anchors[-1].end()
        if not text[:end].endswith("\n"):
            block = "\n" + block
        text = text[:end] + block + text[end:]
    if text == old:
        return f"{report}: links already current"
    return f"{report}: links set to {[p[0] for p in pages]}" + _commit_report(run_dir, report, old_bytes, text)


_NOTES_START = "<!-- nii-run-notes:start -->"
_NOTES_END = "<!-- nii-run-notes:end -->"
_NOTES_BLOCK_RE = re.compile(
    r"[ \t]*" + re.escape(_NOTES_START) + r".*?" + re.escape(_NOTES_END) + r"\n?", re.S)
# The Notes card goes right after the page header, before the Reports card.
_HEADER_END_RE = re.compile(r"</header>\n?")


def inject_run_report_notes(run_dir: Path, notes: tuple[str, ...] | list[str]) -> str:
    """Write `notes` (from publish.yaml) into run_dir/report.html as a Notes
    card between marker comments. Idempotent: an earlier block is replaced,
    and an empty `notes` removes the block.

    report.html is a pipeline page and a sync of a new run overwrites it, so
    bin/sync_reports.sh calls this again (bin/study_pages.py inject-notes)
    before it records the run. run.json's report_sha256 is handled as in
    inject_run_report_links. Returns a one-line status."""
    report = run_dir / "report.html"
    if not report.is_file():
        return f"no {report}; nothing to annotate"
    old_bytes = report.read_bytes()
    old = old_bytes.decode("utf-8")
    text = _NOTES_BLOCK_RE.sub("", old)
    if notes:
        anchor = _HEADER_END_RE.search(text)
        if anchor is None:
            return f"WARNING: {report} has no </header>; notes not added"
        block = (
            f"  {_NOTES_START}\n"
            '  <section class="card">\n'
            '    <h2 class="card-title">Notes</h2>\n'
            + "".join(f'    <p class="card-note">{escape(n)}</p>\n' for n in notes)
            + "  </section>\n"
            f"  {_NOTES_END}\n"
        )
        text = text[:anchor.end()] + block + text[anchor.end():]
    if text == old:
        return f"{report}: notes already current"
    return f"{report}: notes set ({len(notes)})" + _commit_report(run_dir, report, old_bytes, text)


def _commit_report(run_dir: Path, report: Path, old_bytes: bytes, text: str) -> str:
    """Write the edited report. When run.json's report_sha256 matched the
    report before the edit, update it to the new content. When it did not
    match, leave run.json alone. Returns a status suffix for the caller."""
    new_bytes = text.encode("utf-8")
    report.write_bytes(new_bytes)
    status = ""
    meta_path = run_dir / "run.json"
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text())
        if meta.get("report_sha256") == _sha256_bytes(old_bytes):
            meta["report_sha256"] = _sha256_bytes(new_bytes)
            meta_path.write_text(json.dumps(meta, indent=2) + "\n")
            status += "; run.json report_sha256 updated"
        elif "report_sha256" in meta:
            status += "; WARNING: run.json report_sha256 did not match report.html before the edit, not changed"
    return status


def run_summary(r: dict) -> str:
    parts = []
    if r.get("n_families") is not None:
        parts.append(f"{r['n_families']:,} families")
    if r.get("n_strains") is not None:
        parts.append(f"{r['n_strains']:,} strains")
    if r.get("n_species") is not None:
        parts.append(f"{r['n_species']:,} species")
    return ", ".join(parts)


_RUNS_CSS = """
<style>
.crumbs { font-size: 0.85rem; margin-bottom: 1rem; }
.runs { border-collapse: collapse; margin-top: 1rem; font-size: 0.9rem; }
.runs th, .runs td { border: 1px solid #ccc; padding: 0.4rem 0.7rem; text-align: left;
                     vertical-align: top; }
.runs td.params { font-family: ui-monospace, monospace; font-size: 0.8rem;
                  overflow-wrap: anywhere; max-width: 32rem; }
.runs tr.current td { font-weight: 600; }
</style>
"""


def render_run_index(domain: str, study: str, runs: list[dict], current: str | None,
                     extras: dict[str, list[tuple[str, str]]] | None = None) -> str:
    """docs/<domain>/<study>/report.html: one row per run, current first,
    then the others newest first. extras maps a run name to its
    run_extra_links(); each one is linked after the run name."""
    extras = extras or {}
    body = [_RUNS_CSS]
    body.append(
        '<div class="crumbs"><a href="../../index.html">gallery</a> / '
        f'<a href="../index.html">{escape(domain)}</a> / {escape(study)}</div>'
    )
    body.append(f"<h1>{escape(study.replace('_', ' '))}</h1>")
    if not runs:
        body.append("<p>No runs published yet.</p>")
        return _page(study, "\n".join(body))
    if current:
        body.append(
            f'<p>Current result: <a href="{escape(current)}/report.html">{escape(current)}</a>. '
            "Other runs are archived versions (other parameters, dates, or pipeline versions).</p>"
        )
    ordered = sorted(runs, key=lambda r: (r.get("published") or "", r["run"]), reverse=True)
    ordered.sort(key=lambda r: r["run"] != current)
    body.append(
        '<table class="runs"><thead><tr><th>Run</th><th></th><th>Pipeline</th>'
        "<th>Summary</th><th>Parameters</th><th>Published</th></tr></thead><tbody>"
    )
    for r in ordered:
        is_current = r["run"] == current
        params = " ".join(r.get("params") or [])
        extra_links = "".join(
            f' &middot; <a href="{escape(r["run"])}/{escape(fname)}">{escape(label)}</a>'
            for fname, label in extras.get(r["run"], [])
        )
        body.append(
            f'<tr class="{"current" if is_current else ""}">'
            f'<td><a href="{escape(r["run"])}/report.html">{escape(r["run"])}</a>{extra_links}</td>'
            f'<td>{"current" if is_current else "archived"}</td>'
            f"<td>{escape(r.get('pipeline') or '')}</td>"
            f"<td>{escape(run_summary(r))}</td>"
            f'<td class="params">{escape(params)}</td>'
            f"<td>{escape(r.get('published') or '')}</td></tr>"
        )
    body.append("</tbody></table>")
    return _page(study, "\n".join(body))


def rebuild_study_pages(study_docs: Path, domain: str, study: str) -> list[dict]:
    """Rewrite the study-level report.html, index.html, and (when needed)
    study.json from the run folders on disk. If study.json names no run, or a
    run that no longer exists, and exactly one run exists, that run becomes
    current. With several runs and no valid current, no run is marked
    current and index.html points at the run list."""
    runs = load_runs(study_docs)
    names = {r["run"] for r in runs}
    current = read_current(study_docs)
    if current not in names:
        current = None
        if len(runs) == 1:
            current = runs[0]["run"]
            set_current(study_docs, current)
    extras = {r["run"]: run_extra_links(study_docs / r["run"]) for r in runs}
    (study_docs / "report.html").write_text(render_run_index(domain, study, runs, current, extras))
    target = f"{current}/report.html" if current else "report.html"
    (study_docs / "index.html").write_text(render_redirect(target))
    return runs
