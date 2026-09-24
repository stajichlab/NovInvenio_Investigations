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
    return PublishTarget(study=study, run=run)


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


def render_run_index(domain: str, study: str, runs: list[dict], current: str | None) -> str:
    """docs/<domain>/<study>/report.html: one row per run, current first,
    then the others newest first."""
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
        body.append(
            f'<tr class="{"current" if is_current else ""}">'
            f'<td><a href="{escape(r["run"])}/report.html">{escape(r["run"])}</a></td>'
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
    (study_docs / "report.html").write_text(render_run_index(domain, study, runs, current))
    target = f"{current}/report.html" if current else "report.html"
    (study_docs / "index.html").write_text(render_redirect(target))
    return runs
