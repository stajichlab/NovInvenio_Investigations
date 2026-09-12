"""Minimal, self-contained HTML templates for the NII docs/ gallery (DESIGN.md Sec 8).

Deliberately independent of nf_NovInvenio's lib/report_common.py / skins.py -- NII
should not depend on pipeline-repo code for its own publishing (that's the coupling
this whole split exists to avoid). Kept intentionally simple: card grids, no JS
framework, system fonts only, light+dark via prefers-color-scheme.
"""
from __future__ import annotations

import base64
from html import escape
from pathlib import Path

# Embedded as base64 data URIs rather than a relative path to assets/logo/:
# GitHub Pages (.github/workflows/static.yml) publishes only docs/, not the
# repo-root assets/ these files actually live in, so a relative link is broken
# on the live site regardless of path depth -- and depth itself varies (this
# module renders pages at docs/index.html, docs/<domain>/index.html, and,
# indirectly, docs/<domain>/<set>/*.html sit at yet another level). Same fix as
# NovInvenio's own report templates (lib/report_common.py's FAVICON_DATA_URI).
_ASSETS_LOGO_DIR = Path(__file__).resolve().parent.parent / "assets" / "logo"


def _data_uri(filename: str, mime: str) -> str:
    data = (_ASSETS_LOGO_DIR / filename).read_bytes()
    return f"data:{mime};base64,{base64.b64encode(data).decode('ascii')}"


_FAVICON_DATA_URI = _data_uri("NI_logo_favicon.ico", "image/x-icon")
_LOGO_DATA_URI = _data_uri("NI_logo_card-96.png", "image/png")

_BASE_CSS = """
:root { color-scheme: light dark; }
body { margin: 0; padding: 2rem; font-family: system-ui, -apple-system, sans-serif;
       background: #f7f7f5; color: #1a1a1a; }
@media (prefers-color-scheme: dark) {
  body { background: #16181c; color: #e6e6e6; }
  .card { background: #22252b !important; border-color: #33373f !important; }
  a { color: #7db4ff; }
}
h1 { font-size: 1.6rem; margin-bottom: 0.2rem; }
.title-row { display: flex; align-items: center; gap: 0.8rem; }
.title-row img.logo { width: 40px; height: 40px; border-radius: 8px; flex: 0 0 auto; }
.subtitle { color: #666; margin-bottom: 1.5rem; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 1rem; }
.card { border: 1px solid #ddd; border-radius: 10px; padding: 1rem 1.2rem;
        background: #fff; text-decoration: none; color: inherit; display: block; }
.card:hover { border-color: #888; }
.card.empty { opacity: 0.5; border-style: dashed; pointer-events: none; }
.card h2 { font-size: 1.1rem; margin: 0 0 0.3rem; overflow-wrap: anywhere; }
.card .desc { font-size: 0.9rem; color: #555; margin: 0 0 0.5rem; }
.card .meta { font-size: 0.8rem; color: #888; }
footer { margin-top: 2.5rem; font-size: 0.8rem; color: #888; }
.attribution { font-size: 0.75rem; color: #999; margin-top: 0.4rem; }
"""

_ATTRIBUTION = (
    "Includes data from UniProt (https://www.uniprot.org) and the Gene Ontology "
    "Consortium (http://geneontology.org), CC BY 4.0."
)


def _page(title: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<link rel="icon" href="{_FAVICON_DATA_URI}">
<style>{_BASE_CSS}</style>
</head>
<body>
{body}
<footer>{_ATTRIBUTION}</footer>
</body>
</html>
"""


def render_report_redirect() -> str:
    """A study's docs/<domain>/<set>/index.html -- so the bare directory URL
    (e.g. from a gallery card's own address bar, or someone guessing the
    folder path) lands on report.html instead of a static-host directory
    listing or 404. GitHub Pages has no server-side redirect config, so this
    is a client-side meta-refresh + JS fallback stub, not a rename of
    report.html itself -- report.html stays the canonical filename every
    other script (generate_docs.py's own study_status(), sync_reports.sh,
    the publish_*_release.sh scripts, nf_NovInvenio's make_index_report.py)
    already hardcodes."""
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=report.html">
<title>Redirecting…</title>
</head>
<body>
<p>Redirecting to <a href="report.html">report.html</a>…</p>
<script>location.replace("report.html");</script>
</body>
</html>
"""


def render_top_level(domains: list[dict], site_name: str = "NovInvenio Investigations") -> str:
    """domains: [{name, slug, desc, n_studies, status}]

    site_name: this repo's own display name (default matches NII's own,
    preserving today's output byte-for-byte) -- parameterized so a repo
    scaffolded from this one (nf_NovInvenio's bin/ni, issue #76) doesn't
    brand its gallery "NovInvenio Investigations" regardless of what it's
    actually called. bin/generate_docs.py derives it from pixi.toml's own
    `name` field rather than requiring a second place to edit it.
    """
    cards = []
    for d in domains:
        cls = "card" if d["n_studies"] else "card empty"
        href = f"{d['slug']}/index.html" if d["n_studies"] else "#"
        cards.append(f"""
<a class="{cls}" href="{escape(href)}">
  <h2>{escape(d['name'])}</h2>
  <p class="desc">{escape(d['desc'])}</p>
  <p class="meta">{d['n_studies']} study{'ies' if d['n_studies'] != 1 else ''}
     {'' if d['n_studies'] else '(not yet populated)'}</p>
</a>""")
    body = f"""
<div class="title-row"><img class="logo" src="{_LOGO_DATA_URI}" alt=""><h1>{escape(site_name)}</h1></div>
<p class="subtitle">Lineage-specific gene novelty/loss studies, organized by taxonomic domain.</p>
<div class="grid">{''.join(cards)}</div>
"""
    return _page(site_name, body)


_STATUS_TEXT = {
    "pending": " &mdash; not yet started (no data pulled)",
    "staged": " &mdash; data staged, pipeline not yet run",
    "complete": "",
}


def render_domain_index(
    domain_name: str, studies: list[dict], site_name: str = "NovInvenio Investigations"
) -> str:
    """studies: [{name, slug, hypothesis, n_ingroup, n_outgroup, status, updated}]

    site_name: see render_top_level()'s docstring -- same parameterization,
    same default.
    """
    cards = []
    for s in studies:
        complete = s["status"] == "complete"
        cls = "card" if complete else "card empty"
        href = f"{s['slug']}/report.html" if complete else "#"
        meta = (
            f"{s['n_ingroup']} ingroup &middot; {s['n_outgroup']} outgroup"
            if s.get("n_ingroup") is not None else ""
        )
        status_bit = _STATUS_TEXT.get(s["status"], "")
        cards.append(f"""
<a class="{cls}" href="{escape(href)}">
  <h2>{escape(s['name'].replace('_', ' '))}</h2>
  <p class="desc">{escape(s.get('hypothesis', ''))}</p>
  <p class="meta">{meta}{status_bit}</p>
</a>""")
    body = f"""
<div class="title-row"><img class="logo" src="{_LOGO_DATA_URI}" alt=""><h1>{escape(domain_name)}</h1></div>
<p class="subtitle"><a href="../index.html">&larr; all domains</a></p>
<div class="grid">{''.join(cards)}</div>
"""
    return _page(f"{domain_name} — {site_name}", body)
