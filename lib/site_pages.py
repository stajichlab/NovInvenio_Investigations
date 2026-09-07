"""Minimal, self-contained HTML templates for the NII docs/ gallery (DESIGN.md Sec 8).

Deliberately independent of nf_NovInvenio's lib/report_common.py / skins.py -- NII
should not depend on pipeline-repo code for its own publishing (that's the coupling
this whole split exists to avoid). Kept intentionally simple: card grids, no JS
framework, system fonts only, light+dark via prefers-color-scheme.
"""
from __future__ import annotations

from html import escape

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
.subtitle { color: #666; margin-bottom: 1.5rem; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 1rem; }
.card { border: 1px solid #ddd; border-radius: 10px; padding: 1rem 1.2rem;
        background: #fff; text-decoration: none; color: inherit; display: block; }
.card:hover { border-color: #888; }
.card.empty { opacity: 0.5; border-style: dashed; pointer-events: none; }
.card h2 { font-size: 1.1rem; margin: 0 0 0.3rem; }
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
<link rel="icon" href="../assets/logo/NI_logo_favicon.ico">
<style>{_BASE_CSS}</style>
</head>
<body>
{body}
<footer>{_ATTRIBUTION}</footer>
</body>
</html>
"""


def render_top_level(domains: list[dict]) -> str:
    """domains: [{name, slug, desc, n_studies, status}]"""
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
<h1>NovInvenio Investigations</h1>
<p class="subtitle">Lineage-specific gene novelty/loss studies, organized by taxonomic domain.</p>
<div class="grid">{''.join(cards)}</div>
"""
    return _page("NovInvenio Investigations", body)


_STATUS_TEXT = {
    "pending": " &mdash; not yet started (no data pulled)",
    "staged": " &mdash; data staged, pipeline not yet run",
    "complete": "",
}


def render_domain_index(domain_name: str, studies: list[dict]) -> str:
    """studies: [{name, slug, hypothesis, n_ingroup, n_outgroup, status, updated}]"""
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
  <h2>{escape(s['name'])}</h2>
  <p class="desc">{escape(s.get('hypothesis', ''))}</p>
  <p class="meta">{meta}{status_bit}</p>
</a>""")
    body = f"""
<h1>{escape(domain_name)}</h1>
<p class="subtitle"><a href="../index.html">&larr; all domains</a></p>
<div class="grid">{''.join(cards)}</div>
"""
    return _page(f"{domain_name} — NovInvenio Investigations", body)
