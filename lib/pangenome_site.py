"""Site pages for pangenome.nf runs (DESIGN.md Sec 8, pangenome extension).

A pangenome study can have several runs (studies/<domain>/<set>/results/<run>/).
Each run gets its own page, and the study page lists the runs:

  docs/<domain>/<set>/report.html          study run list (lib/study_runs.py)
  docs/<domain>/<set>/<run>/report.html    run report, rendered from the
                                           pipeline's report/report.md (committed)
  docs/<domain>/<set>/<run>/run.json       run metadata + source checksum (committed)
  docs/<domain>/<set>/<run>/figures/       } release asset only, never committed
  docs/<domain>/<set>/<run>/figures_pdf/   } (bin/publish_report_release.sh,
  docs/<domain>/<set>/<run>/archive/       }  merged into docs/ at Pages-deploy
  docs/<domain>/<set>/<run>/island_synteny.html  time by static.yml)

Only report.html and run.json are committed. They are small text files. The
figures are regenerated binary files, and the tables and island_synteny.html
grow with family/island count, so a rerun must not recommit them.

The report is rendered from Markdown with markdown-it-py, with raw HTML
disabled: report.md text comes from pipeline output, and this module does not
control what an upstream table cell quotes.
"""
from __future__ import annotations

from html import escape

from markdown_it import MarkdownIt

from site_pages import _page

_REPORT_CSS = """
<style>
.report { max-width: 72rem; }
.report img { max-width: 100%; height: auto; background: #fff; border-radius: 6px; }
.report table { border-collapse: collapse; font-size: 0.85rem; display: block;
                overflow-x: auto; margin: 1rem 0; }
.report th, .report td { border: 1px solid #ccc; padding: 0.3rem 0.5rem;
                         text-align: left; vertical-align: top; }
.report td { overflow-wrap: anywhere; }
.crumbs { font-size: 0.85rem; margin-bottom: 1rem; }
.downloads li { margin: 0.2rem 0; }
</style>
"""


def render_markdown(md_text: str) -> str:
    """CommonMark + GFM tables, raw HTML disabled (escaped, not passed through)."""
    md = MarkdownIt("commonmark", {"html": False}).enable("table")
    return md.render(md_text)


def _crumbs(parts: list[tuple[str, str | None]]) -> str:
    items = []
    for label, href in parts:
        if href:
            items.append(f'<a href="{escape(href)}">{escape(label)}</a>')
        else:
            items.append(escape(label))
    return '<div class="crumbs">' + " / ".join(items) + "</div>"


def render_run_report(
    md_text: str,
    domain: str,
    set_name: str,
    run: str,
    downloads: list[tuple[str, str]],
    extra_pages: list[tuple[str, str]],
) -> str:
    """The per-run page: docs/<domain>/<set>/<run>/report.html.

    `downloads` / `extra_pages` are (label, relative href) pairs. The caller
    passes only files that exist, so every link on the page resolves once the
    release asset is merged in."""
    body = [_REPORT_CSS, '<div class="report">']
    body.append(_crumbs([
        ("gallery", "../../../index.html"),
        (domain, "../../index.html"),
        (set_name, "../report.html"),
        (run, None),
    ]))
    if extra_pages or downloads:
        body.append("<h2>Files</h2><ul class=\"downloads\">")
        for label, href in extra_pages:
            body.append(f'<li><a href="{escape(href)}">{escape(label)}</a></li>')
        for label, href in downloads:
            body.append(f'<li><a href="{escape(href)}">{escape(label)}</a> (tsv.gz)</li>')
        body.append("</ul>")
    body.append(render_markdown(md_text))
    body.append("</div>")
    return _page(f"{set_name} / {run}", "\n".join(body))
