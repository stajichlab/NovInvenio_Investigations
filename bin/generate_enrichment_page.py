#!/usr/bin/env python3
"""Write a static functional-enrichment page into a published run folder.

Input (written by the study's enrichment analysis, not by this script):
  studies/<domain>/<folder>/analysis/enrichment/outputs/<sub>/proteome_clustered/
      tests.tsv  coverage.tsv  enrichment_domains.tsv.gz  enrichment_go.tsv.gz
      PROVENANCE.yaml  controls.tsv (optional)
  studies/<domain>/<folder>/analysis/enrichment/figures/<folder><suffix>_summary.png
where <sub> is "" for --direction novelty and "loss/" for --direction loss,
and <suffix> is "" for novelty and "_loss" for loss.

Output, in the run folder docs/<domain>/<study>/<run>/ (from publish.yaml):
  enrichment.html / enrichment_loss.html          the page
  enrichment_summary.png / enrichment_loss_summary.png   the figure, scaled down
  enrichment_provenance.yaml / enrichment_loss_provenance.yaml   copy of PROVENANCE.yaml
Then the run report.html gets a link to each enrichment page that exists
(lib/study_runs.inject_run_report_links). The gallery card and the study run
list link the pages when bin/generate_docs.py runs next (file-exists check).

All outputs are class 2 (CLAUDE.md): their size does not grow with the
candidate count. The page holds summary tables and top hits only, never
candidate protein lists.

If the outputs of a direction do not exist, the script writes nothing for
that direction and exits 0.

Usage:
  bin/generate_enrichment_page.py --study fungi/yeast_filamentous [--direction loss]
  bin/generate_enrichment_page.py --links-only --docs-run-dir docs/fungi/x/run
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shlex
import shutil
import sys
from html import escape
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "lib"))
from study_runs import inject_run_report_links, load_publish_target, vocab_meta_tag  # noqa: E402
from te_flags import te_class  # noqa: E402

# candidate_ids cells in the enrichment tables can exceed csv's default 128 KB field limit.
csv.field_size_limit(sys.maxsize)

DIRECTIONS = {
    # direction: (outputs subdir, figure suffix, page name, figure name, provenance copy, text)
    "novelty": ("", "", "enrichment.html", "enrichment_summary.png", "enrichment_provenance.yaml",
                "novelty candidates: gene families present in the IN group and absent from the OUT group"),
    "loss": ("loss/", "_loss", "enrichment_loss.html", "enrichment_loss_summary.png",
             "enrichment_loss_provenance.yaml",
             "loss candidates: gene families present in the OUT group and absent from the IN group"),
}
TOP_DOMAINS = 25
TOP_GO = 15
FIG_MAX_WIDTH = 1600
FDR_CUT = 0.05
LOW_ANNOT = 0.5

# Columns this script reads. They are checked against each file header;
# a missing column stops the script.
DOMAIN_COLS = ["scope", "level", "vocab", "term", "name", "k_cand_with_term", "N_cand_annotated",
               "n_bg_with_term", "M_bg_annotated", "fold", "pvalue", "fdr_bh"]
GO_COLS = DOMAIN_COLS + ["go_evidence", "n_cand_units_non_IEA"]
TESTS_COLS = ["scope", "level", "vocab", "role", "n_cand_units", "n_terms_tested", "n_fdr_lt_0.05", "status"]
COVERAGE_COLS = ["species", "subset", "n_proteins", "frac_pfam", "frac_pfam_clan", "frac_interpro",
                 "frac_go", "frac_any"]

TERM_URL = {
    "pfam": (re.compile(r"^PF\d{5}$"), "https://www.ebi.ac.uk/interpro/entry/pfam/{}/"),
    "interpro": (re.compile(r"^IPR\d{6}$"), "https://www.ebi.ac.uk/interpro/entry/InterPro/{}/"),
    "go": (re.compile(r"^GO:\d{7}$"), "https://www.ebi.ac.uk/QuickGO/term/{}"),
}


# --- reading -----------------------------------------------------------------


def read_tsv(path: Path, need: list[str]) -> list[dict]:
    """Read a TSV (plain or .gz) with a header row. Stop when a needed column is missing."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in need if c not in header]
        if missing:
            raise SystemExit(f"ERROR: {path} has no column(s) {missing}; header is {header}")
        # Keep only the needed columns; candidate_ids cells can be very long.
        return [{c: row[c] for c in need} for row in reader]


def read_controls(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        # Drop the protein ID columns (candidate_ids, near_miss_*_ids): the page
        # holds no candidate protein lists. pfam_ids holds Pfam accessions; keep it.
        cols = [c for c in (reader.fieldnames or []) if c == "pfam_ids" or not c.endswith("_ids")]
    return cols, rows


def fnum(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def results_dir_of(prov: dict, root: Path) -> str | None:
    """The --results argument of the recorded command, relative to the repo
    root. When the path is not under root (e.g. a copy of the repo), the
    part from the last 'results' component on is used."""
    cmd = prov.get("command")
    if not isinstance(cmd, str):
        return None
    toks = shlex.split(cmd)
    for i, tok in enumerate(toks[:-1]):
        if tok == "--results":
            p = Path(toks[i + 1])
            try:
                return str(p.relative_to(root.resolve()))
            except ValueError:
                parts = p.parts
                if "results" in parts:
                    j = len(parts) - 1 - parts[::-1].index("results")
                    return str(Path(*parts[j:]))
                return str(p)
    return None


# --- rendering helpers ---------------------------------------------------------


def fmt_frac(v: str) -> str:
    x = fnum(v)
    return "" if x != x else f"{100 * x:.1f}%"


def fmt_fdr(v: str) -> str:
    x = fnum(v)
    if x != x:
        return ""
    return f"{x:.2e}" if x < 0.001 else f"{x:.3f}"


def term_cell(vocab: str, term: str) -> str:
    rx, url = TERM_URL.get(vocab, (None, None))
    if rx is not None and rx.match(term):
        return f'<a href="{escape(url.format(term))}">{escape(term)}</a>'
    return escape(term)


def table(headers: list[str], rows: list[list[str]], cls: str = "") -> str:
    """rows hold HTML that is already escaped."""
    head = "".join(f"<th>{escape(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return (f'<div class="scroll"><table class="{cls}"><thead><tr>{head}</tr></thead>'
            f"<tbody>{body}</tbody></table></div>")


def top_hits(rows: list[dict], vocab: str, n: int) -> list[dict]:
    sel = [r for r in rows if r["scope"] == "whole_set" and r["level"] == "family" and r["vocab"] == vocab]
    sel.sort(key=lambda r: (fnum(r["fdr_bh"]), fnum(r["pvalue"]), r["term"]))
    return sel[:n]


def hits_table(hits: list[dict], vocab: str, go: bool = False, unavailable: str | None = None) -> str:
    """unavailable: vocabulary label when no protein of the study has a term of it."""
    if unavailable:
        return (f'<p class="note"><b>Not available</b>: no protein of the candidate source group has '
                f"a {escape(unavailable)} annotation in this study (for example Pfam-only HMM "
                f"annotation). No {escape(unavailable)} test was run.</p>")
    if not hits:
        return '<p class="note">No whole-set family-level rows for this vocabulary.</p>'
    headers = ["Term", "Name", "k / N cand.", "n / M bg", "Fold", "FDR (BH)"]
    if go:
        headers += ["GO evidence", "Families with non-IEA"]
    rows = []
    for r in hits:
        sig = fnum(r["fdr_bh"]) < FDR_CUT
        tag = (' <span class="tag-te" title="transposable-element family by name (heuristic)">TE</span>'
               if vocab in ("pfam", "pfam_clan") and te_class(r["name"]) == "TE" else "")
        cells = [
            term_cell(vocab, r["term"]),
            escape(r["name"]) + tag,
            f'{escape(r["k_cand_with_term"])} / {escape(r["N_cand_annotated"])}',
            f'{escape(r["n_bg_with_term"])} / {escape(r["M_bg_annotated"])}',
            escape(r["fold"]),
            ("<b>" if sig else "") + fmt_fdr(r["fdr_bh"]) + ("</b>" if sig else ""),
        ]
        if go:
            cells += [escape(r["go_evidence"]), escape(r["n_cand_units_non_IEA"])]
        rows.append(cells)
    n_sig = sum(fnum(r["fdr_bh"]) < FDR_CUT for r in hits)
    n_te = sum(1 for r in hits if vocab in ("pfam", "pfam_clan") and te_class(r["name"]) == "TE")
    te_note = (f" {n_te} are tagged <span class=\"tag-te\">TE</span>: transposable-element families, "
               "flagged by name (heuristic, not a test)." if n_te else "")
    return (f'<p class="note">{len(hits)} rows shown, ordered by FDR then p-value; '
            f"{n_sig} have FDR &lt; {FDR_CUT} (bold).{te_note}</p>" + table(headers, rows, "hits"))


VOCABS = [("pfam", "Pfam"), ("pfam_clan", "Pfam clan"), ("interpro", "InterPro"), ("go", "GO")]


def tests_table(tests: list[dict]) -> str:
    """One row per scope x level. Per vocabulary: terms with FDR < 0.05 / terms tested."""
    order: list[tuple[str, str]] = []
    by_key: dict[tuple[str, str], dict[str, dict]] = {}
    for r in tests:
        key = (r["scope"], r["level"])
        if key not in by_key:
            order.append(key)
            by_key[key] = {}
        by_key[key][r["vocab"]] = r
    order.sort(key=lambda k: (k[0] != "whole_set", k[0], k[1] != "family", k[1]))
    rows = []
    for key in order:
        cell = by_key[key]
        n_units = next(iter(cell.values()))["n_cand_units"]
        line = [escape(key[0]), escape(key[1]), escape(n_units)]
        for vocab, _ in VOCABS:
            r = cell.get(vocab)
            if r is None:
                line.append("")
                continue
            txt = f'{escape(r["n_fdr_lt_0.05"])} / {escape(r["n_terms_tested"])}'
            if r["n_fdr_lt_0.05"] not in ("0", ""):
                txt = f"<b>{txt}</b>"
            if r["status"] != "ok":
                txt += f' <span class="warn">({escape(r["status"])})</span>'
            line.append(txt)
        rows.append(line)
    headers = ["Scope", "Level", "Cand. units"] + [f"{label} sig / tested" for _, label in VOCABS]
    return table(headers, rows, "tests")


def coverage_table(cov: list[dict]) -> str:
    rows = []
    for r in cov:
        rows.append([escape(r["species"]), escape(r["subset"]), escape(r["n_proteins"])]
                    + [fmt_frac(r[c]) for c in ("frac_pfam", "frac_pfam_clan", "frac_interpro",
                                                  "frac_go", "frac_any")])
    return table(["Species", "Set", "Proteins", "Pfam", "Pfam clan", "InterPro", "GO", "Any"], rows, "cov")


def low_annotation(cov: list[dict], skip: set[str] | frozenset[str] = frozenset()) -> list[str]:
    """Whole-set candidate annotation fractions below LOW_ANNOT. `skip`: vocabularies
    (pfam, interpro, go) to leave out, e.g. ones that the page already lists as
    not available."""
    out = []
    for r in cov:
        if r["subset"] == "candidates" and r["species"].startswith("ALL_"):
            for col, label in (("frac_pfam", "Pfam"), ("frac_interpro", "InterPro"),
                               ("frac_go", "GO"), ("frac_any", "any vocabulary")):
                if col[len("frac_"):] in skip:
                    continue
                x = fnum(r[col])
                if x == x and x < LOW_ANNOT:
                    out.append(f"{label} {100 * x:.1f}%")
    return out


GAP_VOCABS = [("go", "GO"), ("interpro", "InterPro"), ("pfam", "Pfam"), ("pfam_clan", "Pfam clan")]


def vocabulary_gaps(cov: list[dict]) -> list[str]:
    """Caveat lines for vocabularies that some species have no term for.

    A species lacks a vocabulary when its universe row has proteins and a zero fraction
    for that vocabulary (for example a Pfam-only HMM annotation has no GO or InterPro).
    Vocabularies with the same covered species share one line. No gap: empty list, and
    the page is the same as before this check existed."""
    sp_rows = [r for r in cov if r["subset"] == "universe" and not r["species"].startswith("ALL_")
               and fnum(r["n_proteins"]) > 0]
    groups: dict[tuple[str, ...], list[str]] = {}
    for vocab, label in GAP_VOCABS:
        covered = tuple(r["species"] for r in sp_rows if fnum(r[f"frac_{vocab}"]) > 0)
        if len(covered) < len(sp_rows):
            groups.setdefault(covered, []).append(label)
    lines = []
    for covered, labels in groups.items():
        names = labels[0] if len(labels) == 1 else ", ".join(labels[:-1]) + " and " + labels[-1]
        missing = [r["species"] for r in sp_rows if r["species"] not in covered]
        if not covered:
            lines.append(
                f"{escape(names)}: not available for this study. No species has "
                f"{escape(names)} annotation ({escape(', '.join(missing))}), so no "
                f"{escape(names)} test was run.")
            continue
        lines.append(
            f"{escape(names)} {'covers' if len(labels) == 1 else 'cover'} only {len(covered)} of "
            f"{len(sp_rows)} species: "
            f"{escape(', '.join(covered)) if covered else 'none'}. "
            f"No {escape(names)} annotation for {escape(', '.join(missing))}; "
            f"{'this species adds' if len(missing) == 1 else 'these species add'} no units "
            f"to the {escape(names)} tests, so those results describe the covered species only.")
    return lines


def unavailable_vocabs(cov: list[dict]) -> set[str]:
    """Vocabularies with no term in any protein of the candidate source group (the ALL_*
    universe row has fraction 0). Empty set when there is no ALL_* universe row."""
    rows = [r for r in cov if r["subset"] == "universe" and r["species"].startswith("ALL_")]
    if not rows:
        return set()
    return {v for v, _ in GAP_VOCABS if all(fnum(r[f"frac_{v}"]) == 0 for r in rows)}


def available_vocabs(cov: list[dict]) -> list[str]:
    """pfam / interpro / go that have a term in the candidate source group (not in
    unavailable_vocabs). Written into the page as a <meta> tag; lib/study_runs.py reads
    it for the link text."""
    na = unavailable_vocabs(cov)
    return [v for v in ("pfam", "interpro", "go") if v not in na]


def all_sources_hmm(prov: dict) -> bool:
    """True when PROVENANCE.yaml annotation_sources names a Pfam HMM search for every file."""
    src = prov.get("annotation_sources") or {}
    return bool(src) and all("HMM search" in str(v) for v in src.values())


def controls_table(cols: list[str], rows: list[dict]) -> str:
    return table(cols, [[escape(r.get(c, "")) for c in cols] for r in rows], "ctl")


def scale_figure(src: Path, dst: Path) -> str:
    """Write src to dst as a small PNG: at most FIG_MAX_WIDTH px wide, flat on
    white, 256-colour palette. The source figures are about 2,850 px wide RGBA
    (14 in at 200 dpi, panels stacked in rows); this gives about 0.1-0.25 MB.
    The full-size figure stays in the study's analysis/enrichment/figures/.
    FASTOCTREE, not MEDIANCUT: median cut merged small legend patches into
    other colours (the #7b5aa6 purple came out about 43 RGB units off on
    the chytrid_vs_mucoromycota loss figure; octree: under 9)."""
    from PIL import Image
    with Image.open(src) as im:
        im = im.convert("RGBA")
        w, h = im.size
        if w > FIG_MAX_WIDTH:
            im = im.resize((FIG_MAX_WIDTH, round(h * FIG_MAX_WIDTH / w)), Image.LANCZOS)
        flat = Image.new("RGB", im.size, (255, 255, 255))
        flat.paste(im, mask=im.split()[3])
        pal = flat.quantize(colors=256, method=Image.Quantize.FASTOCTREE)
        buf = io.BytesIO()
        pal.save(buf, format="PNG", optimize=True)
    dst.write_bytes(buf.getvalue())
    return f"{w} px wide source, written {min(w, FIG_MAX_WIDTH)} px wide, {len(buf.getvalue()):,} bytes"


CSS = """
:root { color-scheme: light dark;
  --bg: #f7f7f5; --fg: #1a1a1a; --muted: #666; --line: #d8d8d4; --card: #ffffff;
  --link: #0b5cad; --warn-bg: #fff6e0; --warn-line: #e0b84c; --head: #efefeb; }
@media (prefers-color-scheme: dark) { :root {
  --bg: #16181c; --fg: #e6e6e6; --muted: #a0a4ab; --line: #33373f; --card: #22252b;
  --link: #7db4ff; --warn-bg: #2e2716; --warn-line: #8a6d1f; --head: #2a2e35; } }
body { margin: 0; padding: 1.5rem 1rem 3rem; background: var(--bg); color: var(--fg);
       font-family: system-ui, -apple-system, sans-serif; line-height: 1.45; }
main { max-width: 1100px; margin: 0 auto; }
a { color: var(--link); }
h1 { font-size: 1.5rem; margin: 0.2rem 0; }
h2 { font-size: 1.15rem; margin: 2rem 0 0.5rem; }
.crumbs, .note, .method { font-size: 0.88rem; color: var(--muted); }
.method { color: var(--fg); }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; font-size: 0.82rem; background: var(--card); }
th, td { border: 1px solid var(--line); padding: 0.25rem 0.5rem; text-align: left; vertical-align: top; }
th { background: var(--head); }
td { font-variant-numeric: tabular-nums; }
table.hits td:nth-child(2) { min-width: 18rem; }
.caveats { background: var(--warn-bg); border: 1px solid var(--warn-line); border-radius: 8px;
           padding: 0.6rem 1rem; margin: 1rem 0; font-size: 0.9rem; }
.caveats ul { margin: 0.3rem 0; padding-left: 1.2rem; }
.warn { color: #b06000; }
.tag-te { background: #f3d9b8; color: #5a3200; border-radius: 3px; padding: 0 0.3em; font-size: 0.75rem; font-weight: 600; }
@media (prefers-color-scheme: dark) { .tag-te { background: #5a3a14; color: #f3d9b8; } }
img.fig { max-width: 100%; height: auto; background: #fff; border: 1px solid var(--line); border-radius: 6px; }
"""


def render(study: str, direction: str, run: str, prov: dict, cov: list[dict], tests: list[dict],
           domains: list[dict], go: list[dict], controls: tuple[list[str], list[dict]] | None,
           fig_name: str | None, prov_name: str, src_rel: str) -> str:
    _, _, _, _, _, dir_text = DIRECTIONS[direction]
    params = prov.get("parameters") or {}
    bgc = prov.get("background_clustering") or {}
    method = (
        f"One-sided hypergeometric test per term; Benjamini-Hochberg FDR within each scope, level "
        f"and vocabulary. Candidates are the {escape(dir_text)}. Universe: "
        f"{escape(str(prov.get('universe', 'not recorded')))}. Background clustered: "
        f"{escape(str(params.get('cluster_background', 'not recorded')))}"
        + (f" ({escape(str(bgc.get('n_sequences')))} sequences in "
           f"{escape(str(bgc.get('n_clusters')))} clusters, mmseqs easy-cluster --min-seq-id 0.3 -c 0.8)"
           if bgc.get("n_clusters") else "")
        + ". Family level: one unit per candidate family and per background cluster."
    )
    na = unavailable_vocabs(cov)
    # GO and InterPro both not available: the vocabulary-gap line already says so;
    # their 0% would only repeat it in the low-annotation line.
    low = low_annotation(cov, skip={"go", "interpro"} if {"go", "interpro"} <= na else set())
    caveats = [
        ("Annotation: Pfam-A HMM search (hmmsearch, Pfam gathering thresholds) for every species; "
         "no UniProt annotation. Each annotation file has a provenance sidecar (domain table, HMMER "
         "version, Pfam release); PROVENANCE.yaml lists the files and their checksums."
         if all_sources_hmm(prov) else
         "Annotation provenance is not fully recorded: PROVENANCE.yaml lists the per-species "
         "annotation files and their checksums, but not the UniProt release they come from."),
        ("GO: not available for this study (no GO annotation in any species)." if "go" in na else
         "GO results are secondary. GO annotations include computational evidence (IEA, IBA); "
         "the GO evidence column shows if a term has non-IEA support."),
        "Candidate families (pipeline clustering of candidates) and background clusters (mmseqs "
        "clustering of the background proteome) are made by different runs, so family-level counts "
        "are not fully comparable.",
    ]
    if low:
        caveats.append("Low annotation: the candidate set has less than 50% annotation for "
                       + ", ".join(low) + ". Enrichment for these vocabularies has low power.")
    caveats.extend(vocabulary_gaps(cov))
    pfam_rel = (prov.get("pfam_dat") or {}).get("path", "")
    go_ver = (prov.get("go_obo") or {}).get("data_version", "not recorded")
    ipr = prov.get("interpro_entries", "")
    sources = (f"Pfam: {escape(str(Path(str(pfam_rel)).parent.name) or 'not recorded')}; "
               f"InterPro entries: {escape(str(Path(str(ipr)).parent.name) or 'not recorded')}; "
               f"GO: {escape(str(go_ver))}; analysis generated {escape(str(prov.get('generated', '')))}.")

    parts = [f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
{vocab_meta_tag(available_vocabs(cov))}
<title>{escape(study)} enrichment ({escape(direction)})</title>
<style>{CSS}</style>
</head>
<body>
<main>
<div class="crumbs"><a href="../../../index.html">gallery</a> / <a href="../../index.html">domain</a> /
<a href="../report.html">{escape(study)}</a> / <a href="report.html">{escape(run)}</a> / enrichment ({escape(direction)})</div>
<h1>{escape(study.replace('_', ' '))}: functional enrichment ({escape(direction)})</h1>
<p class="method">{method}</p>
<p class="note">Sources: {sources} Source tables: <code>{escape(src_rel)}</code>.
Full record: <a href="{escape(prov_name)}">PROVENANCE</a>. Run report: <a href="report.html">report.html</a>.</p>
<div class="caveats"><b>Caveats</b><ul>{''.join(f'<li>{c}</li>' for c in caveats)}</ul></div>
"""]
    if fig_name:
        parts.append(f'<h2>Summary figure</h2><img class="fig" src="{escape(fig_name)}" '
                     f'alt="Enrichment summary figure for {escape(study)} ({escape(direction)})">')
    parts.append("<h2>Annotation coverage</h2>"
                 '<p class="note">Fraction of proteins with at least one term per vocabulary. '
                 "Universe = all proteins of the species; candidates = the candidate set.</p>"
                 + coverage_table(cov))
    parts.append("<h2>Tests summary</h2>"
                 f'<p class="note">Per vocabulary: terms with FDR &lt; {FDR_CUT} / terms tested. '
                 "Protein level counts proteins; family level counts candidate families and background "
                 "clusters.</p>" + tests_table(tests))
    parts.append(f"<h2>Top {TOP_DOMAINS} Pfam terms (whole set, family level)</h2>"
                 + hits_table(top_hits(domains, "pfam", TOP_DOMAINS), "pfam",
                              unavailable="Pfam" if "pfam" in na else None))
    parts.append(f"<h2>Top {TOP_DOMAINS} InterPro terms (whole set, family level)</h2>"
                 + hits_table(top_hits(domains, "interpro", TOP_DOMAINS), "interpro",
                              unavailable="InterPro" if "interpro" in na else None))
    parts.append(f"<h2>Top {TOP_GO} GO terms (whole set, family level) &mdash; secondary, "
                 "computational evidence (IEA/IBA) included</h2>"
                 + hits_table(top_hits(go, "go", TOP_GO), "go", go=True,
                              unavailable="GO" if "go" in na else None))
    if controls is not None:
        cols, rows = controls
        parts.append("<h2>Actin cytoskeleton controls</h2>"
                     '<p class="note">Known actin-related families and how many of their proteins pass '
                     "each filter. Protein ID columns are left out.</p>" + controls_table(cols, rows))
    # No write date on the page: the same inputs give the same bytes.
    parts.append('<p class="note">Page written by bin/generate_enrichment_page.py.</p>\n'
                 "</main>\n</body>\n</html>\n")
    return "\n".join(parts)


# --- main --------------------------------------------------------------------


def split2(value: str) -> tuple[str, str]:
    parts = value.split("/")
    if len(parts) != 2 or not all(parts):
        raise SystemExit(f"ERROR: --study must be <domain>/<folder>, got {value!r}")
    return parts[0], parts[1]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", help="<domain>/<folder> under studies/")
    ap.add_argument("--direction", choices=sorted(DIRECTIONS), default="novelty")
    ap.add_argument("--docs-run-dir", type=Path,
                    help="published run folder (default: docs/<domain>/<study>/<run>/ from publish.yaml)")
    ap.add_argument("--out", type=Path, help="page path (default: <docs-run-dir>/enrichment[_loss].html)")
    ap.add_argument("--enrichment-dir", type=Path,
                    help="default: studies/<domain>/<folder>/analysis/enrichment")
    ap.add_argument("--figure", type=Path, help="default: <enrichment-dir>/figures/<folder>[_loss]_summary.png")
    ap.add_argument("--links-only", action="store_true",
                    help="only update the links in <docs-run-dir>/report.html (for bin/sync_reports.sh)")
    ap.add_argument("--no-links", action="store_true", help="do not edit report.html")
    ap.add_argument("--allow-source-mismatch", action="store_true",
                    help="write the page even when run.json source_dir differs from the --results "
                         "of the enrichment command")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    args = ap.parse_args(argv)
    root = args.repo_root

    if args.links_only:
        if args.docs_run_dir is None:
            raise SystemExit("ERROR: --links-only needs --docs-run-dir")
        print(inject_run_report_links(args.docs_run_dir), file=sys.stderr)
        return 0
    if not args.study:
        raise SystemExit("ERROR: --study is required")

    domain, folder = split2(args.study)
    sub, fig_suffix, page_name, fig_name, prov_name, _ = DIRECTIONS[args.direction]
    enr = args.enrichment_dir or root / "studies" / domain / folder / "analysis" / "enrichment"
    src = enr / "outputs" / sub / "proteome_clustered"
    needed = ["tests.tsv", "coverage.tsv", "enrichment_domains.tsv.gz", "enrichment_go.tsv.gz", "PROVENANCE.yaml"]
    missing = [n for n in needed if not (src / n).is_file()]
    if missing:
        print(f"== {args.study} {args.direction}: {src} lacks {missing}; no page written ==", file=sys.stderr)
        return 0

    run_dir = args.docs_run_dir
    if run_dir is None:
        target = load_publish_target(root / "studies" / domain / folder)
        if target is None or target.run is None:
            print(f"== {args.study}: no publish.yaml with a run; not published, no page written ==",
                  file=sys.stderr)
            return 0
        run_dir = root / "docs" / domain / target.study / target.run
    run_meta = run_dir / "run.json"
    if not run_meta.is_file():
        print(f"== {run_dir}/run.json not found; run not published, no page written ==", file=sys.stderr)
        return 0
    study_name = run_dir.parent.name

    prov = yaml.safe_load((src / "PROVENANCE.yaml").read_text()) or {}
    source_dir = json.loads(run_meta.read_text()).get("source_dir")
    results = results_dir_of(prov, root)
    if source_dir and results and source_dir != results:
        msg = f"run.json source_dir {source_dir!r} differs from enrichment --results {results!r}"
        if not args.allow_source_mismatch:
            raise SystemExit(f"ERROR: {msg}; use --allow-source-mismatch to write the page anyway")
        print(f"WARNING: {msg}", file=sys.stderr)

    cov = read_tsv(src / "coverage.tsv", COVERAGE_COLS)
    tests = read_tsv(src / "tests.tsv", TESTS_COLS)
    domains = read_tsv(src / "enrichment_domains.tsv.gz", DOMAIN_COLS)
    go = read_tsv(src / "enrichment_go.tsv.gz", GO_COLS)
    controls = read_controls(src / "controls.tsv") if (src / "controls.tsv").is_file() else None

    fig_src = args.figure or enr / "figures" / f"{folder}{fig_suffix}_summary.png"
    fig_out = None
    if fig_src.is_file():
        print(f"figure: {scale_figure(fig_src, run_dir / fig_name)}", file=sys.stderr)
        fig_out = fig_name
    else:
        print(f"WARNING: {fig_src} not found; page has no figure", file=sys.stderr)
    shutil.copyfile(src / "PROVENANCE.yaml", run_dir / prov_name)

    try:
        src_rel = str(src.resolve().relative_to(root.resolve()))
    except ValueError:
        src_rel = str(src)
    html = render(study_name, args.direction, run_dir.name, prov, cov, tests, domains, go, controls,
                  fig_out, prov_name, src_rel)
    out = args.out or run_dir / page_name
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out} ({len(html.encode()):,} bytes; sha256 "
          f"{hashlib.sha256(html.encode()).hexdigest()[:12]})", file=sys.stderr)
    if not args.no_links:
        print(inject_run_report_links(run_dir), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
