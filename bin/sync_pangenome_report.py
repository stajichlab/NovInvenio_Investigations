#!/usr/bin/env python3
"""Stage one pangenome.nf run's report into the docs/ gallery.

pangenome.nf writes report/report.md + figures, not the novelties/core/losses
HTML set that bin/sync_reports.sh handles. This script is the pangenome
counterpart. Layout and data classes: see lib/pangenome_site.py's docstring.

The folder must have a publish.yaml naming its study (spec:
notes/superpowers/specs/2026-09-24-study-run-site-layout-design.md); without
one the run is not published. The run name is --run (one folder holds several
runs under results/).

Steps, for --study <domain>/<folder> --run <run>:
  1. Find the run's pangenome output dir:
     studies/<domain>/<set>/results/<run>/*/pangenome/ (the `*` is the
     pipeline's --outdir project subdir). Exactly one match is required;
     pass --pangenome-dir to override.
  2. Replace docs/<domain>/<set>/<run>/{figures,figures_pdf,archive,clinker}/
     and island_synteny.html / assembly_quality.html with this run's copies.
     clinker/ holds the per-locus clinker pages island_synteny.html loads in
     its "Synteny (clinker)" panel (NovInvenio spec 2026-09-24 section 8).
     These are gitignored (release asset only). Only the top
     --clinker_publish_top loci' pages (default 25; 0 = none; -1 = all) are
     staged into clinker/, in ascending L### order; the rest stay in the
     pipeline's own clinker/ dir (never copied anywhere in NII) and
     island_synteny.html gets a small injected script
     (window.CLINKER_PUBLISHED / window.CLINKER_FULL_DIR) so an
     unpublished locus shows a note pointing at that full folder instead
     of a broken iframe (N3, NovInvenio#202).
     Optional pages, only when their dir is passed (never auto-detected):
     --locus_summary_dir -> locus_summary.html + locus_summary/;
     --starship_overlap_dir (+ --starship_sm_dir) -> starship_overlap.html +
     starship_overlap/ (stage_starship_overlap()). Both are release assets.
  3. Render docs/<domain>/<set>/<run>/report.html from report/report.md, and
     write run.json (source path, sha256 of report.md, counts, and, when
     clinker/ exists, clinker_published/clinker_total/clinker_full_dir).
  4. Rebuild the study pages (run list, index redirect, study.json) from
     every docs/<domain>/<study>/*/run.json. --current marks this run current.

It does not publish anything. Run bin/generate_docs.py after it to refresh the
gallery, then bin/publish_report_release.sh <domain>/<set>/<run> to upload the
release asset.

Usage:
  bin/sync_pangenome_report.py --study fungi/coccidioides_pangenome \\
      --run rescue_structural_genus_vs_ureesii
"""
from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from pangenome_site import render_markdown, render_run_report  # noqa: E402
from site_pages import _page, render_report_redirect  # noqa: E402
from study_runs import load_publish_target, rebuild_study_pages, set_current, validate_slug  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent

# (path relative to the pangenome dir, archive file name). Fixed allowlist:
# every entry is an aggregate or per-island table. Per-gene and per-pair
# tables (gene_positions, cooccurring_pairs, pair_classification, the
# presence matrices) are hundreds of MB to GB and are not published.
ARCHIVE_TABLES = [
    ("frequency_table.tsv", "frequency_table.tsv.gz"),
    ("clade_assignments.tsv", "clade_assignments.tsv.gz"),
    ("strain_inventory.tsv", "strain_inventory.tsv.gz"),
    ("significant_islands.tsv", "significant_islands.tsv.gz"),
    ("island_pfam_enrichment.tsv", "island_pfam_enrichment.tsv.gz"),
    ("report/pangenome_openness.tsv", "pangenome_openness.tsv.gz"),
    ("report_tables/per_strain_summary.tsv", "per_strain_summary.tsv.gz"),
    ("report_tables/classification_counts.tsv", "classification_counts.tsv.gz"),
    ("report_tables/island_size_distribution.tsv", "island_size_distribution.tsv.gz"),
    ("report_tables/marker_summary.tsv", "marker_summary.tsv.gz"),
    ("report_tables/islands_with_domains.tsv", "islands_with_domains.tsv.gz"),
    ("assembly_quality_vs_content.tsv", "assembly_quality_vs_content.tsv.gz"),
    ("assembly_quality_correlations.tsv", "assembly_quality_correlations.tsv.gz"),
    ("diagnostics/diagnostics.tsv", "diagnostics.tsv.gz"),
    ("report_tables/group_class_overlap.tsv", "group_class_overlap.tsv.gz"),  # NovInvenio #212 PR 2
]

# Directories and pages copied or rendered into the run dir. All gitignored.
# clinker/ is staged separately (stage_clinker()): only the top
# --clinker_publish_top loci' pages are copied, not the whole directory.
ASSET_DIRS = [("report/figures", "figures"), ("report/figures_pdf", "figures_pdf")]
RELEASE_ONLY = ["figures", "figures_pdf", "archive", "clinker", "island_synteny.html",
                "assembly_quality.html", "locus_summary", "locus_summary.html",
                "starship_overlap", "starship_overlap.html"]

# Files of a bin/starship_island_overlap.py output dir (--starship_overlap_dir)
# and of an SM table dir (--starship_sm_dir, default: the same dir) that are
# copied into run_docs/starship_overlap/. Small aggregate tables plus the
# per-island/per-cluster gz tables (< 1 MB each for the runs of 2026-10-17).
STARSHIP_OVERLAP_FILES = ["summary.tsv", "enrichment.tsv", "per_starship.tsv", "class_breakdown.tsv",
                          "coord_permutation.tsv", "per_island.tsv.gz", "annotated_families.tsv.gz",
                          "PROVENANCE.yaml"]
STARSHIP_SM_FILES = ["sm_by_type.tsv", "sm_starship_overlap.tsv", "sm_clusters.tsv.gz",
                     "sm_permutation.tsv", "sm_family_class_direction.tsv", "sm_core_families.tsv.gz",
                     "PROVENANCE_sm.yaml"]
# Only when --starship_sm_dir is a separate dir (e.g. full_v072_antismash/): its
# own summary/provenance files, renamed so they do not overwrite the overlap
# dir's files of the same name. (source name, published name).
STARSHIP_SM_RENAMED = [("summary.tsv", "sm_summary.tsv"), ("PROVENANCE.yaml", "PROVENANCE_sm.yaml"),
                       ("PROVENANCE_analysis.yaml", "PROVENANCE_sm_analysis.yaml")]
# sm summary metrics shown in the SM key-metrics table when present (both schemas).
STARSHIP_SM_KEY_METRICS = [
    "antismash_version", "strains_in_run", "strains_with_antismash", "sm_clusters",
    "sm_clusters_overlap_ge1gene", "sm_clusters_overlap_ge50pct",
    "sm_clusters_overlapping_starship_ge1gene", "sm_clusters_overlapping_starship_ge50pct",
    "sm_clusters_with_captain", "sm_clusters_with_captain_gene", "sm_cluster_genes", "sm_core_genes",
    "elements_with_sm_cluster_ge1gene", "elements_with_sm_cluster_ge50pct",
    "hypergeom_sm_cluster_genes_in_starships_k", "hypergeom_sm_cluster_genes_in_starships_expected",
    "hypergeom_sm_cluster_genes_in_starships_p_greater",
    "sm_core_families", "sm_core_families_noncore", "sm_core_families_noncore_in_island",
    "hypergeom_core_families_in_islands_noncore_expected",
    "hypergeom_core_families_in_islands_noncore_p_greater",
    "hypergeom_cluster_gene_families_in_islands_noncore_k",
    "hypergeom_cluster_gene_families_in_islands_noncore_K",
    "hypergeom_cluster_gene_families_in_islands_noncore_expected",
    "hypergeom_cluster_gene_families_in_islands_noncore_p_greater",
    "captain_genes_in_sm_clusters", "expected_captain_genes_in_sm_clusters", "p_hypergeom_less"]
STARSHIP_SETS = [("starship_all", "Starship families, all elements"),
                 ("starship_high", "Starship families, high-confidence elements"),
                 ("captain", "Captain-gene (DUF3435) families")]
STARSHIP_ELEMENT_COLS = ["element_id", "label", "confidence", "contig", "start", "end", "length_bp",
                         "n_genes", "n_families", "n_noncore_families", "n_families_in_island",
                         "frac_families_in_island", "n_islands_hit_family_level",
                         "n_island_instances_overlapping", "frac_genes_inside_island_instances", "bins"]
STARSHIP_COORD_METRICS = ["n_island_instances_in_annotated_strains",
                          "n_islands_overlapping_element_in_annotated_strain",
                          "n_high_elements_hit_family_level", "n_high_elements_overlapped_by_instance",
                          "island_gene_frac_in_starship_sweep_metric",
                          "noncore_gene_frac_in_starship_sweep_metric"]
STARSHIP_PAGE_CSS = """<style>
.sso table { border-collapse: collapse; margin: 0.5rem 0 1.2rem; font-size: 0.9rem; }
.sso th, .sso td { border: 1px solid #ccc; padding: 0.25rem 0.5rem; text-align: left; vertical-align: top; }
.sso th { background: #eee; }
.sso table.sort th { cursor: pointer; }
.sso td.n { text-align: right; font-variant-numeric: tabular-nums; }
.sso .limits { border: 1px solid #c9a227; background: #fff8e1; border-radius: 8px; padding: 0.6rem 1rem; }
.sso .na { color: #8a5a00; }
.sso .wrap { overflow-x: auto; max-width: 100%; }
@media (prefers-color-scheme: dark) {
  .sso th, .sso td { border-color: #444; }
  .sso th { background: #2a2d33; }
  .sso .limits { background: #2b2614; border-color: #8a6d12; }
  .sso .na { color: #e0b252; }
}
</style>"""
_SORT_JS = ("<script>document.querySelectorAll('.sso table.sort').forEach(t=>t.querySelectorAll('th')"
            ".forEach((th,i)=>th.onclick=()=>{const b=t.tBodies[0],d=th.dataset.d=th.dataset.d==='1'?-1:1;"
            "[...b.rows].sort((x,y)=>{const a=x.cells[i].textContent,c=y.cells[i].textContent,"
            "n=parseFloat(a)-parseFloat(c);return d*(isNaN(n)?a.localeCompare(c):n)})"
            ".forEach(r=>b.appendChild(r))}));</script>")

# Files of a locus-summary fig/ dir (studies/*/analysis/locus_summary/make_locus_figures.py)
# that are published. The per-locus PNGs (about 0.5 MB each) are not; the PDF has them all.
LOCUS_SUMMARY_FILES = ["locus_overview.png", "starship_overlay.png", "locus_strips_selected.png",
                       "locus_strips.pdf", "idiomorph_candidates.tsv"]
LOCUS_SUMMARY_CAPTIONS = {
    "locus_overview.png": "All loci: carrier fraction against locus type.",
    "starship_overlay.png": "Af293 loci with Starship spans (Af293 only).",
    "locus_strips_selected.png": "Selected loci: one row per strain, Mash-clade order.",
}

DEFAULT_CLINKER_PUBLISH_TOP = 25
CLINKER_KEY_RE = re.compile(r"L(\d+)\.html")


def stage_clinker(pangenome_dir: Path, run_docs: Path, publish_top: int) -> dict | None:
    """Copy up to `publish_top` clinker/<L###>.html pages (in ascending L###
    order) into run_docs/clinker/. `publish_top < 0` means all; 0 means none.
    Non-locus files (anything not matching L<digits>.html) are ignored, not
    counted or copied. Returns None if the pipeline produced no clinker/ dir
    at all (N4), else {published, total, full_dir, keys} for the caller to
    fold into run.json and, if published < total, inject into
    island_synteny.html."""
    src = pangenome_dir / "clinker"
    if not src.is_dir():
        return None
    keyed = []
    for p in sorted(src.iterdir()):
        m = CLINKER_KEY_RE.fullmatch(p.name)
        if m:
            keyed.append((int(m.group(1)), p))
    keyed.sort(key=lambda t: t[0])
    total = len(keyed)
    chosen = keyed if publish_top < 0 else keyed[:publish_top]
    if chosen:
        dest = run_docs / "clinker"
        dest.mkdir(exist_ok=True)
        for _, p in chosen:
            shutil.copyfile(p, dest / p.name)
    return {"published": len(chosen), "total": total, "full_dir": str(src.resolve()),
            "keys": [p.stem for _, p in chosen]}


def inject_clinker_published(html: str, keys: list[str], full_dir: str) -> str:
    """island_synteny.html's clinker panel (lib/island_locus_template.py's
    N3) reads window.CLINKER_PUBLISHED/window.CLINKER_FULL_DIR to show a
    note instead of an iframe for a locus that exists in clinker/ but wasn't
    published. A no-op if already present (idempotent on re-sync) -- though
    in practice island_synteny.html is always freshly copied from the
    pipeline's own (never-injected) copy on every sync, so this only
    guards against a future caller reusing an already-staged page as the
    source."""
    keys_json = json.dumps(keys).replace("</", "<\\/")
    dir_json = json.dumps(full_dir).replace("</", "<\\/")
    script = ("<script>window.CLINKER_PUBLISHED=" + keys_json +
              ";window.CLINKER_FULL_DIR=" + dir_json + ";</script>")
    if script in html:
        return html
    return html.replace("</head>", script + "</head>", 1)


def find_pangenome_dir(study_dir: Path, run: str) -> Path:
    run_dir = study_dir / "results" / run
    matches = sorted(p.parent.parent for p in run_dir.glob("*/pangenome/report/report.md"))
    if len(matches) != 1:
        raise SystemExit(
            f"ERROR: expected exactly one {run_dir}/*/pangenome/report/report.md, "
            f"found {len(matches)}: {[str(m) for m in matches]}. "
            "Pass --pangenome-dir to choose one."
        )
    return matches[0]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_copy(src: Path, dest: Path) -> None:
    """Deterministic gzip (mtime=0, no embedded file name), so an unchanged
    table gives a byte-identical .gz on every sync."""
    with open(src, "rb") as fin, open(dest, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fout:
            shutil.copyfileobj(fin, fout)


def count_data_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    with open(path) as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def stage_locus_summary(fig_dir: Path, run_docs: Path) -> bool:
    """Copy the published subset of a locus-summary fig/ dir into
    run_docs/locus_summary/ and write run_docs/locus_summary.html. Returns
    False when fig_dir lacks locus_table.tsv."""
    import csv
    import html
    from collections import Counter
    table = fig_dir / "locus_table.tsv"
    if not table.is_file():
        return False
    dest = run_docs / "locus_summary"
    dest.mkdir(exist_ok=True)
    present = [f for f in LOCUS_SUMMARY_FILES if (fig_dir / f).is_file()]
    for f in present:
        shutil.copy2(fig_dir / f, dest / f)
    with open(table, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = ["key", "label", "type", "size", "n_strains", "content_n", "content_frac", "n_states",
            "tier", "starship", "pair"]
    cols = [c for c in cols if rows and c in rows[0]]
    types = Counter(r["type"] for r in rows)
    body = ['<p><a href="report.html">back to run report</a></p>',
            f"<p>{len(rows)} loci, one focal species. Exploratory: the <code>unresolved</code> "
            "type is the majority, and <code>swap_candidate</code> is verified on one locus only. "
            "Starship overlap uses Af293 coordinates. Strip colours: gene in place; DNA at site, "
            "no gene; DNA checked, absent; not DNA-checked; elsewhere or contig break.</p>",
            "<p>Types: " + ", ".join(f"{html.escape(t)} {n}" for t, n in types.most_common()) + "</p>"]
    for f in present:
        if f in LOCUS_SUMMARY_CAPTIONS:
            body.append(f'<h3>{html.escape(LOCUS_SUMMARY_CAPTIONS[f])}</h3>'
                        f'<p><img src="locus_summary/{f}" style="max-width:100%"></p>')
    if "locus_strips.pdf" in present:
        body.append('<p><a href="locus_summary/locus_strips.pdf">Strips for all loci (PDF)</a></p>')
    if "idiomorph_candidates.tsv" in present:
        body.append('<p><a href="locus_summary/idiomorph_candidates.tsv">Swap candidates (TSV)</a></p>')
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    trs = "".join("<tr>" + "".join(f"<td>{html.escape(r.get(c, ''))}</td>" for c in cols) + "</tr>"
                  for r in rows)
    body.append(f'<table id="lt"><thead><tr>{head}</tr></thead><tbody>{trs}</tbody></table>')
    body.append("<script>document.querySelectorAll('#lt th').forEach((th,i)=>th.onclick=()=>{"
                "const b=document.querySelector('#lt tbody'),d=th.dataset.d=th.dataset.d==='1'?-1:1;"
                "[...b.rows].sort((x,y)=>{const a=x.cells[i].textContent,c=y.cells[i].textContent,"
                "n=parseFloat(a)-parseFloat(c);return d*(isNaN(n)?a.localeCompare(c):n)})"
                ".forEach(r=>b.appendChild(r))});</script>")
    (run_docs / "locus_summary.html").write_text(_page(
        f"{run_docs.parent.name} / {run_docs.name} / locus summary", "\n".join(body)))
    return True


def _read_tsv(path: Path) -> list[dict] | None:
    import csv
    if not path.is_file():
        return None
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


_NUM_RE = re.compile(r"-?[\d.]+(e[-+]?\d+)?|inf|NA")


def _html_table(cols: list[str], rows: list[list[str]], sortable: bool = False,
                headers: list[str] | None = None) -> str:
    import html
    cls = ' class="sort"' if sortable else ""
    head = "".join(f"<th>{html.escape(h)}</th>" for h in (headers or cols))

    def td(v) -> str:
        v = "" if v is None else str(v)
        if v and _NUM_RE.fullmatch(v):
            return f'<td class="n">{html.escape(v)}</td>'
        return f"<td>{html.escape(v)}</td>"
    trs = "".join("<tr>" + "".join(td(v) for v in r) + "</tr>" for r in rows)
    return f'<div class="wrap"><table{cls}><thead><tr>{head}</tr></thead><tbody>{trs}</tbody></table></div>'


def _has_cols(rows: list[dict] | None, path: Path, cols: tuple[str, ...]) -> bool:
    """True when the table's header has every column in `cols`; otherwise
    print a warning naming the columns found (an unexpected schema must not
    crash a sync). An empty table (header only) has no row to check: True."""
    if not rows:
        return True
    have = list(rows[0].keys())
    if all(c in have for c in cols):
        return True
    print(f"WARNING: {path}: expected columns {list(cols)}, found {have}", file=sys.stderr)
    return False


def _find_named_figure(overlap_dir: Path, suffix: str) -> Path | None:
    """make_figures.py writes <figures>/<label>_<suffix> next to (or above) the
    overlap dir, with <label> = the overlap dir's path below that parent joined
    by '_' (full_v072 -> full_v072_starship_islands.png;
    full_v072/with_unplaced_af293_coords -> full_v072_with_unplaced_af293_coords_sm.png).
    Exact names only, so full_v072 never picks up the with_unplaced figure.
    A file named <suffix> inside the dir itself wins."""
    if (overlap_dir / suffix).is_file():
        return overlap_dir / suffix
    for anc in list(overlap_dir.parents)[:3]:
        cand = anc / "figures" / ("_".join(overlap_dir.relative_to(anc).parts) + "_" + suffix)
        if cand.is_file():
            return cand
    return None


def stage_starship_overlap(overlap_dir: Path, run_docs: Path, sm_dir: Path | None = None) -> bool:
    """Copy a bin/starship_island_overlap.py output dir's small tables and PNGs
    into run_docs/starship_overlap/ and write run_docs/starship_overlap.html.
    SM tables (STARSHIP_SM_FILES, figures) come from `sm_dir` (default:
    overlap_dir). A separate sm_dir also contributes its summary.tsv
    (metric/value/source) and PROVENANCE*.yaml under new names
    (STARSHIP_SM_RENAMED). sm_starship_overlap.tsv may be a metric table
    (metric/value/definition) or a per-element table (element_id, ...).
    Every input except summary.tsv is optional; a missing one, or one with
    unexpected columns, gives a warning and a 'not available' line, never an
    exception. Returns False (no page) when summary.tsv is missing or has no
    metric/value columns."""
    import html
    summary_rows = _read_tsv(overlap_dir / "summary.tsv")
    if summary_rows is None:
        print(f"WARNING: {overlap_dir}/summary.tsv not found -- no starship_overlap page",
              file=sys.stderr)
        return False
    if not _has_cols(summary_rows, overlap_dir / "summary.tsv", ("metric", "value")):
        print(f"WARNING: {overlap_dir}/summary.tsv has no metric/value columns -- no starship_overlap page",
              file=sys.stderr)
        return False
    sm_dir = sm_dir or overlap_dir
    separate_sm = sm_dir.resolve() != overlap_dir.resolve()
    dest = run_docs / "starship_overlap"
    dest.mkdir(exist_ok=True)
    copied: set[str] = set()
    pairs = [(overlap_dir / f, f) for f in STARSHIP_OVERLAP_FILES] + [(sm_dir / f, f) for f in STARSHIP_SM_FILES]
    if separate_sm:
        pairs += [(sm_dir / s, d) for s, d in STARSHIP_SM_RENAMED]
    for src, name in pairs:
        if src.is_file() and name not in copied:
            shutil.copyfile(src, dest / name)
            copied.add(name)
    inputs_prov = next((a / "inputs" / "PROVENANCE.yaml" for a in [overlap_dir, *list(overlap_dir.parents)[:2]]
                        if (a / "inputs" / "PROVENANCE.yaml").is_file()), None)
    if inputs_prov is not None:
        shutil.copyfile(inputs_prov, dest / "inputs_PROVENANCE.yaml")
        copied.add("inputs_PROVENANCE.yaml")
    isl_fig = _find_named_figure(overlap_dir, "starship_islands.png")
    if isl_fig is not None:
        shutil.copyfile(isl_fig, dest / "starship_islands.png")
    sm_figs: list[str] = []
    named_sm = _find_named_figure(sm_dir, "sm.png")
    if named_sm is not None:
        shutil.copyfile(named_sm, dest / "sm.png")
        sm_figs.append("sm.png")
    if separate_sm:
        for p in sorted([*sm_dir.glob("*.png"), *sm_dir.glob("figures/*.png")]):
            if p.name not in ("sm.png", "starship_islands.png") and p.name not in sm_figs:
                shutil.copyfile(p, dest / p.name)
                sm_figs.append(p.name)

    S = {r["metric"]: r["value"] for r in summary_rows}
    D = {r["metric"]: r.get("definition") or "" for r in summary_rows}
    enr = _read_tsv(overlap_dir / "enrichment.tsv") or []
    if enr and not _has_cols(enr, overlap_dir / "enrichment.tsv", ("annotation", "universe")):
        enr = []
    E = {(r["annotation"], r["universe"]): r for r in enr}
    sets = [(k, lab) for k, lab in STARSHIP_SETS if f"n_{k}_families" in S]
    run_label = summary_rows[0].get("run", "") if summary_rows else ""

    def na(text: str) -> str:
        return f'<p class="na">{text}</p>'

    def link(f: str) -> str:
        return f'<a href="starship_overlap/{f}">{html.escape(f)}</a>'

    def headline(keys: list[tuple[str, str]]) -> str:
        cols = ["set", "families", "in >= 1 island", "in >= 1 multi-strain island",
                "islands with >= 1 family of the set", "non-core: set families in islands",
                "non-core: all families in islands (background)", "fold",
                "hypergeometric p (greater)"]
        rows = []
        for k, lab in keys:
            e = E.get((k, "noncore"), {})
            rows.append([lab, S.get(f"n_{k}_families", ""), S.get(f"n_{k}_families_in_island", ""),
                         S.get(f"n_{k}_families_in_multi_strain_island", ""),
                         f"{S.get(f'n_islands_with_{k}_family', 'NA')} of {S.get('n_islands', 'NA')}",
                         f"{e.get('k_both', 'NA')} of {e.get('K_annotated', 'NA')} "
                         f"({e.get('frac_annotated_in_island', 'NA')})",
                         f"{e.get('n_island', 'NA')} of {e.get('N_universe', 'NA')} "
                         f"({e.get('frac_universe_in_island', 'NA')})",
                         e.get("fold", "NA"), e.get("p_hypergeom_greater", "NA")])
        return _html_table(cols, rows)

    def enrichment_rows(keys: list[str]) -> str:
        rows = [r for r in enr if r["annotation"] in keys]
        if not rows:
            return na("enrichment.tsv: not available.")
        cols = ["annotation", "universe", "N_universe", "K_annotated", "n_island", "k_both", "expected_k",
                "fold", "odds_ratio", "p_hypergeom_greater"]
        return ("<details><summary>All universes (enrichment.tsv)</summary>"
                + _html_table(cols, [[r.get(c, "") for c in cols] for r in rows]) + "</details>")

    body = [STARSHIP_PAGE_CSS, '<div class="sso">',
            '<p><a href="report.html">back to run report</a></p>',
            "<h1>Starships, islands and secondary-metabolite clusters</h1>",
            f"<p>Source: <code>bin/starship_island_overlap.py</code> output for run "
            f"<code>{html.escape(run_label)}</code>. A family is 'in an island' when it is a member of "
            "&gt;= 1 row of significant_islands.tsv. Background: all non-core families of the run "
            "(enrichment.tsv, universe <code>noncore</code>). Exploratory analysis; read the "
            "'Data and limits' box.</p>"]
    run_cols = [("n_families", "families"), ("n_islands", "islands"),
                ("n_islands_multi_strain", "islands with >= 2 strains"),
                ("n_island_families", "families in >= 1 island"),
                ("n_island_families_multi_strain", "families in >= 1 multi-strain island")]
    present_cols = [(c, h) for c, h in run_cols if c in S]
    if present_cols:
        body.append(_html_table([c for c, _ in present_cols], [[S[c] for c, _ in present_cols]],
                                headers=[h for _, h in present_cols]))

    body.append("<h2>Starship families in islands</h2>")
    ship_sets = [s for s in sets if s[0].startswith("starship")]
    if ship_sets:
        body.append(f"<p>Starship family: a family with &gt;= 1 gene of "
                    f"{html.escape(S.get('annotated_strains', 'NA'))} whose midpoint lies inside one of "
                    f"{html.escape(S.get('n_elements', 'NA'))} element(s) "
                    f"({html.escape(S.get('n_elements_high', 'NA'))} high confidence).</p>")
        body.append(headline(ship_sets))
        body.append(enrichment_rows([k for k, _ in ship_sets]))
        coord = [[m, S[m], D.get(m, "")] for m in STARSHIP_COORD_METRICS if m in S]
        if coord:
            body.append("<h3>Coordinate checks in the annotated strain</h3>")
            body.append(_html_table(["metric", "value", "definition"], coord))
        perm = _read_tsv(overlap_dir / "coord_permutation.tsv")
        if perm:
            cols = list(perm[0].keys())
            body.append("<h3>Coordinate permutation</h3>"
                        + _html_table(cols, [[r[c] for c in cols] for r in perm]))
        per = _read_tsv(overlap_dir / "per_starship.tsv")
        body.append("<h3>Per element</h3>")
        if per:
            cols = [c for c in STARSHIP_ELEMENT_COLS if c in per[0]]
            body.append(_html_table(cols, [[r.get(c, "") for c in cols] for r in per], sortable=True))
        else:
            body.append(na("per_starship.tsv: not available."))
    else:
        body.append(na("Not available: this run has no Starship element coordinates "
                       "(summary.tsv has no starship_* metrics). The captain-gene families below are "
                       "the only Starship proxy."))

    body.append("<h2>Captain-gene families (non-core) in islands</h2>")
    if any(k == "captain" for k, _ in sets):
        body.append(f"<p>Captain gene: {html.escape(D.get('captain_hits') or 'DUF3435 hit')}; "
                    f"{html.escape(S.get('captain_hits', 'NA'))} hits, "
                    f"{html.escape(S.get('captain_hits_unmapped', 'NA'))} not mapped to a family "
                    "of this run.</p>")
        body.append(headline([s for s in sets if s[0] == "captain"]))
        body.append(enrichment_rows(["captain"]))
    else:
        body.append(na("Not available: summary.tsv has no captain metrics."))
    cb = _read_tsv(overlap_dir / "class_breakdown.tsv")
    if cb:
        cols = list(cb[0].keys())
        body.append("<details><summary>Frequency bin and gain/loss direction (class_breakdown.tsv)"
                    "</summary>" + _html_table(cols, [[r[c] for c in cols] for r in cb], sortable=True)
                    + "</details>")
    if isl_fig is not None:
        cap = "Families in islands by annotation set (per element when element coordinates exist)."
        body.append(f'<h3>{cap}</h3><p><img src="starship_overlap/starship_islands.png" '
                    f'alt="{cap}" style="max-width:100%"></p>')

    body.append("<h2>Secondary-metabolite (SM) clusters</h2>")
    # SM metrics, two schemas: (a) sm_starship_overlap.tsv with metric/value/definition
    # (Coccidioides), (b) a separate SM dir whose summary.tsv has metric/value/source and
    # whose sm_starship_overlap.tsv is a per-element table (Af293 antiSMASH).
    M: dict[str, tuple[str, str]] = {}
    per_element: list[dict] | None = None
    smo_path = sm_dir / "sm_starship_overlap.tsv"
    smo = _read_tsv(smo_path) or []
    smo_note = ""
    if smo:
        cols0 = list(smo[0].keys())
        if "metric" in cols0:
            M.update({r["metric"]: (r.get("value") or "", r.get("definition") or "") for r in smo})
        elif "element_id" in cols0:
            per_element = smo
        else:
            print(f"WARNING: {smo_path}: unexpected columns {cols0}", file=sys.stderr)
            smo_note = ("Not available: sm_starship_overlap.tsv has unexpected columns "
                        + ", ".join(cols0) + ".")
    if separate_sm:
        sm_sum = _read_tsv(sm_dir / "summary.tsv") or []
        if sm_sum and _has_cols(sm_sum, sm_dir / "summary.tsv", ("metric", "value")):
            M.update({r["metric"]: (r.get("value") or "", r.get("source") or r.get("definition") or "")
                      for r in sm_sum})
    prov_sm = (dest / "PROVENANCE_sm.yaml").read_text() if "PROVENANCE_sm.yaml" in copied else ""
    has_family_map = (sm_dir / "sm_family_class_direction.tsv").is_file()
    sm = _read_tsv(sm_dir / "sm_by_type.tsv")
    if sm is None:
        sm_kind = "not available (no sm_by_type.tsv)"
        body.append(na("Not available: no sm_by_type.tsv for this run."))
    else:
        cols = list(sm[0].keys()) if sm else []
        if "n_clusters" in cols:
            if M.get("antismash_version", ("", ""))[0]:
                ver = "antiSMASH " + M["antismash_version"][0]
            else:
                m_ver = re.search(r"antiSMASH \d[\d.]*\d", prov_sm + " " + " ".join(d for _, d in M.values()))
                ver = m_ver.group(0) if m_ver else "antiSMASH"
            sm_kind = ver + " regions"
        elif "af293_genes" in cols:
            sm_kind = ("proxy: PKS KS (PF00109) and NRPS C (PF00668) backbone genes only, "
                       "not antiSMASH clusters")
        else:
            sm_kind = "unknown source (see the SM provenance file)"
        body.append(f"<p>Source: {html.escape(sm_kind)}.</p>")
        key = [[k, M[k][0][:300], M[k][1]] for k in STARSHIP_SM_KEY_METRICS if k in M]
        if key:
            body.append("<h3>SM key metrics</h3>" + _html_table(["metric", "value", "definition / source"], key))
        if not any("starship" in c.lower() or "overlap" in c.lower() for c in cols):
            ge1 = M.get("sm_clusters_overlapping_starship_ge1gene") or M.get("sm_clusters_overlap_ge1gene")
            ge50 = M.get("sm_clusters_overlapping_starship_ge50pct") or M.get("sm_clusters_overlap_ge50pct")
            if ge1 or ge50:
                body.append(na("Overlap with a Starship: &gt;= 1 gene rule: "
                               f"{html.escape(' '.join(ge1 or ('NA', '')).strip())}; &gt;= 50% rule: "
                               f"{html.escape(' '.join(ge50 or ('NA', '')).strip())}."))
            else:
                body.append(na("Overlap with a Starship: not available for this SM table."))
        body.append("<h3>SM clusters by type</h3>")
        if sm:
            body.append(_html_table(cols, [[r.get(c, "") for c in cols] for r in sm], sortable=True))
        else:
            body.append(na("sm_by_type.tsv has no rows."))
    if per_element:
        cols = list(per_element[0].keys())
        body.append("<h3>Table S7 elements and SM clusters</h3>"
                    + _html_table(cols, [[r.get(c, "") for c in cols] for r in per_element], sortable=True))
    if smo_note:
        body.append(na(html.escape(smo_note)))
    perm_sm = _read_tsv(sm_dir / "sm_permutation.tsv")
    if perm_sm:
        cols = list(perm_sm[0].keys())
        body.append("<h3>SM cluster position permutation</h3>"
                    + _html_table(cols, [[r.get(c, "") for c in cols] for r in perm_sm]))
    fam = _read_tsv(sm_dir / "sm_family_class_direction.tsv")
    if fam:
        cols = list(fam[0].keys())
        body.append("<details><summary>SM core-gene families by class, direction and island "
                    "(sm_family_class_direction.tsv)</summary>"
                    + _html_table(cols, [[r.get(c, "") for c in cols] for r in fam], sortable=True)
                    + "</details>")
    if M:
        body.append("<details><summary>All SM summary metrics</summary>"
                    + _html_table(["metric", "value", "definition / source"],
                                  [[k, v[:300], d] for k, (v, d) in M.items()]) + "</details>")
    for f in sm_figs:
        body.append(f'<p><img src="starship_overlap/{html.escape(f)}" alt="SM clusters figure" '
                    'style="max-width:100%"></p>')

    lim = ["<h2>Data and limits</h2>", '<div class="limits"><ul>']
    if ship_sets:
        lim.append(f"<li>Starship element coordinates exist only for "
                   f"{html.escape(S.get('annotated_strains', 'NA'))} "
                   f"({html.escape(S.get('n_elements', 'NA'))} elements). Family-level results transfer "
                   "these calls to all strains through family membership; coordinate-level results apply "
                   "to that strain only.</li>")
    else:
        lim.append("<li>No genome of this run has Starship element coordinates. Captain-gene (DUF3435) "
                   "families are a proxy: a captain gene marks a possible Starship, not its span.</li>")
    lim.append(f"<li>SM annotation: {html.escape(sm_kind)}.</li>")
    m_src = re.search(r"^source_annotation:\s*(.+)$", prov_sm, re.M)
    if m_src:
        lim.append("<li>SM annotation source (SM provenance): "
                   f"<code>{html.escape(m_src.group(1).strip())}</code>.</li>")
    if has_family_map:
        mapped = next(((k, v) for k, (v, _) in M.items() if k.endswith("_cds_mapped_to_member")), None)
        total = next((v for k, (v, _) in M.items() if k.endswith("_cds_in_antismash")), None)
        extra = (f" ({html.escape(mapped[0])}: {html.escape(mapped[1])} of {html.escape(total or 'NA')} CDS)"
                 if mapped else "")
        lim.append("<li>SM cluster genes are mapped to pangenome families of this run"
                   f"{extra}; see sm_family_class_direction.tsv.</li>")
    elif m_src:
        lim.append("<li>This page does not map SM cluster genes to the pangenome run's protein IDs "
                   "or to islands.</li>")
    if "strains_with_antismash" in M:
        lim.append(f"<li>Strains with SM annotation: {html.escape(M['strains_with_antismash'][0])} of "
                   f"{html.escape(M.get('strains_in_run', ('NA', ''))[0])}.</li>")
    prov = [f for f in ("PROVENANCE.yaml", "PROVENANCE_sm.yaml", "PROVENANCE_sm_analysis.yaml",
                        "inputs_PROVENANCE.yaml") if f in copied]
    lim.append("<li>Provenance: " + (", ".join(link(f) for f in prov) if prov else "not available")
               + ".</li>")
    lim.append("</ul></div>")
    body.extend(lim)
    published = STARSHIP_OVERLAP_FILES + STARSHIP_SM_FILES + [d for _, d in STARSHIP_SM_RENAMED]
    tables = [f for f in dict.fromkeys(published) if f in copied and not f.endswith(".yaml")]
    body.append("<h2>Tables</h2><ul>" + "".join(f"<li>{link(f)}</li>" for f in tables) + "</ul>")
    body.append("</div>")
    body.append(_SORT_JS)
    (run_docs / "starship_overlap.html").write_text(_page(
        f"{run_docs.parent.name} / {run_docs.name} / Starships, islands and SM clusters", "\n".join(body)))
    return True


def stage_run(pangenome_dir: Path, run_docs: Path, domain: str, set_name: str, run: str,
              source_label: str, today: str, run_record: Path | None = None,
              clinker_publish_top: int = DEFAULT_CLINKER_PUBLISH_TOP,
              locus_summary_dir: Path | None = None,
              starship_overlap_dir: Path | None = None,
              starship_sm_dir: Path | None = None) -> dict:
    run_docs.mkdir(parents=True, exist_ok=True)
    for name in RELEASE_ONLY:
        target = run_docs / name
        if target.is_dir():
            shutil.rmtree(target)
        elif target.exists():
            target.unlink()

    for src_rel, dest_name in ASSET_DIRS:
        src = pangenome_dir / src_rel
        if src.is_dir():
            shutil.copytree(src, run_docs / dest_name)

    clinker_info = stage_clinker(pangenome_dir, run_docs, clinker_publish_top)

    downloads: list[tuple[str, str]] = []
    archive = run_docs / "archive"
    for src_rel, archive_name in ARCHIVE_TABLES:
        src = pangenome_dir / src_rel
        # 0 bytes = the pipeline's "not computed" convention: no dead download.
        if src.is_file() and src.stat().st_size > 0:
            archive.mkdir(exist_ok=True)
            gzip_copy(src, archive / archive_name)
            downloads.append((archive_name.removesuffix(".tsv.gz"), f"archive/{archive_name}"))

    extra_pages: list[tuple[str, str]] = []
    synteny = pangenome_dir / "island_synteny.html"
    if synteny.is_file():
        html = synteny.read_text()
        if clinker_info is not None and clinker_info["published"] < clinker_info["total"]:
            html = inject_clinker_published(html, clinker_info["keys"], clinker_info["full_dir"])
        (run_docs / "island_synteny.html").write_text(html)
        extra_pages.append(("Island synteny viewer", "island_synteny.html"))
    if locus_summary_dir is not None and stage_locus_summary(locus_summary_dir, run_docs):
        extra_pages.append(("Locus summary", "locus_summary.html"))
    if starship_overlap_dir is not None and stage_starship_overlap(starship_overlap_dir, run_docs,
                                                                   starship_sm_dir):
        extra_pages.append(("Starships, islands and SM clusters", "starship_overlap.html"))
    aq_md = pangenome_dir / "assembly_quality_report.md"
    if aq_md.is_file():
        (run_docs / "assembly_quality.html").write_text(_page(
            f"{set_name} / {run} / assembly quality QC",
            '<p><a href="report.html">back to run report</a></p>'
            + render_markdown(aq_md.read_text()),
        ))
        extra_pages.append(("Assembly quality QC", "assembly_quality.html"))

    report_md = pangenome_dir / "report" / "report.md"
    (run_docs / "report.html").write_text(render_run_report(
        report_md.read_text(), domain, set_name, run, downloads, extra_pages,
    ))
    (run_docs / "index.html").write_text(render_report_redirect())

    meta = {
        "run": run,
        "published": today,
        "pipeline": "pangenome.nf",
        "source_dir": source_label,
        "report_md_sha256": sha256_of(report_md),
        "n_families": count_data_rows(pangenome_dir / "frequency_table.tsv"),
        "n_strains": count_data_rows(pangenome_dir / "report_tables" / "per_strain_summary.tsv"),
    }
    if run_record is not None and run_record.is_file():
        rec = json.loads(run_record.read_text())
        meta["pipeline_repo"] = rec.get("pipeline")
        meta["pipeline_commit"] = rec.get("pipeline_commit")
    # report/rerender.json: written when report steps were re-run offline at a
    # newer pipeline commit than the run itself (NovInvenio #212 rollout).
    rerender = pangenome_dir / "report" / "rerender.json"
    if rerender.is_file():
        meta["report_rerender"] = json.loads(rerender.read_text())
    if clinker_info is not None:
        meta["clinker_published"] = clinker_info["published"]
        meta["clinker_total"] = clinker_info["total"]
        meta["clinker_full_dir"] = clinker_info["full_dir"]
    (run_docs / "run.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="<domain>/<folder>, e.g. fungi/coccidioides_pangenome")
    ap.add_argument("--run", required=True, help="run name under studies/<domain>/<folder>/results/")
    ap.add_argument("--current", action="store_true", help="mark this run as the study's current run")
    ap.add_argument("--pangenome-dir", type=Path, default=None,
                    help="override the auto-detected <run>/*/pangenome/ dir")
    ap.add_argument("--clinker_publish_top", type=int, default=DEFAULT_CLINKER_PUBLISH_TOP,
                    help="publish only the top N loci' clinker/<L###>.html pages "
                         f"(default: {DEFAULT_CLINKER_PUBLISH_TOP}; 0 = publish none; "
                         "-1 = publish all)")
    ap.add_argument("--locus_summary_dir", type=Path, default=None,
                    help="fig/ dir written by make_locus_figures.py; adds a Locus summary page "
                         "(release asset only). Not auto-detected, so a stale dir is never published.")
    ap.add_argument("--starship_overlap_dir", type=Path, default=None,
                    help="bin/starship_island_overlap.py --out-dir of this run; adds a 'Starships, "
                         "islands and SM clusters' page (release asset only). Not auto-detected.")
    ap.add_argument("--starship_sm_dir", type=Path, default=None,
                    help="dir with sm_by_type.tsv / sm_starship_overlap.tsv / sm_clusters.tsv.gz / "
                         "PROVENANCE_sm.yaml / *.png for that page (default: --starship_overlap_dir)")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    domain, _, folder = args.study.partition("/")
    if not domain or not folder or "/" in folder:
        raise SystemExit(f"ERROR: --study must be <domain>/<folder> (got {args.study!r})")
    try:
        validate_slug(args.run, "--run")
    except ValueError as e:
        raise SystemExit(f"ERROR: {e}")
    study_dir = args.repo_root / "studies" / domain / folder
    if not (study_dir / "species.csv").exists():
        raise SystemExit(f"ERROR: {study_dir}/species.csv not found -- not a defined study")
    target = load_publish_target(study_dir)
    if target is None:
        print(f"== {study_dir}/publish.yaml not found -- not published. Add one "
              "(study: <name>) to publish this folder. ==", file=sys.stderr)
        return 0
    set_name = target.study

    pangenome_dir = args.pangenome_dir or find_pangenome_dir(study_dir, args.run)
    if not (pangenome_dir / "report" / "report.md").is_file():
        raise SystemExit(f"ERROR: {pangenome_dir}/report/report.md not found")

    try:
        source_label = str(pangenome_dir.resolve().relative_to(args.repo_root.resolve()))
    except ValueError:
        source_label = str(pangenome_dir.resolve())

    study_docs = args.repo_root / "docs" / domain / set_name
    today = datetime.datetime.now(datetime.UTC).date().isoformat()
    meta = stage_run(pangenome_dir, study_docs / args.run, domain, set_name, args.run,
                     source_label, today,
                     run_record=study_dir / ".nf_launch" / args.run / "ni_run.json",
                     clinker_publish_top=args.clinker_publish_top,
                     locus_summary_dir=args.locus_summary_dir,
                     starship_overlap_dir=args.starship_overlap_dir,
                     starship_sm_dir=args.starship_sm_dir)
    if args.current:
        set_current(study_docs, args.run)
    runs = rebuild_study_pages(study_docs, domain, set_name)
    print(f"staged {args.study}/{args.run} from {source_label} "
          f"({meta['n_families']} families, {meta['n_strains']} strains); "
          f"study index lists {len(runs)} run(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
