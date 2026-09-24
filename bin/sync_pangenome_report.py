#!/usr/bin/env python3
"""Stage one pangenome.nf run's report into the docs/ gallery.

pangenome.nf writes report/report.md + figures, not the novelties/core/losses
HTML set that bin/sync_reports.sh handles. This script is the pangenome
counterpart. Layout and data classes: see lib/pangenome_site.py's docstring.

Steps, for --study <domain>/<set> --run <run>:
  1. Find the run's pangenome output dir:
     studies/<domain>/<set>/results/<run>/*/pangenome/ (the `*` is the
     pipeline's --outdir project subdir). Exactly one match is required;
     pass --pangenome-dir to override.
  2. Replace docs/<domain>/<set>/<run>/{figures,figures_pdf,archive}/ and
     island_synteny.html / assembly_quality.html with this run's copies.
     These are gitignored (release asset only).
  3. Render docs/<domain>/<set>/<run>/report.html from report/report.md, and
     write run.json (source path, sha256 of report.md, counts).
  4. Rebuild docs/<domain>/<set>/report.html (the run list) from every
     docs/<domain>/<set>/*/run.json.

It does not publish anything. Run bin/generate_docs.py after it to refresh the
gallery, then bin/publish_report_release.sh <domain>/<set>/<run> to upload the
release asset.

Usage:
  bin/sync_pangenome_report.py --study fungi/coccidioides_pangenome \\
      --run rescue_structural_genus_vs_ureesii
"""
from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from pangenome_site import render_markdown, render_run_index, render_run_report  # noqa: E402
from site_pages import _page, render_report_redirect  # noqa: E402

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
]

# Directories and pages copied or rendered into the run dir. All gitignored.
ASSET_DIRS = [("report/figures", "figures"), ("report/figures_pdf", "figures_pdf")]
RELEASE_ONLY = ["figures", "figures_pdf", "archive", "island_synteny.html", "assembly_quality.html"]


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


def stage_run(pangenome_dir: Path, run_docs: Path, domain: str, set_name: str, run: str,
              source_label: str, today: str) -> dict:
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

    downloads: list[tuple[str, str]] = []
    archive = run_docs / "archive"
    for src_rel, archive_name in ARCHIVE_TABLES:
        src = pangenome_dir / src_rel
        if src.is_file():
            archive.mkdir(exist_ok=True)
            gzip_copy(src, archive / archive_name)
            downloads.append((archive_name.removesuffix(".tsv.gz"), f"archive/{archive_name}"))

    extra_pages: list[tuple[str, str]] = []
    synteny = pangenome_dir / "island_synteny.html"
    if synteny.is_file():
        shutil.copyfile(synteny, run_docs / "island_synteny.html")
        extra_pages.append(("Island synteny viewer", "island_synteny.html"))
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
        "source_dir": source_label,
        "report_md_sha256": sha256_of(report_md),
        "n_families": count_data_rows(pangenome_dir / "frequency_table.tsv"),
        "n_strains": count_data_rows(pangenome_dir / "report_tables" / "per_strain_summary.tsv"),
    }
    (run_docs / "run.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def rebuild_study_index(study_docs: Path, domain: str, set_name: str) -> list[dict]:
    runs = []
    for meta_path in sorted(study_docs.glob("*/run.json")):
        with open(meta_path) as fh:
            runs.append(json.load(fh))
    (study_docs / "report.html").write_text(render_run_index(domain, set_name, runs))
    return runs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="<domain>/<set_name>, e.g. fungi/coccidioides_pangenome")
    ap.add_argument("--run", required=True, help="run name under studies/<domain>/<set>/results/")
    ap.add_argument("--pangenome-dir", type=Path, default=None,
                    help="override the auto-detected <run>/*/pangenome/ dir")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    domain, _, set_name = args.study.partition("/")
    if not domain or not set_name or "/" in set_name or "/" in args.run:
        raise SystemExit(f"ERROR: --study must be <domain>/<set> and --run a single name "
                         f"(got {args.study!r}, {args.run!r})")
    study_dir = args.repo_root / "studies" / domain / set_name
    if not (study_dir / "species.csv").exists():
        raise SystemExit(f"ERROR: {study_dir}/species.csv not found -- not a defined study")

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
                     source_label, today)
    runs = rebuild_study_index(study_docs, domain, set_name)
    print(f"staged {args.study}/{args.run} from {source_label} "
          f"({meta['n_families']} families, {meta['n_strains']} strains); "
          f"study index lists {len(runs)} run(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
