#!/usr/bin/env python3
"""Generate the two-tier docs/ gallery (DESIGN.md Sec 8): top-level domain cards,
each linking to a per-domain page of study cards.

Study/run layout (notes/superpowers/specs/2026-09-24-study-run-site-layout-design.md):
a gallery card is one *study*, which can have several runs at
docs/<domain>/<study>/<run>/ (each with a run.json written by a sync script).

Cards come from two places:
  - docs/<domain>/<study>/*/run.json -- "complete" (at least one published
    run). This script rebuilds that study's run list, index.html redirect,
    and study.json (lib/study_runs.py). The card links to the current run.
  - studies/<domain>/<folder>/ with species.csv AND publish.yaml, grouped
    by publish.yaml's `study`. Folders without publish.yaml are not
    published and get no card. A study with no published run is:
      "pending" -- no config.csv yet (bin/build_study_config.py not run)
      "staged"  -- config.csv/data_dir exist, no run published yet
    Species counts and the hypothesis line come from the folder whose name
    equals the study name when there is one, else the first folder.

Re-run any time a study is added or its status changes (data staged, pipeline run).
Does not touch a study's own report.html/novelties.html/etc -- those are published
by copying nf_NovInvenio's COLLATE_REPORTS output into docs/<domain>/<set>/, same
place this script expects to find them.

The domain registry itself (slug/name/desc for each tab) is curated data, not code
-- see conf/domains.yaml, edited directly rather than here.
"""
import csv
import sys
import tomllib
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from site_pages import render_domain_index, render_top_level  # noqa: E402
from study_runs import load_publish_target, read_current, rebuild_study_pages  # noqa: E402

REPO_ROOT = Path(__file__).parent.parent
STUDIES_ROOT = REPO_ROOT / "studies"
DOCS_ROOT = REPO_ROOT / "docs"
DOMAINS_YAML = REPO_ROOT / "conf" / "domains.yaml"
PIXI_TOML = REPO_ROOT / "pixi.toml"


def load_site_name() -> str:
    """This repo's own display name for the gallery header/title, derived
    from pixi.toml's [workspace] name field (underscores read as spaces)
    rather than hardcoded -- so a repo scaffolded from this one
    (nf_NovInvenio's bin/ni, issue #76) brands its own gallery correctly
    instead of inheriting "NovInvenio Investigations" regardless of what
    it's actually called. Falls back to that same default if pixi.toml is
    missing/unparseable, matching lib/site_pages.py's own default."""
    try:
        with open(PIXI_TOML, "rb") as fh:
            name = tomllib.load(fh)["workspace"]["name"]
        return name.replace("_", " ")
    except (OSError, KeyError, tomllib.TOMLDecodeError):
        return "NovInvenio Investigations"


def load_domains() -> list[tuple[str, str, str]]:
    """Return [(slug, name, desc), ...] from conf/domains.yaml, in file order."""
    with open(DOMAINS_YAML) as fh:
        rows = yaml.safe_load(fh)
    return [(row["slug"], row["name"], row["desc"]) for row in rows]


def species_counts(study_dir: Path) -> tuple[int, int] | tuple[None, None]:
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        return None, None
    n_in = n_out = 0
    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("Group") == "IN":
                n_in += 1
            elif row.get("Group") == "OUT":
                n_out += 1
    return n_in, n_out


def hypothesis_text(study_dir: Path) -> str:
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        return ""
    taxa = []
    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            t = row.get("TaxonGroup")
            if t and t not in taxa:
                taxa.append(t)
    return ", ".join(taxa)


def find_studies(domain_slug: str) -> list[dict]:
    """One dict per study: published studies from docs/, plus unpublished
    studies/ folders that have a publish.yaml."""
    by_study: dict[str, dict] = {}

    domain_docs = DOCS_ROOT / domain_slug
    if domain_docs.exists():
        for study_docs in sorted(domain_docs.iterdir()):
            if not study_docs.is_dir() or not any(study_docs.glob("*/run.json")):
                continue
            runs = rebuild_study_pages(study_docs, domain_slug, study_docs.name)
            current = read_current(study_docs)
            by_study[study_docs.name] = {
                "name": study_docs.name, "slug": study_docs.name,
                "status": "complete", "n_runs": len(runs),
                "href": f"{study_docs.name}/{current}/report.html" if current
                        else f"{study_docs.name}/report.html",
                "hypothesis": "", "n_ingroup": None, "n_outgroup": None,
                "folders": [],
            }

    domain_dir = STUDIES_ROOT / domain_slug
    if domain_dir.exists():
        for folder in sorted(domain_dir.iterdir()):
            if not folder.is_dir() or not (folder / "species.csv").exists():
                continue
            try:
                target = load_publish_target(folder)
            except ValueError as e:
                print(f"WARNING: skipping {folder}: {e}", file=sys.stderr)
                continue
            if target is None:
                continue
            entry = by_study.setdefault(target.study, {
                "name": target.study, "slug": target.study, "status": None,
                "n_runs": 0, "href": None, "hypothesis": "",
                "n_ingroup": None, "n_outgroup": None, "folders": [],
            })
            entry["folders"].append(folder)

    out = []
    for name in sorted(by_study):
        entry = by_study[name]
        folders = entry.pop("folders")
        if folders:
            src = next((f for f in folders if f.name == name), folders[0])
            entry["n_ingroup"], entry["n_outgroup"] = species_counts(src)
            entry["hypothesis"] = hypothesis_text(src)
            if entry["status"] is None:
                staged = any((f / "config.csv").exists() and (f / "data_dir").exists() for f in folders)
                entry["status"] = "staged" if staged else "pending"
        out.append(entry)
    return out


def main() -> int:
    site_name = load_site_name()
    domains = load_domains()
    domain_data = {slug: find_studies(slug) for slug, _, _ in domains}

    DOCS_ROOT.mkdir(parents=True, exist_ok=True)
    (DOCS_ROOT / "index.html").write_text(render_top_level([
        {"name": name, "slug": slug, "desc": desc, "n_studies": len(domain_data[slug])}
        for slug, name, desc in domains
    ], site_name=site_name))
    print(f"Wrote {DOCS_ROOT / 'index.html'}", file=sys.stderr)

    for slug, name, _ in domains:
        studies = domain_data[slug]
        if not studies:
            continue
        domain_docs_dir = DOCS_ROOT / slug
        domain_docs_dir.mkdir(parents=True, exist_ok=True)
        (domain_docs_dir / "index.html").write_text(
            render_domain_index(name, studies, site_name=site_name)
        )
        print(f"Wrote {domain_docs_dir / 'index.html'} ({len(studies)} studies)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
