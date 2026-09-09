#!/usr/bin/env python3
"""Generate the two-tier docs/ gallery (DESIGN.md Sec 8): top-level domain cards,
each linking to a per-domain page of study cards.

Scans studies/<domain>/<set_name>/ for species.csv (existence = a defined study).
A study's status is:
  - "pending"  -- species.csv exists but studies/<domain>/<set>/config.csv doesn't
                  (bin/build_study_config.py hasn't been run yet)
  - "staged"   -- config.csv/data_dir exist but docs/<domain>/<set>/report.html
                  doesn't (data pulled, pipeline not yet run/published)
  - "complete" -- docs/<domain>/<set>/report.html exists (a real nf_NovInvenio run
                  has been published here -- not built by this script, this script
                  only creates the gallery pages that link to it)

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


def study_status(domain_slug: str, set_name: str) -> str:
    study_dir = STUDIES_ROOT / domain_slug / set_name
    docs_dir = DOCS_ROOT / domain_slug / set_name
    if (docs_dir / "report.html").exists():
        return "complete"
    if (study_dir / "config.csv").exists() and (study_dir / "data_dir").exists():
        return "staged"
    return "pending"


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
    domain_dir = STUDIES_ROOT / domain_slug
    if not domain_dir.exists():
        return []
    out = []
    for study_dir in sorted(domain_dir.iterdir()):
        if not study_dir.is_dir() or not (study_dir / "species.csv").exists():
            continue
        n_in, n_out = species_counts(study_dir)
        out.append({
            "name": study_dir.name,
            "slug": study_dir.name,
            "hypothesis": hypothesis_text(study_dir),
            "n_ingroup": n_in,
            "n_outgroup": n_out,
            "status": study_status(domain_slug, study_dir.name),
        })
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
        for s in studies:
            (domain_docs_dir / s["slug"]).mkdir(parents=True, exist_ok=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
