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
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from site_pages import render_domain_index, render_top_level  # noqa: E402

REPO_ROOT = Path(__file__).parent.parent
STUDIES_ROOT = REPO_ROOT / "studies"
DOCS_ROOT = REPO_ROOT / "docs"

# Domain registry (DESIGN.md Sec 3/8) -- fixed set, populated over time.
DOMAINS = [
    ("fungi", "Fungal", "Lineage-specific gene novelty/loss in fungi (Ascomycota, Basidiomycota, early-diverging lineages)."),
    ("animal", "Animal", "Not yet populated."),
    ("plant", "Plant", "Not yet populated."),
    ("bacteria", "Bacteria", "Not yet populated."),
    ("other", "Other", "Not yet populated."),
]


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
    domain_data = {slug: find_studies(slug) for slug, _, _ in DOMAINS}

    DOCS_ROOT.mkdir(parents=True, exist_ok=True)
    (DOCS_ROOT / "index.html").write_text(render_top_level([
        {"name": name, "slug": slug, "desc": desc, "n_studies": len(domain_data[slug])}
        for slug, name, desc in DOMAINS
    ]))
    print(f"Wrote {DOCS_ROOT / 'index.html'}", file=sys.stderr)

    for slug, name, _ in DOMAINS:
        studies = domain_data[slug]
        if not studies:
            continue
        domain_docs_dir = DOCS_ROOT / slug
        domain_docs_dir.mkdir(parents=True, exist_ok=True)
        (domain_docs_dir / "index.html").write_text(render_domain_index(name, studies))
        print(f"Wrote {domain_docs_dir / 'index.html'} ({len(studies)} studies)", file=sys.stderr)
        for s in studies:
            (domain_docs_dir / s["slug"]).mkdir(parents=True, exist_ok=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
