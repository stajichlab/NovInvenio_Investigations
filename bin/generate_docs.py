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
import re
import sys
import tomllib
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from site_pages import (  # noqa: E402
    GROUP_PAGE_MARKER, render_domain_index, render_group_index, render_top_level,
)
from study_runs import (  # noqa: E402
    load_publish_target, read_current, rebuild_study_pages, run_extra_links, validate_related,
    validate_slug,
)

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


def load_domain_pages() -> dict[str, list[dict]]:
    """slug -> hand-made page cards from conf/domains.yaml (`pages:`, optional)."""
    with open(DOMAINS_YAML) as fh:
        rows = yaml.safe_load(fh)
    out = {}
    for row in rows:
        pages = row.get("pages") or []
        for pg in pages:
            if not (isinstance(pg, dict) and pg.get("title") and pg.get("href")):
                raise SystemExit(f"{DOMAINS_YAML}: {row['slug']}: each page needs a title and href, got {pg!r}")
        out[row["slug"]] = pages
    return out


def load_domain_groups() -> dict[str, list[dict]]:
    """slug -> study groups from conf/domains.yaml (`groups:`, optional): [{slug, title, desc,
    member_label, page: {title, href, desc}}]. A study is shown inside a group only when its
    domain lists the group here."""
    with open(DOMAINS_YAML) as fh:
        rows = yaml.safe_load(fh)
    out = {}
    for row in rows:
        groups = []
        for g in row.get("groups") or []:
            if not isinstance(g, dict):
                raise SystemExit(f"{DOMAINS_YAML}: {row['slug']}: each group must be a mapping, got {g!r}")
            try:
                slug = validate_slug(g.get("slug"), "group slug")
            except ValueError as e:
                raise SystemExit(f"{DOMAINS_YAML}: {row['slug']}: {e}")
            if any(x["slug"] == slug for x in groups):
                raise SystemExit(f"{DOMAINS_YAML}: {row['slug']}: group {slug!r} is listed twice")
            page = g.get("page")
            if page is not None:
                if not (isinstance(page, dict) and page.get("title") and page.get("href")):
                    raise SystemExit(f"{DOMAINS_YAML}: {row['slug']}/{slug}: page needs a title and href, got {page!r}")
                try:  # a path under docs/<domain>/, like a publish.yaml `related:` link
                    validate_related([page], f"{DOMAINS_YAML}: {row['slug']}/{slug}: page")
                except ValueError as e:
                    raise SystemExit(str(e))
            groups.append({"slug": slug, "title": g.get("title") or slug.replace("_", " "),
                           "desc": g.get("desc") or "", "member_label": g.get("member_label") or "studies",
                           "page": page})
        out[row["slug"]] = groups
    return out


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


def study_folders(domain_dir: Path) -> list[Path]:
    """Folders under studies/<domain>/ that hold a species.csv: direct children,
    plus the children of a grouping folder that has none itself (e.g.
    studies/fungi/yeast_convergence/<contrast>/)."""
    out = []
    for folder in sorted(p for p in domain_dir.iterdir() if p.is_dir()):
        if (folder / "species.csv").exists():
            out.append(folder)
        else:
            out.extend(sorted(c for c in folder.iterdir()
                              if c.is_dir() and (c / "species.csv").exists()))
    return out


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
            # Optional pages of the current run (e.g. enrichment.html) that
            # exist on disk; the card links to each one.
            extras = [(f"{study_docs.name}/{current}/{fname}", label)
                      for fname, label in run_extra_links(study_docs / current)] if current else []
            by_study[study_docs.name] = {
                "name": study_docs.name, "slug": study_docs.name,
                "status": "complete", "n_runs": len(runs), "extras": extras,
                "href": f"{study_docs.name}/{current}/report.html" if current
                        else f"{study_docs.name}/report.html",
                "hypothesis": "", "n_ingroup": None, "n_outgroup": None,
                "folders": [],
            }

    domain_dir = STUDIES_ROOT / domain_slug
    if domain_dir.exists():
        for folder in study_folders(domain_dir):
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
            # Group: publish.yaml `group:` (null = none), else the grouping folder that
            # holds this folder (studies/<domain>/<group>/<folder>/), else none.
            nested = folder.parent != domain_dir
            entry.setdefault("_groups", set()).add(
                target.group if target.group_set else (folder.parent.name if nested else None))

    out = []
    for name in sorted(by_study):
        entry = by_study[name]
        folders = entry.pop("folders")
        groups = entry.pop("_groups", set())
        if len(groups) > 1:
            raise SystemExit(f"study {name!r}: its folders give different groups: "
                             f"{sorted(str(g) for g in groups)}. Set the same `group:` in each publish.yaml.")
        entry["group"] = next(iter(groups)) if groups else None
        if folders:
            src = next((f for f in folders if f.name == name), folders[0])
            entry["n_ingroup"], entry["n_outgroup"] = species_counts(src)
            entry["hypothesis"] = hypothesis_text(src)
            if entry["status"] is None:
                staged = any((f / "config.csv").exists() and (f / "data_dir").exists() for f in folders)
                entry["status"] = "staged" if staged else "pending"
        out.append(entry)
    return out


def split_groups(domain_slug: str, studies: list[dict], groups_cfg: dict[str, dict],
                 pages: list[dict]) -> tuple[list[dict], list[dict], list[tuple[dict, list[dict]]]]:
    """Split a domain's studies into cards for the domain page.

    Returns (studies shown directly, group cards, [(group, members)] to render as pages).
    A group page is written only when the domain lists the group in conf/domains.yaml and it
    has at least two members. A group with one member shows that member directly. A group with
    none adds its `page:` to `pages` so the link is not lost. A study whose group is not listed
    stays ungrouped (warning)."""
    names = {s["slug"] for s in studies}
    members: dict[str, list[dict]] = {}
    shown = []
    for s in studies:
        g = s.get("group")
        if g in groups_cfg:
            members.setdefault(g, []).append(s)
        else:
            if g:
                print(f"WARNING: {domain_slug}/{s['slug']}: group {g!r} is not listed under `groups:` "
                      f"in {DOMAINS_YAML.name}; shown as its own card.", file=sys.stderr)
            shown.append(s)
    cards, group_pages = [], []
    for slug, g in groups_cfg.items():
        if slug in names:
            raise SystemExit(f"{domain_slug}: group {slug!r} has the same name as a study; rename one.")
        mem = members.get(slug, [])
        if len(mem) >= 2:
            cards.append({**g, "n_members": len(mem)})
            group_pages.append((g, mem))
            continue
        print(f"WARNING: {domain_slug}: group {slug!r} has {len(mem)} member(s); no group page.",
              file=sys.stderr)
        shown.extend(mem)
        if g.get("page"):
            pages.append(g["page"])
    shown.sort(key=lambda s: s["name"])
    return shown, cards, group_pages


def _title(path: Path) -> str | None:
    m = re.search(r"<title>(.*?)</title>", path.read_text(errors="replace"), re.S)
    return m.group(1).strip() if m else None


def check_group_target(domain_docs_dir: Path, group: dict) -> None:
    """Refuse to overwrite a file that is not a generated group page. The one exception is the
    old page of the group (for example the comparison page) at the same address: it may be
    replaced once the group's `page:` file exists at its new address with the same <title>."""
    index = domain_docs_dir / group["slug"] / "index.html"
    if not index.exists() or GROUP_PAGE_MARKER in index.read_text(errors="replace"):
        return
    page = group.get("page")
    moved = False
    if page and page["href"] != f"{group['slug']}/index.html":
        new = domain_docs_dir / page["href"]
        moved = new.exists() and _title(new) is not None and _title(new) == _title(index)
    if not moved:
        raise SystemExit(f"{index} exists and is not a generated group page"
                         + (f"; to replace the old page, write {page['href']} first (same <title>)" if page else "") + ".")


def warn_stale_group_pages(domain_docs_dir: Path, rendered: set[str], configured: set[str]) -> None:
    for index in sorted(domain_docs_dir.glob("*/index.html")):
        name = index.parent.name
        if name in rendered or GROUP_PAGE_MARKER not in index.read_text(errors="replace"):
            continue
        why = ("the group now has fewer than 2 members" if name in configured
               else "the group is no longer in conf/domains.yaml")
        print(f"WARNING: {index} is an old group page ({why}); delete it.", file=sys.stderr)


def main() -> int:
    site_name = load_site_name()
    domains = load_domains()
    domain_data = {slug: find_studies(slug) for slug, _, _ in domains}
    domain_pages = load_domain_pages()
    domain_groups = load_domain_groups()

    # Decide and check everything for all domains first, so that a refusal stops the run
    # before any group or domain page is written.
    plans = []
    for slug, name, _ in domains:
        studies = domain_data[slug]
        pages = list(domain_pages.get(slug) or [])
        groups_cfg = {g["slug"]: g for g in domain_groups.get(slug, [])}
        shown, group_cards, group_pages = split_groups(slug, studies, groups_cfg, pages)
        if not studies and not pages:
            continue
        domain_docs_dir = DOCS_ROOT / slug
        for g, _members in group_pages:
            check_group_target(domain_docs_dir, g)
        plans.append((slug, name, studies, pages, shown, group_cards, group_pages, set(groups_cfg)))

    DOCS_ROOT.mkdir(parents=True, exist_ok=True)
    (DOCS_ROOT / "index.html").write_text(render_top_level([
        {"name": name, "slug": slug, "desc": desc, "n_studies": len(domain_data[slug])}
        for slug, name, desc in domains
    ], site_name=site_name))
    print(f"Wrote {DOCS_ROOT / 'index.html'}", file=sys.stderr)

    for slug, name, studies, pages, shown, group_cards, group_pages, configured in plans:
        domain_docs_dir = DOCS_ROOT / slug
        domain_docs_dir.mkdir(parents=True, exist_ok=True)
        for g, members in group_pages:
            gdir = domain_docs_dir / g["slug"]
            gdir.mkdir(parents=True, exist_ok=True)
            (gdir / "index.html").write_text(render_group_index(name, g, members, site_name=site_name))
            print(f"Wrote {gdir / 'index.html'} ({len(members)} members)", file=sys.stderr)
        warn_stale_group_pages(domain_docs_dir, {g["slug"] for g, _ in group_pages}, configured)
        (domain_docs_dir / "index.html").write_text(
            render_domain_index(name, shown, site_name=site_name, pages=pages, groups=group_cards)
        )
        print(f"Wrote {domain_docs_dir / 'index.html'} ({len(studies)} studies, {len(group_cards)} groups)",
              file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
