#!/usr/bin/env python3
"""Move flat docs/<domain>/<folder>/ report sets into the study/run layout.

Spec: notes/superpowers/specs/2026-09-24-study-run-site-layout-design.md
(Migration, steps 1-2). Local only: this script touches no release and
triggers no deploy.

For every studies/<domain>/<folder>/publish.yaml that names a run, and whose
docs/<domain>/<folder>/report.html exists but has not been migrated yet:
  1. Move the report files (KNOWN_ENTRIES) into docs/<domain>/<study>/<run>/.
     git-tracked files move with `git mv`; gitignored files with a plain move.
  2. Write that run's run.json (bin/study_pages.py record-run). `published` is
     the date of the last commit that touched the old report.html, and a
     `note` records the migration and that params were read from
     run_params.txt at migration time (they may differ from the params the
     original run used, if run_params.txt changed since).
  3. Write redirect stubs at the old URLs (index.html, report.html,
     novelties.html, core.html, losses.html, alignment.html) pointing at the
     moved files, and `git add -f` them (the set-level .gitignore rules still
     apply to real reports).
  4. Rebuild the study pages.

Re-runnable. When the run is already migrated (run.json exists) but the old
folder still holds real report files -- e.g. gitignored files in a checkout
that received the migration through a git merge -- only those files are
moved; redirect stubs are recognised and left in place.

Default is a dry run that prints the plan. Pass --apply to act. Set the
current run afterwards with bin/study_pages.py set-current.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from study_runs import (  # noqa: E402
    REDIRECT_MARKER, REDIRECT_MAX_BYTES, load_publish_target, rebuild_study_pages, render_redirect,
)

REPO_ROOT = Path(__file__).resolve().parent.parent

KNOWN_ENTRIES = [
    "report.html", "alignment.html", "novelties.html", "core.html", "losses.html",
    "summary.pdf", "alignments", "loss_alignments", "archive",
    "go_enrichment.html", "domain_enrichment.html",
]
STUB_NAMES = ["report.html", "novelties.html", "core.html", "losses.html", "alignment.html"]


def git(root: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(root), *args], check=check,
                          capture_output=True, text=True)


def is_tracked(root: Path, path: Path) -> bool:
    r = git(root, "ls-files", "--error-unmatch", str(path.relative_to(root)), check=False)
    return r.returncode == 0


def last_commit_date(root: Path, path: Path) -> str | None:
    r = git(root, "log", "-1", "--format=%cs", "--", str(path.relative_to(root)), check=False)
    return r.stdout.strip() or None


def is_stub(path: Path) -> bool:
    return (path.is_file() and path.stat().st_size <= REDIRECT_MAX_BYTES
            and REDIRECT_MARKER in path.read_text(errors="replace"))


def plan(root: Path) -> list[dict]:
    items = []
    for publish in sorted((root / "studies").glob("*/*/publish.yaml")):
        folder_dir = publish.parent
        target = load_publish_target(folder_dir)
        if target is None or target.run is None:
            continue
        domain, folder = folder_dir.parent.name, folder_dir.name
        old = root / "docs" / domain / folder
        new = root / "docs" / domain / target.study / target.run
        if not old.is_dir() or (old / "run.json").exists():
            continue
        entries = [e for e in KNOWN_ENTRIES
                   if (old / e).exists() and not is_stub(old / e)
                   and not (old == new.parent and e == "report.html" and (new / "run.json").exists())]
        leftover = (new / "run.json").exists()
        if not entries or (not leftover and "report.html" not in entries):
            continue
        unknown = sorted(p.name for p in old.iterdir()
                         if p.name not in KNOWN_ENTRIES and p.name != "index.html"
                         and not (p.is_dir() and (p / "run.json").exists()))
        items.append({"domain": domain, "folder": folder, "study": target.study,
                      "run": target.run, "old": old, "new": new,
                      "entries": entries, "unknown": unknown, "leftover": leftover})
    return items


def apply_item(root: Path, it: dict) -> None:
    old, new = it["old"], it["new"]
    if it["leftover"]:
        rel = os.path.relpath(new, old)
        stubs = []
        for name in it["entries"]:
            src, dest = old / name, new / name
            if dest.exists():
                print(f"  WARNING: {dest} exists -- leaving {src} in place", file=sys.stderr)
                continue
            shutil.move(str(src), str(dest))
            if name in STUB_NAMES:
                (old / name).write_text(render_redirect(f"{rel}/{name}"))
                stubs.append(old / name)
        if stubs:
            git(root, "add", "-f", *[str(s.relative_to(root)) for s in stubs])
        return
    published = last_commit_date(root, old / "report.html")
    new.mkdir(parents=True, exist_ok=True)
    for name in it["entries"]:
        src, dest = old / name, new / name
        if is_tracked(root, src) and src.is_file():
            git(root, "mv", str(src.relative_to(root)), str(dest.relative_to(root)))
        else:
            shutil.move(str(src), str(dest))

    cmd = [sys.executable, str(root / "bin" / "study_pages.py"), "--repo-root", str(root),
           "record-run", f"{it['domain']}/{it['folder']}", "--report", str(new / "report.html"),
           "--note", f"migrated from docs/{it['domain']}/{it['folder']}/ on 2026-09-24; "
                     "params read from run_params.txt at migration time"]
    if published:
        cmd += ["--published", published]
    subprocess.run(cmd, check=True)

    same_dir = old == new.parent
    rel = os.path.relpath(new, old)
    stubs = []
    for name in STUB_NAMES:
        if same_dir and name == "report.html":
            continue  # the old set page becomes the study's run list
        if (new / name).exists():
            (old / name).write_text(render_redirect(f"{rel}/{name}"))
            stubs.append(old / name)
    if not same_dir:
        (old / "index.html").write_text(render_redirect(f"{rel}/report.html"))
        stubs.append(old / "index.html")
    if stubs:
        git(root, "add", "-f", *[str(s.relative_to(root)) for s in stubs])

    study_docs = new.parent
    rebuild_study_pages(study_docs, it["domain"], it["study"])
    to_add = [study_docs / "report.html", study_docs / "index.html", new / "run.json"]
    if (study_docs / "study.json").exists():
        to_add.append(study_docs / "study.json")
    git(root, "add", *[str(p.relative_to(root)) for p in to_add])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="act (default: print the plan only)")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    args = ap.parse_args(argv)
    root = args.repo_root.resolve()

    items = plan(root)
    if not items:
        print("nothing to migrate", file=sys.stderr)
        return 0
    for it in items:
        mode = " (leftover files only)" if it["leftover"] else ""
        print(f"{it['domain']}/{it['folder']} -> {it['domain']}/{it['study']}/{it['run']}{mode}: "
              f"{', '.join(it['entries'])}", file=sys.stderr)
        if it["unknown"]:
            print(f"  WARNING: left in place (not a known report file): {', '.join(it['unknown'])}",
                  file=sys.stderr)
    if not args.apply:
        print("dry run -- pass --apply to act", file=sys.stderr)
        return 0
    for it in items:
        apply_item(root, it)
    print(f"migrated {len(items)} folder(s); set current runs with bin/study_pages.py set-current",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
