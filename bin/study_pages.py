#!/usr/bin/env python3
"""Study/run page helpers for the docs/ gallery (lib/study_runs.py).

Subcommands:
  target <domain>/<folder>
      Print "<study> <run>" from studies/<domain>/<folder>/publish.yaml.
      Exit 3 when publish.yaml is missing (folder not published), exit 4 when
      it names no run (main.nf folders need one). For bin/sync_reports.sh.
  record-run <domain>/<folder> --report <file> [--current]
      Write docs/<domain>/<study>/<run>/run.json for a main.nf folder
      (params from run_params.txt, species count from species.csv, sha256 of
      --report), then rebuild the study pages.
  set-current <domain>/<study> <run>
      Mark <run> as the study's current run and rebuild the study pages.
  rebuild <domain>/<study>
      Rebuild docs/<domain>/<study>/{report,index}.html from its run.json files.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from study_runs import load_publish_target, rebuild_study_pages, set_current  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent


def split2(value: str, what: str) -> tuple[str, str]:
    parts = value.split("/")
    if len(parts) != 2 or not all(parts):
        raise SystemExit(f"ERROR: {what} must be <domain>/<name>, got {value!r}")
    return parts[0], parts[1]


def read_params(path: Path) -> list[str]:
    """Non-comment, non-blank lines of a run_params.txt, stripped."""
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def count_species(path: Path) -> int | None:
    if not path.exists():
        return None
    with open(path, newline="") as fh:
        return sum(1 for _ in csv.DictReader(fh))


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_target(args, root: Path) -> int:
    domain, folder = split2(args.folder, "folder")
    target = load_publish_target(root / "studies" / domain / folder)
    if target is None:
        return 3
    if target.run is None:
        return 4
    print(f"{target.study} {target.run}")
    return 0


def cmd_record_run(args, root: Path) -> int:
    domain, folder = split2(args.folder, "folder")
    folder_dir = root / "studies" / domain / folder
    target = load_publish_target(folder_dir)
    if target is None or target.run is None:
        raise SystemExit(f"ERROR: {folder_dir}/publish.yaml missing or has no run")
    study_docs = root / "docs" / domain / target.study
    run_docs = study_docs / target.run
    report = Path(args.report)
    if not report.is_file():
        raise SystemExit(f"ERROR: --report {report} not found")
    meta = {
        "run": target.run,
        "published": args.published or datetime.datetime.now(datetime.UTC).date().isoformat(),
        "pipeline": "main.nf",
        "folder": f"studies/{domain}/{folder}",
        "source_dir": args.source_dir or f"results/{folder}",
        "params": read_params(folder_dir / "run_params.txt"),
        "n_species": count_species(folder_dir / "species.csv"),
        "report_sha256": sha256_of(report),
    }
    if args.note:
        meta["note"] = args.note
    run_docs.mkdir(parents=True, exist_ok=True)
    (run_docs / "run.json").write_text(json.dumps(meta, indent=2) + "\n")
    if args.current:
        set_current(study_docs, target.run)
    runs = rebuild_study_pages(study_docs, domain, target.study)
    print(f"recorded {domain}/{target.study}/{target.run}; study has {len(runs)} run(s)", file=sys.stderr)
    return 0


def cmd_set_current(args, root: Path) -> int:
    domain, study = split2(args.study, "study")
    study_docs = root / "docs" / domain / study
    try:
        set_current(study_docs, args.run)
    except ValueError as e:
        raise SystemExit(f"ERROR: {e}")
    rebuild_study_pages(study_docs, domain, study)
    return 0


def cmd_rebuild(args, root: Path) -> int:
    domain, study = split2(args.study, "study")
    rebuild_study_pages(root / "docs" / domain / study, domain, study)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("target")
    p.add_argument("folder")
    p = sub.add_parser("record-run")
    p.add_argument("folder")
    p.add_argument("--report", required=True)
    p.add_argument("--current", action="store_true")
    p.add_argument("--published", help="override the publish date (YYYY-MM-DD)")
    p.add_argument("--source-dir", help="override source_dir (default results/<folder>)")
    p.add_argument("--note", help="free-text note stored in run.json")
    p = sub.add_parser("set-current")
    p.add_argument("study")
    p.add_argument("run")
    p = sub.add_parser("rebuild")
    p.add_argument("study")
    args = ap.parse_args(argv)
    return {"target": cmd_target, "record-run": cmd_record_run,
            "set-current": cmd_set_current, "rebuild": cmd_rebuild}[args.cmd](args, args.repo_root)


if __name__ == "__main__":
    raise SystemExit(main())
