#!/usr/bin/env python3
"""Dry run of the whole yeast_convergence setup: print what each step would
create or change, per contrast, without writing, fetching or submitting anything.

Steps covered (PLAN.md section 3):
  render   bin/yc_batches.py render: species.csv, run_params.txt, publish.yaml
           (rendered to a temporary directory and compared with what is on disk)
  resolve  bin/ni resolve: rows with a blank Protein_Source; for each, whether
           the lab's local genome lists hold a copy (use as local_faa if
           resolve finds no public proteome)
  fetch    bin/ni fetch (bin/build_study_config.py): config.csv,
           DATA_MANIFEST.yaml, data_dir/{pep,dna,gff3}/, annotations/, using
           that script's own file-naming rules
  run      bin/run_all.sh: submit or skip, results/<id>/, launch dir, logs
  enrich   bin/run_enrichment_all.sh: analysis/enrichment/outputs/{novelty,loss}/
  publish  bin/sync_reports.sh: docs/fungi/<id>/pairwise-diamond/

Usage:
  studies/fungi/yeast_convergence/bin/setup_dryrun.py [--ids yc_a yc_b] [--markdown OUT.md]
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
import tempfile
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER_REL = Path("studies/fungi/yeast_convergence")
LOCAL_LISTS = ("config_support/1KFG_samples.csv", "config_support/Chaetothyriales_samples.csv")
RENDERED = ("species.csv", "run_params.txt", "publish.yaml")


def load_yc(nii_root: Path):
    path = nii_root / FOLDER_REL / "bin" / "yc_batches.py"
    spec = importlib.util.spec_from_file_location("yc_batches", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def local_copies(nii_root: Path) -> dict[str, list[str]]:
    """Species -> ['<list file>: <Protein filename>', ...] from the lab's genome lists."""
    found: dict[str, list[str]] = {}
    for rel in LOCAL_LISTS:
        path = nii_root / rel
        if not path.exists():
            continue
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("Protein"):
                    found.setdefault(row["Species"], []).append(f"{Path(rel).name}: {row['Protein']}")
    return found


def render_plan(yc, batches, folder: Path, known) -> dict[str, dict[str, str]]:
    """Per contrast, per rendered file: create / unchanged / differs (render keeps it)."""
    plan: dict[str, dict[str, str]] = {}
    with tempfile.TemporaryDirectory() as tmp:
        yc.render(batches, Path(tmp), known, force=True)
        for sub in sorted(Path(tmp).glob("yc_*")):
            plan[sub.name] = {}
            for fname in RENDERED:
                new, old = sub / fname, folder / sub.name / fname
                if not old.exists():
                    plan[sub.name][fname] = "create"
                elif old.read_text() == new.read_text():
                    plan[sub.name][fname] = "unchanged"
                else:
                    plan[sub.name][fname] = "differs from YAML (render keeps it; --force to overwrite)"
    return plan


def fetch_files(row: dict) -> list[str]:
    """Files build_study_config.py writes for one resolved row (its naming rules)."""
    src, acc, taxid, short = row["Protein_Source"], row["Protein_Accession"], row["Taxon_ID"], row["Short"]
    out = []
    if src == "uniprot":
        stem = f"{acc}_{taxid}"
        out += [f"data_dir/pep/{stem}.pep.fa", f"annotations/{acc}.tsv.gz"]
    else:
        stem = short
        out.append(f"data_dir/pep/{stem}.pep.fa")
    if row["Genome_Source"] in ("ncbi", "local_genome"):
        out.append(f"data_dir/dna/{stem}.dna.fa")
    if row["Genome_Source"] == "ncbi" and row.get("GFF3_Source") != "none":
        out.append(f"data_dir/gff3/{stem}.gff3 (if the NCBI package has one)")
    return out


def cache_paths(row: dict) -> list[str]:
    out = []
    if row["Protein_Source"] == "uniprot":
        out.append(f"data/uniprot/{row['Protein_Accession']}/")
    if row["Genome_Source"] == "ncbi":
        out.append(f"data/ncbi/{row['Genome_Accession']}/")
    return out


def plan_contrast(sid: str, folder: Path, render: dict[str, str], local: dict[str, list[str]]) -> dict:
    with open(folder / sid / "species.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    resolved = [r for r in rows if r["Protein_Source"]]
    pending = [r for r in rows if not r["Protein_Source"]]
    files, caches = ["config.csv", "DATA_MANIFEST.yaml"], set()
    for r in resolved:
        files += fetch_files(r)
        caches.update(cache_paths(r))
    if any(r["Protein_Source"] not in ("", "uniprot") for r in rows):
        files.append("analysis/pfam_hmmsearch/ (Pfam jobs for the non-UniProt proteomes)")
    elif pending:
        files.append("analysis/pfam_hmmsearch/ (only if a pending row resolves to NCBI or a local file)")
    return {
        "id": sid, "render": render, "rows": rows, "resolved": resolved, "pending": pending,
        "pending_local": {r["Species"]: local.get(r["Species"], []) for r in pending},
        "fetch_files": files, "caches": sorted(caches),
        "run": "submit" if not pending else f"skip ({len(pending)} unresolved)",
    }


def report(plans: list[dict], index: dict[str, dict]) -> str:
    out = ["# yeast_convergence setup: dry run", "",
           "Nothing below has been written, fetched or submitted. Paths are relative to",
           "`studies/fungi/yeast_convergence/<id>/` unless they start with `results/`,",
           "`data/`, `docs/` or `logs/` (repository root).", ""]
    n_rows = sum(len(p["rows"]) for p in plans)
    n_pend = sum(len(p["pending"]) for p in plans)
    n_local = sum(1 for p in plans for v in p["pending_local"].values() if v)
    n_submit = sum(p["run"] == "submit" for p in plans)
    out += [f"- {len(plans)} contrasts, {n_rows} species rows ({n_rows - n_pend} resolved, {n_pend} need `ni resolve`; "
            f"{n_local} of those have a local copy in the lab's genome lists)",
            f"- `run_all.sh` would submit {n_submit} job(s) now and skip {len(plans) - n_submit} until resolved",
            "", "| id | role | rendered files | species resolved | run_all.sh |", "|---|---|---|---|---|"]
    for p in plans:
        rstat = ", ".join(sorted(set(p["render"].values())))
        out.append(f"| `{p['id']}` | {index[p['id']]['role']} | {rstat} | "
                   f"{len(p['resolved'])}/{len(p['rows'])} | {p['run']} |")
    for p in plans:
        sid = p["id"]
        out += ["", f"## `{sid}`", "", "```", f"studies/fungi/yeast_convergence/{sid}/"]
        for fname, st in p["render"].items():
            out.append(f"├── {fname:<16} [render: {st}]")
        out.append("│   (after ni fetch)")
        for f in p["fetch_files"]:
            out.append(f"├── {f}")
        if p["pending"]:
            out.append(f"├── ... plus pep/dna/annotation files for {len(p['pending'])} rows once resolved")
        out += ["│   (after run_enrichment_all.sh)",
                "└── analysis/enrichment/outputs/{novelty,loss}/proteome_clustered/", "",
                "repository root:",
                f"  results/{sid}/                          (pipeline output)",
                f"  .nf_launch/fungi/yeast_convergence/{sid}/  (Nextflow launch dir)",
                f"  logs/slurm/nf_{sid}_<jobid>.out|err",
                f"  docs/fungi/{sid}/pairwise-diamond/       (report.html, run.json)"]
        for c in p["caches"]:
            out.append(f"  {c:<38} (download cache, gitignored)")
        out.append("```")
        if p["pending"]:
            out += ["", "Needs `ni resolve`:", ""]
            for sp, hits in p["pending_local"].items():
                where = f"local copy: {'; '.join(hits[:2])}" + (f" (+{len(hits) - 2} more)" if len(hits) > 2 else "") \
                    if hits else "no local copy listed"
                out.append(f"- *{sp}*: {where}")
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nii-root", type=Path, default=NII_ROOT)
    ap.add_argument("--ids", nargs="*", help="only these contrasts (default: all in study_index.tsv)")
    ap.add_argument("--markdown", type=Path, help="also write the report here")
    args = ap.parse_args(argv)

    yc = load_yc(args.nii_root)
    folder = args.nii_root / FOLDER_REL
    batches = yc.load_batches(folder / "batches")
    problems = yc.check(batches, yc.load_growth_forms(args.nii_root / "config_support/traits/trait_definitions.yaml"))
    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        return 1
    known = yc.known_sources(args.nii_root / "studies", exclude=folder)
    renders = render_plan(yc, batches, folder, known)
    with open(folder / "study_index.tsv", newline="") as fh:
        index = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    ids = args.ids or list(index)
    unknown = [i for i in ids if i not in index]
    if unknown:
        print(f"ERROR: unknown contrast id(s): {unknown}", file=sys.stderr)
        return 1
    local = local_copies(args.nii_root)
    plans = [plan_contrast(i, folder, renders[i], local) for i in ids]
    text = report(plans, index)
    print(text)
    if args.markdown:
        args.markdown.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
