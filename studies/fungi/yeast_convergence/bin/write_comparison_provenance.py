#!/usr/bin/env python3
"""Write analysis/convergence/PROVENANCE.yaml for a convergence_compare.py run.

Records, for the tables in analysis/convergence/: the SHA-256 of every input
(enrichment tables, study_index.tsv, focus_terms.tsv), the script and its git
commit, the command and parameters, the workflow versions (the nf_NovInvenio checkout and
Nextflow version named in the run reports, Pfam), and the SHA-256 of every output. Also gzips term_recurrence.tsv
(gzip, mtime 0) so the tracked copy is term_recurrence.tsv.gz.

Run after convergence_compare.py:  python3.12 bin/write_comparison_provenance.py
Absolute paths on purpose (see NII CLAUDE.md, study-specific scripts).
"""
import glob, gzip, hashlib, json, re, shutil, subprocess, sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER = NII_ROOT / "studies/fungi/yeast_convergence"
OUT = FOLDER / "analysis/convergence"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from yc_run_info import pipeline_info  # noqa: E402
sys.path.insert(0, str(NII_ROOT))
import yaml  # noqa: E402
from lib.provenance import sha256_of, now_utc_iso  # noqa: E402


def git(root, *a):
    return subprocess.run(["git", "-C", str(root), *a], capture_output=True, text=True, check=True).stdout.strip()


def main():
    plain = OUT / "term_recurrence.tsv"
    gz = OUT / "term_recurrence.tsv.gz"
    if plain.exists():
        with open(plain, "rb") as fi, open(gz, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
    inputs = {}
    ids = [l.split("\t")[0] for l in (FOLDER / "study_index.tsv").read_text().splitlines()[1:]]
    for cid in ids:
        for d in ("novelty", "loss"):
            for f in sorted(glob.glob(str(FOLDER / cid / "analysis/enrichment/outputs" / d / "proteome_clustered/*"))):
                p = Path(f)
                if p.name.startswith("enrichment_") or p.name == "PROVENANCE.yaml":
                    inputs[str(p.relative_to(NII_ROOT))] = sha256_of(p)
    for f in ("study_index.tsv", "focus_terms.tsv"):
        inputs[str((FOLDER / f).relative_to(NII_ROOT))] = sha256_of(FOLDER / f)
    record_path = OUT / "run_record.json"
    if not record_path.is_file():
        raise SystemExit(f"ERROR: {record_path} missing; run convergence_compare.py first")
    run = json.loads(record_path.read_text())      # written by convergence_compare.py itself
    script = Path(run["script"])                   # the script that ran (may be a branch checkout)
    top = Path(git(script.parent, "rev-parse", "--show-toplevel"))
    rel = str(script.relative_to(top))
    pipe = pipeline_info()                          # from the Nextflow reports of the runs
    # Pfam and enrichment script come from the enrichment runs' own PROVENANCE.yaml files.
    enr = yaml.safe_load(open(next(iter(sorted(glob.glob(str(
        FOLDER / "yc_*/analysis/enrichment/outputs/novelty/proteome_clustered/PROVENANCE.yaml")))))))
    outputs = {p.name: sha256_of(p) for p in sorted(OUT.iterdir())
               if p.is_file() and p.name not in ("PROVENANCE.yaml", "term_recurrence.tsv")}
    rec = {
        "generated": now_utc_iso(),
        "class": "derived table (tracked, class 2); inputs are ephemeral and not tracked",
        "source_url": "none; derived from this repo's yeast_convergence enrichment output",
        "license": "derived; inputs follow their own licenses (UniProt CC BY 4.0, BFD lab annotations unpublished)",
        "derived_by": {
            "script": rel,
            "script_sha256": sha256_of(script),
            "script_last_commit": git(top, "log", "-1", "--format=%H", "--", rel),
            "repo_head_when_recorded": git(top, "rev-parse", "HEAD"),
            "script_checkout": str(top),
            "command": " ".join([run["executable"], *run["argv"]]),
            "parameters": run["parameters"],
            "run_finished_utc": run["finished_utc"],
            "python_of_run": run["python"],
            "source": "analysis/convergence/run_record.json, written by convergence_compare.py during the run",
        },
        "workflow_versions": {
            "pipeline_main_nf_roots": pipe["main_nf_roots"],          # named in the Nextflow reports
            "pipeline_checkout_head_now": pipe["heads_now"],
            "pipeline_dir_name_matches_head_now": pipe["dir_name_matches_head"],
            "nextflow": pipe["nextflow"],
            "pipeline_run_mode": "pairwise diamond (--run_tool diamond)",
            "enrichment_script": enr.get("script", {}).get("path"),
            "enrichment_script_sha256": enr.get("script", {}).get("sha256"),
            "pfam_annotation_file": (enr.get("pfam_dat") or {}).get("path"),
            "python_of_this_record": sys.version.split()[0],
        },
        "inputs_sha256": inputs,
        "outputs_sha256": outputs,
        "output_note": "term_recurrence.tsv.gz is gzip of term_recurrence.tsv (mtime 0); the plain file is ignored by git. Its sha256 is of the .gz.",
        "scope_note": "All contrasts that have enrichment tables (see run_record.json missing_tables). GO/InterPro only over all-UniProt contrasts.",
    }
    with open(OUT / "PROVENANCE.yaml", "w") as fh:
        yaml.safe_dump(rec, fh, sort_keys=False, width=120)
    print(f"{len(inputs)} inputs, {len(outputs)} outputs, nextflow={pipe['nextflow']}")


if __name__ == "__main__":
    main()
