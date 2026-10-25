#!/usr/bin/env python3
"""Write analysis/convergence/PROVENANCE.yaml for a convergence_compare.py run.

Records, for the tables in analysis/convergence/: the SHA-256 of every input
(enrichment tables, study_index.tsv, focus_terms.tsv), the script and its git
commit, the command and parameters, the workflow versions (nf_NovInvenio commit,
Nextflow, Pfam), and the SHA-256 of every output. Also gzips term_recurrence.tsv
(gzip, mtime 0) so the tracked copy is term_recurrence.tsv.gz.

Run after convergence_compare.py:  python3.12 bin/write_comparison_provenance.py
Absolute paths on purpose (see NII CLAUDE.md, study-specific scripts).
"""
import glob, gzip, re, shutil, subprocess, sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER = NII_ROOT / "studies/fungi/yeast_convergence"
OUT = FOLDER / "analysis/convergence"
NOVINVENIO_ROOT = Path("/rhome/jstajich/projects/NovInvenio-worktrees/main-d03f5a2")
PFAM_DAT = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat"
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
    nf_versions = set()
    ids = [l.split("\t")[0] for l in (FOLDER / "study_index.tsv").read_text().splitlines()[1:]]
    for cid in ids:
        for d in ("novelty", "loss"):
            for f in sorted(glob.glob(str(FOLDER / cid / "analysis/enrichment/outputs" / d / "proteome_clustered/*"))):
                p = Path(f)
                if p.name.startswith("enrichment_") or p.name == "PROVENANCE.yaml":
                    inputs[str(p.relative_to(NII_ROOT))] = sha256_of(p)
        for rep in glob.glob(str(NII_ROOT / "results" / cid / "nextflow_log/*report.html")):
            m = re.search(r"Nextflow version</dt>\s*<dd[^>]*>version ([^,]+),", Path(rep).read_text())
            if m:
                nf_versions.add(m.group(1))
    for f in ("study_index.tsv", "focus_terms.tsv"):
        inputs[str((FOLDER / f).relative_to(NII_ROOT))] = sha256_of(FOLDER / f)
    script = FOLDER / "bin/convergence_compare.py"
    rel = str(script.relative_to(NII_ROOT))
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
            "script_last_commit": git(NII_ROOT, "log", "-1", "--format=%H", "--", rel),
            "repo_head_when_recorded": git(NII_ROOT, "rev-parse", "HEAD"),
            "command": f"/usr/bin/python3.12 {script} (no arguments)",
            "parameters": {
                "vocab": ["pfam", "pfam_clan", "interpro", "go"], "scope": "whole_set", "level": "family",
                "fdr": 0.05, "min_k": 3, "alpha": 0.05, "min_origins": 3, "min_in_novelty": 2,
                "include_mixed_annotation": False,
                "note": "defaults of convergence_compare.py at the commit above",
            },
        },
        "workflow_versions": {
            "nf_NovInvenio_commit": git(NOVINVENIO_ROOT, "rev-parse", "HEAD"),
            "nf_NovInvenio_path": str(NOVINVENIO_ROOT),
            "nextflow": sorted(nf_versions),
            "pipeline_run_mode": "pairwise diamond (--run_tool diamond)",
            "enrichment_script": "studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py",
            "enrichment_script_sha256": sha256_of(NII_ROOT / "studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py"),
            "pfam_annotation_file": PFAM_DAT,
            "python": sys.version.split()[0],
        },
        "inputs_sha256": inputs,
        "outputs_sha256": outputs,
        "output_note": "term_recurrence.tsv.gz is gzip of term_recurrence.tsv (mtime 0); the plain file is ignored by git. Its sha256 is of the .gz.",
        "scope_note": "Excludes yc_pseudo_lichen (no run). GO/InterPro only over all-UniProt contrasts.",
    }
    with open(OUT / "PROVENANCE.yaml", "w") as fh:
        yaml.safe_dump(rec, fh, sort_keys=False, width=120)
    print(f"{len(inputs)} inputs, {len(outputs)} outputs, nextflow={sorted(nf_versions)}")


if __name__ == "__main__":
    main()
