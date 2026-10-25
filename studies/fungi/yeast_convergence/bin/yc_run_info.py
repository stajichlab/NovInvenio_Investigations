"""Pipeline facts for the yeast_convergence provenance and figure titles, read from the runs.

Each contrast's Nextflow report (results/<id>/nextflow_log/*report.html) names the main.nf
that ran (for example .../NovInvenio-worktrees/main-d03f5a2/main.nf) and the Nextflow
version. The commit comes from `git rev-parse HEAD` of that checkout; the directory name
also carries a short hash, and `agrees` says whether the two match now.
"""
from __future__ import annotations

import csv
import glob
import re
import subprocess
from pathlib import Path

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER = NII / "studies/fungi/yeast_convergence"
MAIN_NF = re.compile(r"(/[^\s\"'<>,]*?/main\.nf)")
NF_VERSION = re.compile(r"Nextflow version</dt>\s*<dd[^>]*>version ([^,<]+),")


def contrast_ids() -> list[str]:
    return [r["id"] for r in csv.DictReader(open(FOLDER / "study_index.tsv"), delimiter="\t")]


def pipeline_info() -> dict:
    mains, versions = set(), set()
    for cid in contrast_ids():
        for rep in glob.glob(str(NII / "results" / cid / "nextflow_log" / "*report.html")):
            text = Path(rep).read_text(errors="replace")
            mains.update(m for m in MAIN_NF.findall(text) if "NovInvenio" in m)
            versions.update(NF_VERSION.findall(text))
    roots = sorted({str(Path(m).parent) for m in mains})
    heads, agree = {}, {}
    for root in roots:
        r = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True)
        heads[root] = r.stdout.strip() or None
        short = re.search(r"-([0-9a-f]{7,40})$", Path(root).name)
        agree[root] = bool(heads[root] and short and heads[root].startswith(short.group(1)))
    return {"main_nf_roots": roots, "heads_now": heads, "dir_name_matches_head": agree,
            "nextflow": sorted(versions)}


def label() -> str:
    """Short text for figure titles, from the directory name of the main.nf that ran."""
    info = pipeline_info()
    names = sorted({re.sub(r"^main-", "", Path(r).name) for r in info["main_nf_roots"]})
    return "NovInvenio " + "/".join(names) if names else "NovInvenio (version not recorded)"
