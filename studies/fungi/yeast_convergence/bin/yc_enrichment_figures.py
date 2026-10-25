#!/usr/bin/env python3
"""Summary figures for every yeast_convergence contrast's enrichment output.

For each contrast in study_index.tsv (that has enrichment output): links
outputs/proteome_clustered -> novelty/proteome_clustered (the name
bin/generate_enrichment_page.py expects for the novelty direction), then runs
zoosporic_dikarya's make_summary_figure.py for novelty and loss. Writes
<id>/analysis/enrichment/figures/<id>[_loss]_summary.{png,svg,pdf} (ephemeral).
A direction whose figure fails (for example an empty novelty table) is reported and skipped.

Usage: pixi run python bin/yc_enrichment_figures.py [id ...]
Absolute paths (study-specific script, NII CLAUDE.md).
"""
import csv, os, subprocess, sys
from pathlib import Path

import yaml

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER = NII / "studies/fungi/yeast_convergence"
# YC_MAKE_FIGURE overrides the script (for example a branch checkout with the TE colouring, PR #46).
MAKE = Path(os.environ.get("YC_MAKE_FIGURE", NII / "studies/fungi/zoosporic_dikarya/analysis/enrichment/make_summary_figure.py"))
PY = NII / ".pixi/envs/default/bin/python3"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from yc_run_info import label as pipeline_label  # noqa: E402

PIPE = pipeline_label()   # read from the Nextflow reports of the runs
FOOT = ("Annotation: {ann}. Universe: all proteins of the target-group proteomes. Family units: candidate "
        "families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). "
        "One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. "
        "Labels: k = candidate units with term / all units with term.")


def focal_species():
    """id -> focal species, from the batch YAMLs (study_index.tsv has no focal column)."""
    out = {}
    for f in sorted((FOLDER / "batches").glob("*.yaml")):
        for s in yaml.safe_load(open(f)).get("studies", []):
            out[s["id"]] = s.get("focal", "")
    return out


def lineage_label(r, focal):
    """Readable lineage for a figure title. Controls have origin 'none (control)' in
    study_index.tsv, so name the control type (from `contrast`) and the focal species."""
    if r["role"] == "pseudo":
        kind = r["contrast"].replace("pseudo_origin_", "").replace("_", " ").replace("pcwde", "plant cell wall enzyme")
        return f"{kind} control, focal {focal} ({r['split_rank']}-level)"
    return f"{r['origin']}"


ROLE_TEXT = {"origin": "yeast origin", "pseudo": "pseudo-origin control", "reversal": "reversal"}


def main():
    index = list(csv.DictReader(open(FOLDER / "study_index.tsv"), delimiter="\t"))
    focal = focal_species()
    want = set(sys.argv[1:])
    for r in index:
        cid = r["id"]
        if want and cid not in want:
            continue
        enr = FOLDER / cid / "analysis/enrichment"
        if not (enr / "outputs/novelty/proteome_clustered").is_dir():
            print(f"skip {cid}: no enrichment output")
            continue
        link = enr / "outputs/proteome_clustered"
        if not link.exists() and not link.is_symlink():
            os.symlink("novelty/proteome_clustered", link)
        sp = list(csv.DictReader(open(FOLDER / cid / "species.csv")))
        non_up = [f"{s['Species']} ({s['Group']})" for s in sp if s["Protein_Source"] != "uniprot"]
        if non_up:
            ann = "UniProt (GO, Pfam, InterPro) for UniProt proteomes; Pfam only (HMMER hmmsearch) for " + ", ".join(non_up)
            annot_label = "GO/Pfam/InterPro (UniProt) or Pfam (hmmsearch)"
        else:
            ann = "UniProt reference proteomes (Pfam, InterPro, GO)"
            annot_label = "UniProt GO/Pfam/InterPro"
        (enr / "figures").mkdir(exist_ok=True)
        for direction, suffix, extra in (("novelty", "", []),
                                         ("loss", "_loss", ["--group", "OUT", "--candidate-label", "loss candidate proteins"])):
            title = (f"{lineage_label(r, focal.get(cid, ''))} [{ROLE_TEXT[r['role']]}]: "
                     f"{r['n_in']} IN vs {r['n_out']} OUT, {direction} enrichment (diamond pairwise, {PIPE})")
            cmd = [str(PY), str(MAKE), "--enrich", str(enr / "outputs" / direction / "proteome_clustered"),
                   "--config", str(FOLDER / cid / "config.csv"), "--title", title,
                   "--footnote", FOOT.format(ann=ann), "--annot-label", annot_label,
                   "--out", str(enr / "figures" / f"{cid}{suffix}_summary"), *extra]
            p = subprocess.run(cmd, capture_output=True, text=True)
            print(f"{cid} {direction}: {'ok' if p.returncode == 0 else 'FAILED'}")
            if p.returncode:
                print("   ", (p.stderr.strip().splitlines() or ["?"])[-1][:200])


if __name__ == "__main__":
    main()
