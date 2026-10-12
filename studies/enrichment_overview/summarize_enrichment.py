#!/usr/bin/env python3
"""Cross-study summary of enrich_novelties_v2.py primary runs (whole_set, family level).

Writes enrichment_summary_all_studies.tsv (one row per study x vocabulary) and
top_pfam_all_studies.tsv (top N Pfam rows per study), and prints both.
Repeat/TE flag = same Pfam-ID-prefix regex as make_summary_figure.py (heuristic, not a test).
Run: .pixi/envs/default/bin/python3 -I summarize_enrichment.py --outdir DIR STUDY_OUTDIR...
STUDY_OUTDIR = .../analysis/enrichment/outputs/proteome_clustered
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

import yaml

csv.field_size_limit(10**9)
REPEAT = re.compile(r"^(Ank|LRR|WD40|TPR|HEAT|Kelch|Sel1|PPR|RVT|rve|Transposase|DDE|MULE|HTH_Tnp|Retrotrans|"
                    r"zf-CCHC|PDDEXK|Helitron|Tc5|CENP-B|DUF4219|Chromo|Exo_endo_phos|Mito_fiss_reg)", re.I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("dirs", nargs="+", type=Path)
    a = ap.parse_args()
    summ, top = [], []
    for d in a.dirs:
        prov = yaml.safe_load(open(d / "PROVENANCE.yaml"))
        study = prov["study"]
        cc = prov.get("candidate_check", {})
        tests = [r for r in csv.DictReader(open(d / "tests.tsv"), delimiter="\t")
                 if r["scope"] == "whole_set" and r["level"] == "family"]
        rows = {"domains": [], "go": []}
        for role in rows:
            with gzip.open(d / f"enrichment_{role}.tsv.gz", "rt") as fh:
                rows[role] = [r for r in csv.DictReader(fh, delimiter="\t")
                              if r["scope"] == "whole_set" and r["level"] == "family"]
        for t in tests:
            v = t["vocab"]
            sig = [r for r in rows["go" if v == "go" else "domains"] if r["vocab"] == v and float(r["fdr_bh"]) < 0.05]
            sig_k_lt3 = sum(1 for r in sig if int(r["k_cand_with_term"]) < 3)
            prot_ids = set()
            for r in sig:
                prot_ids.update(x for x in r["candidate_ids"].split(";") if x)
            rec = {"study": study, "vocab": v, "n_cand_units": t["n_cand_units"],
                   "N_cand_annotated": t["N_cand_annotated"], "frac_cand_units_annotated": t["frac_cand_annotated"],
                   "M_units_annotated": t["M_units_annotated"], "n_terms_tested": t["n_terms_tested"],
                   "n_fdr_lt_0.05": t["n_fdr_lt_0.05"], "n_fdr_lt_0.05_k_ge_3": t["n_fdr_lt_0.05_k_ge_3"],
                   "n_sig_k_lt_3": sig_k_lt3, "n_candidate_proteins_in_sig_terms": len(prot_ids),
                   "n_sig_repeat_TE_like": (sum(1 for r in sig if REPEAT.match(r["name"].split(":")[0]))
                                            if v == "pfam" else ""),
                   "candidate_check_symdiff": cc.get("symmetric_difference", ""),
                   "candidates_txt": cc.get("candidates_txt", "")}
            summ.append(rec)
        pf = sorted((r for r in rows["domains"] if r["vocab"] == "pfam"), key=lambda r: float(r["pvalue"]))
        for i, r in enumerate(pf[:a.top], 1):
            ids = [x for x in r["candidate_ids"].split(";") if x]
            top.append({"study": study, "rank": i, "term": r["term"], "name": r["name"],
                        "k_cand_units": r["k_cand_with_term"], "N_cand": r["N_cand_annotated"],
                        "n_bg_units": r["n_bg_with_term"], "M": r["M_bg_annotated"], "fold": r["fold"],
                        "fdr_bh": f"{float(r['fdr_bh']):.3g}", "n_candidate_proteins": len(ids) if ids else "",
                        "n_species": len({x.split('::')[0] for x in ids}) if ids else "",
                        "repeat_TE_like": bool(REPEAT.match(r["name"].split(":")[0]))})
    a.outdir.mkdir(parents=True, exist_ok=True)
    for name, data in (("enrichment_summary_all_studies.tsv", summ), ("top_pfam_all_studies.tsv", top)):
        with open(a.outdir / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(data)
        print(f"== {name}")
        w = csv.DictWriter(sys.stdout, fieldnames=list(data[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(data)


if __name__ == "__main__":
    main()
