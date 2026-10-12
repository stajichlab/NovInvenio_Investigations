#!/usr/bin/env python3
"""One-slide summary of the neocallimastigomycota_vs_early_fungi (6 IN) novelty candidates.

Panels: A candidates per IN species, split by UniProt annotation (Pfam/InterPro);
B top Pfam families enriched in the whole candidate set (family-level test);
C expected gut-fungus protein families: share of annotated IN proteins that are candidates.

Inputs: enrich_novelties.py outputs (coverage.tsv, enrichment.tsv.gz), the study
results dir (presence_matrix.uniprot.tsv, candidates.txt). Run with the NII pixi env python.
  make_summary_figure.py --results RESULTS_DIR --enrich OUT_DIR --out PREFIX
"""
import argparse
import csv
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True, type=Path)
ap.add_argument("--enrich", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path, help="output prefix, writes .png .svg .pdf")
ap.add_argument("--top", type=int, default=12)
a = ap.parse_args()

IN = ["Pirofinn", "Anaerobu", "Neocalif", "Neolanat", "Neoconst", "PiroE2sp"]
LABEL = {"Pirofinn": "Piromyces finnis", "Anaerobu": "Anaeromyces robustus", "Neocalif": "Neocallimastix californiae",
         "Neolanat": "N. lanati", "Neoconst": "N. sp. constans", "PiroE2sp": "Piromyces sp. E2"}
SRC = {"Pirofinn": "UniProt", "Anaerobu": "UniProt"}

# --- A: candidates per species ---
cov = {}
for r in csv.DictReader(open(a.enrich / "coverage.tsv"), delimiter="\t"):
    if r["subset"] == "candidates":
        cov[r["species"]] = (int(r["n_proteins"]), int(r["n_any"]))
total = sum(cov[s][0] for s in IN)
n_annot = sum(cov[s][1] for s in IN)

# --- B: top family-level Pfam ---
rows = []
with gzip.open(a.enrich / "enrichment.tsv.gz", "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["scope"] == "whole_set" and r["level"] == "family" and r["vocab"] == "pfam":
            rows.append(r)
rows.sort(key=lambda r: float(r["pvalue"]))
top = rows[:a.top]
REPEAT = ("Ank", "LRR", "AAA-ATPase", "PDDEXK", "Retrotran", "MULE", "LOR")
def kind(r):
    if r["term"] == "PF02013":
        return "dockerin"
    if any(r["name"].startswith(p) for p in REPEAT):
        return "repeat/TE-like"
    return "other"
COL = {"dockerin": "#d95f02", "repeat/TE-like": "#9a9a9a", "other": "#2a78d6"}

# --- C: expected families, annotated IN proteins ---
CTRL = [("PF02013", "Dockerin (CBM_10)"), ("PF00150", "Cellulase (GH5)"), ("PF01341", "GH6"),
        ("PF00759", "GH9"), ("PF02011", "GH48"), ("PF02901", "Pyruvate formate lyase")]
cands = {l.strip() for l in open(a.results / "candidates.txt") if l.strip()}
have = {c: [0, 0, set()] for c, _ in CTRL}
with open(a.results / "presence_matrix.uniprot.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        accs = {x.strip() for x in (r.get("uniprot_pfam_ids") or "").replace(";", "|").replace(",", "|").split("|") if x.strip()}
        key = f'{r["source_proteome"]}::{r["protein_id"]}'
        for c, _ in CTRL:
            if c in accs or any(x.split(".")[0] == c for x in accs):
                have[c][0] += 1
                have[c][1] += key in cands
fig = plt.figure(figsize=(14, 6.2), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.0])

ax = fig.add_subplot(gs[0])
ys = list(range(len(IN)))[::-1]
for y, s in zip(ys, IN):
    n, na = cov[s]
    ax.barh(y, na, color="#2a78d6")
    ax.barh(y, n - na, left=na, color="#c9c9c9")
    ax.text(n + total * 0.005, y, f"{n:,}", va="center", fontsize=9)
ax.set_yticks(ys)
ax.set_yticklabels([f"{LABEL[s]}\n({SRC.get(s, 'NCBI/JGI')})" for s in IN], fontsize=9)
ax.set_xlabel("novelty candidate proteins")
ax.set_title(f"A  {total:,} candidates in 6 IN species", loc="left", fontsize=11, fontweight="bold")
ax.legend(handles=[Patch(color="#2a78d6", label=f"with UniProt annotation ({n_annot:,}, {n_annot/total:.0%})"),
                   Patch(color="#c9c9c9", label="no annotation")], loc="upper right", fontsize=8, frameon=False)
ax.set_xlim(0, max(cov[s][0] for s in IN) * 1.18)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

ax = fig.add_subplot(gs[1])
import math
ys = list(range(len(top)))[::-1]
for y, r in zip(ys, top):
    q = float(r["fdr_bh"]); k = r["k_cand_with_term"]
    ax.barh(y, -math.log10(max(q, 1e-300)), color=COL[kind(r)])
    ax.text(-math.log10(max(q, 1e-300)) + 1, y, f'{r["fold"]}x  (n={k})', va="center", fontsize=8)
ax.set_yticks(ys)
ax.set_yticklabels([f'{r["name"]}' for r in top], fontsize=9)
ax.set_xlabel("-log10 FDR (BH)")
ax.set_title("B  Pfam families enriched in candidates\n(family-level test)", loc="left", fontsize=11, fontweight="bold")
ax.legend(handles=[Patch(color=c, label=k_) for k_, c in COL.items()], loc="lower right", fontsize=8, frameon=False)
ax.set_xlim(0, max(-math.log10(max(float(r["fdr_bh"]), 1e-300)) for r in top) * 1.3)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

ax = fig.add_subplot(gs[2])
ys = list(range(len(CTRL)))[::-1]
for y, (c, lab) in zip(ys, CTRL):
    n, k = have[c][0], have[c][1]
    frac = k / n if n else 0
    ax.barh(y, frac * 100, color="#d95f02" if c == "PF02013" else "#2a78d6")
    ax.text(frac * 100 + 1.5, y, f"{k} of {n}", va="center", fontsize=9)
ax.set_yticks(ys)
ax.set_yticklabels([lab for _, lab in CTRL], fontsize=9)
ax.set_xlim(0, 100)
ax.set_xlabel("% of annotated IN proteins with the domain\nthat are novelty candidates")
ax.set_title("C  Expected gut-fungus families", loc="left", fontsize=11, fontweight="bold")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

fig.suptitle("Neocallimastigomycota vs early-diverging fungi: novelty candidates (6 IN, 13 OUT; diamond, UniProt-derived annotation)",
             fontsize=12)
fig.text(0.01, -0.02,
         "Annotated = Pfam/InterPro from UniProt. The 4 NCBI/JGI proteomes have none, so B and C use Piromyces finnis and Anaeromyces robustus. "
         "A family is a candidate only if absent from all 13 OUT species; one OUT homolog removes it (GH48, PFL).",
         fontsize=8, color="#444")
a.out.parent.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=200, bbox_inches="tight")
for c, lab in CTRL:
    print(c, lab, have[c][:2])
