#!/usr/bin/env python3
"""One-slide summary of novelty-candidate enrichment for one study (enrich_novelties_v2.py output).

Adapted from studies/fungi/neocallimastigomycota_vs_early_fungi/analysis/enrichment/make_summary_figure.py.
Panels:
  A  candidates per IN species, split by UniProt annotation (any GO/Pfam/InterPro)
  B  top Pfam families in the whole candidate set, family-level test (PRIMARY result)
  C  top GO terms, whole set, family-level test (SECONDARY; IEA-only terms marked)

  make_summary_figure.py --enrich OUT_DIR --config CONFIG_CSV --title TEXT --out PREFIX
Run with the NII pixi env python.
"""
import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# Candidate-ID list columns in the enrichment tables can exceed the default 128 KB field limit.
csv.field_size_limit(1 << 30)

ap = argparse.ArgumentParser()
ap.add_argument("--enrich", required=True, type=Path)
ap.add_argument("--config", required=True, type=Path)
ap.add_argument("--title", required=True)
ap.add_argument("--footnote", default="")
ap.add_argument("--out", required=True, type=Path, help="output prefix; writes .png .svg .pdf")
ap.add_argument("--top", type=int, default=12)
ap.add_argument("--top-go", type=int, default=8)
ap.add_argument("--min-k", type=int, default=3, help="display filter: candidate units with the term")
ap.add_argument("--group", choices=["IN", "OUT"], default="IN",
                help="config group the candidates come from (panel A rows): IN for novelty "
                     "(default), OUT for loss-direction enrichment (enrich_novelties_v2.py --direction loss)")
ap.add_argument("--candidate-label", default="novelty candidate proteins", help="panel A x-axis label")
a = ap.parse_args()

cfg = list(csv.DictReader(open(a.config)))
# IN = the candidates' source group (named IN for the novelty default; OUT rows with --group OUT).
IN = [r["Short"] for r in cfg if r["GROUP"] == a.group]
LABEL = {r["Short"]: r["Species"] for r in cfg}

cov = {}
for r in csv.DictReader(open(a.enrich / "coverage.tsv"), delimiter="\t"):
    cov[(r["species"], r["subset"])] = r
total = sum(int(cov[(s, "candidates")]["n_proteins"]) for s in IN)
n_annot = sum(int(cov[(s, "candidates")]["n_any"]) for s in IN)


def load(fn, vocab):
    rows = []
    with gzip.open(a.enrich / fn, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r["scope"] == "whole_set" and r["level"] == "family" and r["vocab"] == vocab
                    and int(r["k_cand_with_term"]) >= a.min_k and float(r["fdr_bh"]) < 0.05):
                rows.append(r)
    rows.sort(key=lambda r: float(r["pvalue"]))
    return rows


pf_all = load("enrichment_domains.tsv.gz", "pfam")
go_all = load("enrichment_go.tsv.gz", "go")
pf = pf_all[:a.top]
# GO display: skip a term whose candidate-unit set equals a better-ranked shown term
# (parent/child terms with identical members). Display only; no test changes.
go, seen = [], set()
for r in go_all:
    key = r["candidate_ids"] or r["term"]
    if key in seen:
        continue
    seen.add(key)
    go.append(r)
    if len(go) == a.top_go:
        break
tests = {(r["scope"], r["level"], r["vocab"]): r for r in csv.DictReader(open(a.enrich / "tests.tsv"), delimiter="\t")}

# Heuristic label only (Pfam ID prefix). Not a test. Repeats and mobile-element domains.
REPEAT = re.compile(r"^(Ank|LRR|WD40|TPR|HEAT|Kelch|Sel1|PPR|RVT|rve|Transposase|DDE|MULE|HTH_Tnp|Retrotrans|"
                    r"zf-CCHC|PDDEXK|Helitron|Tc5|CENP-B|DUF4219|Chromo|Exo_endo_phos|Mito_fiss_reg)", re.I)


def kind(r):
    return "repeat / TE-like (name heuristic)" if REPEAT.match(r["name"].split(":")[0]) else "other"


COL = {"other": "#2a78d6", "repeat / TE-like (name heuristic)": "#9a9a9a"}

fig = plt.figure(figsize=(16, 6.6), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 1.0])

# --- A ---
ax = fig.add_subplot(gs[0])
ys = list(range(len(IN)))[::-1]
xmax = max(int(cov[(s, "candidates")]["n_proteins"]) for s in IN)
for y, s in zip(ys, IN):
    c = cov[(s, "candidates")]
    n, na = int(c["n_proteins"]), int(c["n_any"])
    ax.barh(y, na, color="#2a78d6", height=0.6)
    ax.barh(y, n - na, left=na, color="#c9c9c9", height=0.6)
    u = cov[(s, "universe")]
    ax.text(n + xmax * 0.02, y, f"{n:,} ({na / n:.0%})\nproteome {float(u['frac_any']):.0%}",
            va="center", fontsize=8)
ax.set_yticks(ys)
ax.set_yticklabels([f"{LABEL[s]}" for s in IN], fontsize=9, style="italic")
ax.set_xlabel(a.candidate_label)
ax.set_xlim(0, xmax * 1.7)
ax.set_ylim(-0.6, len(IN) - 0.4)
ax.set_title(f"A  {total:,} candidates in {len(IN)} {a.group} species", loc="left", fontsize=11, fontweight="bold")
ax.legend(handles=[Patch(color="#2a78d6", label=f"UniProt GO/Pfam/InterPro ({n_annot:,}, {n_annot / total:.0%})"),
                   Patch(color="#c9c9c9", label="no annotation")],
          loc="upper left", bbox_to_anchor=(0, -0.09), fontsize=8, frameon=False)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)


def barpanel(ax, rows, colorf, labelf, title, test):
    ys = list(range(len(rows)))[::-1]
    vals = [-math.log10(max(float(r["fdr_bh"]), 1e-300)) for r in rows]
    vmax = max(vals) if vals else 1
    for y, r, v in zip(ys, rows, vals):
        col, hatch = colorf(r)
        ax.barh(y, v, color=col, hatch=hatch, edgecolor="white" if not hatch else "#555", linewidth=0.5, height=0.7)
        ax.text(v + vmax * 0.02, y, f'{float(r["fold"]):.1f}x  k={r["k_cand_with_term"]}/{r["n_bg_with_term"]}',
                va="center", fontsize=7.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([labelf(r) for r in rows], fontsize=8.5)
    ax.set_xlim(0, vmax * 1.45)
    ax.set_xlabel("-log10 FDR (BH within vocabulary)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=14)
    ax.annotate(f"{test['n_terms_tested']} terms tested; {test['n_fdr_lt_0.05']} at FDR<0.05; "
                f"annotated units: {test['N_cand_annotated']} cand. / {test['M_units_annotated']} all",
                xy=(0, 1), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                fontsize=7.5, color="#444", va="bottom")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def short(s, n=46):
    return s if len(s) <= n else s[:n - 1] + "…"


ax = fig.add_subplot(gs[1])
barpanel(ax, pf, lambda r: (COL[kind(r)], None),
         lambda r: short(r["name"].replace(": ", " – "), 42),
         "B  Pfam families enriched in candidates\n(primary; family units on both sides)",
         tests[("whole_set", "family", "pfam")])
ax.legend(handles=[Patch(color=c, label=k_) for k_, c in COL.items()], loc="upper left",
          bbox_to_anchor=(0, -0.09), ncol=2, fontsize=8, frameon=False)

ax = fig.add_subplot(gs[2])
COMPUTATIONAL = {"IEA", "IBA", "ISS", "ISO", "ISA", "ISM", "IGC", "RCA", "NA"}


def go_style(r):
    codes = {x.split(":")[0] for x in r.get("cand_evidence_codes", "").split(";") if x}
    if r["go_evidence"] == "IEA_only":
        return ("#b8b8b8", "///")
    if codes <= COMPUTATIONAL:
        return ("#c7b3e0", None)
    return ("#7b5aa6", None)


barpanel(ax, go, go_style,
         lambda r: short(r["name"].split(":", 1)[1], 44) + f" [{r['name'][0].upper()}]" if ":" in r["name"] else r["term"],
         "C  GO terms (secondary; UniProt\ncomputational evidence)",
         tests[("whole_set", "family", "go")])
ax.legend(handles=[Patch(facecolor="#b8b8b8", hatch="///", edgecolor="#555", label="IEA only"),
                   Patch(color="#c7b3e0", label="IEA + IBA/ISS (computational)"),
                   Patch(color="#7b5aa6", label="incl. curated/experimental code")],
          loc="upper left", bbox_to_anchor=(0, -0.09), ncol=1, fontsize=8, frameon=False)

fig.suptitle(a.title, fontsize=12.5)
if a.footnote:
    import textwrap
    fig.text(0.005, -0.06, "\n".join(textwrap.wrap(a.footnote, 210)), fontsize=8, color="#444", va="top")
a.out.parent.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=200, bbox_inches="tight")
print(f"Pfam rows shown {len(pf)} of {len(pf_all)} (FDR<0.05, k>={a.min_k}); GO shown {len(go)} of {len(go_all)}")
