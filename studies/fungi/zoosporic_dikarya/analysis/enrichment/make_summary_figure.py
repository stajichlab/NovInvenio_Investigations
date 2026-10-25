#!/usr/bin/env python3
"""One-slide summary of novelty-candidate enrichment for one study (enrich_novelties_v2.py output).

Adapted from studies/fungi/neocallimastigomycota_vs_early_fungi/analysis/enrichment/make_summary_figure.py.
Panels, stacked in rows (one full-width panel per row; the figure height grows with the rows):
  A  candidates per IN species, split by annotation (any GO/Pfam/InterPro; legend text --annot-label)
  B  top Pfam families in the whole candidate set, family-level test (PRIMARY result)
  C  top GO terms, whole set, family-level test (SECONDARY; IEA-only terms marked)

Panel A rows are labelled with the config Species. When two or more rows would get the
same label, those rows get ' (<Short>)' appended; unique labels stay unchanged. Panel C's
heading reads 'GO: not available (Pfam-only annotation)' when no unit of the study has a
GO term (GO test M_units_annotated == 0 or no GO test row).

Figure width is 14 in. At the saved size, labels and tick labels are 14 pt, panel titles
17 pt, the figure title 18 pt and the footnote 14 pt. Long term names wrap to 2 lines.
A panel with no term to show prints a message in place of an empty axis.

  make_summary_figure.py --enrich OUT_DIR --config CONFIG_CSV --title TEXT --out PREFIX
Run with the NII pixi env python.
"""
import argparse
import csv
import gzip
import math
import re
import sys
import textwrap
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
ap.add_argument("--annot-label", default="UniProt GO/Pfam/InterPro",
                help="panel A legend text for annotated candidates (default: UniProt GO/Pfam/InterPro). "
                     "Use e.g. 'GO/Pfam/InterPro (UniProt) or Pfam (hmmsearch)' when some species have "
                     "Pfam-only HMM annotation")
a = ap.parse_args()

# ---- sizes (points at the saved 14 in width; inches for layout) ----
FIG_W = 14.0
FS = 14          # labels, tick labels, bar annotations, legends, footnote
FS_PANEL = 17    # panel titles (bold)
FS_TITLE = 18    # figure title
LINE = 1.3 / 72  # inches per point of font size for one text line
LABEL_WRAP = 50  # characters per line for term names (max 2 lines)
plt.rcParams.update({"font.size": FS, "axes.labelsize": FS, "xtick.labelsize": FS,
                     "ytick.labelsize": FS, "legend.fontsize": FS})

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

# Heuristic label only (Pfam short name and description). Not a test. TE families and
# generic repeat families are separate classes; the rules are in lib/te_flags.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "lib"))
from te_flags import te_class  # noqa: E402

KIND_LABEL = {"TE": "transposable element (name heuristic)", "repeat": "repeat family (name heuristic)",
              "other": "other"}


def kind(r):
    return KIND_LABEL[te_class(r["name"])]


COL = {"other": "#2a78d6", "transposable element (name heuristic)": "#e08a1e",
       "repeat family (name heuristic)": "#9a9a9a"}
COMPUTATIONAL = {"IEA", "IBA", "ISS", "ISO", "ISA", "ISM", "IGC", "RCA", "NA"}


def go_style(r):
    codes = {x.split(":")[0] for x in r.get("cand_evidence_codes", "").split(";") if x}
    if r["go_evidence"] == "IEA_only":
        return ("#b8b8b8", "///")
    if codes <= COMPUTATIONAL:
        return ("#c7b3e0", None)
    return ("#7b5aa6", None)


def wrap(s, width=LABEL_WRAP, max_lines=2):
    return textwrap.wrap(s, width, max_lines=max_lines, placeholder="…", break_on_hyphens=False) or [s]


def pfam_label(r):
    return "\n".join(wrap(r["name"].replace(": ", " – ")))


def go_label(r):
    if ":" not in r["name"]:
        return r["term"]
    lines = wrap(r["name"].split(":", 1)[1])
    lines[-1] += f" [{r['name'][0].upper()}]"
    return "\n".join(lines)


def not_available(test):
    """True when no unit has a term of the vocabulary (e.g. GO with Pfam-only HMM annotation)."""
    return test is None or int(test["M_units_annotated"]) == 0


def stats_line(test):
    if not_available(test):
        return "not available: no protein in this study has a term of this vocabulary; no test run"
    return (f"{test['n_terms_tested']} terms tested; {test['n_fdr_lt_0.05']} at FDR<0.05; "
            f"annotated units: {test['N_cand_annotated']} cand. / {test['M_units_annotated']} all")


# ---- text blocks and row heights (inches) ----
title_lines = textwrap.wrap(a.title, 96, break_on_hyphens=False) or [""]
foot_lines = textwrap.wrap(a.footnote, 124, break_on_hyphens=False) if a.footnote else []
pf_labels = [pfam_label(r) for r in pf]
go_labels = [go_label(r) for r in go]


def header_h(n_title, n_sub):
    return n_title * FS_PANEL * LINE + n_sub * FS * LINE + 0.12


def bars_h(labels):
    if not labels:
        return 1.1
    nl = max(s.count("\n") + 1 for s in labels)
    pitch = 0.40 if nl == 1 else 0.64
    return len(labels) * pitch + 0.75 + 0.5  # rows + ticks/xlabel + legend row


A_TITLE = f"A  {total:,} candidates in {len(IN)} {a.group} species"
B_TITLE = textwrap.wrap("B  Pfam families enriched in candidates (primary; family units on both sides)", 95)
GO_TEST = tests.get(("whole_set", "family", "go"))
C_TITLE = textwrap.wrap("C  GO: not available (Pfam-only annotation)" if not_available(GO_TEST)
                        else "C  GO terms (secondary; UniProt computational evidence)", 95)
B_SUB = textwrap.wrap(stats_line(tests.get(("whole_set", "family", "pfam"))), 120)
C_SUB = textwrap.wrap(stats_line(tests.get(("whole_set", "family", "go"))), 120)

blocks = [  # (kind, height_in)
    ("title", len(title_lines) * FS_TITLE * LINE + 0.15),
    ("A_head", header_h(1, 0)),
    ("A", len(IN) * (0.40 if len(IN) <= 12 else 0.34) + 0.75 + 0.5),
    ("B_head", header_h(len(B_TITLE), len(B_SUB)) + 0.15),
    ("B", bars_h(pf_labels)),
    ("C_head", header_h(len(C_TITLE), len(C_SUB)) + 0.15),
    ("C", bars_h(go_labels)),
]
if foot_lines:
    blocks.append(("foot", len(foot_lines) * FS * LINE + 0.25))
heights = [h for _, h in blocks]

fig = plt.figure(figsize=(FIG_W, sum(heights)), layout="constrained")
fig.get_layout_engine().set(w_pad=0.04, h_pad=0.02, hspace=0, wspace=0)
subs = dict(zip([k for k, _ in blocks], fig.subfigures(len(blocks), 1, height_ratios=heights, hspace=0)))
block_h = dict(blocks)


def text_block(name, lines_spec):
    """lines_spec: list of (text, fontsize, weight, color). Drawn top-down at the left edge of
    the text-only row `name`, whose height (inches) was sized for these lines."""
    sf, h_in = subs[name], block_h[name]
    y = 1.0 - 0.05 / h_in
    for text, fs, weight, color in lines_spec:
        sf.text(0.005, y, text, ha="left", va="top", fontsize=fs, fontweight=weight, color=color,
                transform=sf.transSubfigure)
        y -= fs * LINE / h_in


text_block("title", [(ln, FS_TITLE, "normal", "black") for ln in title_lines])
text_block("A_head", [(A_TITLE, FS_PANEL, "bold", "black")])
text_block("B_head", [(ln, FS_PANEL, "bold", "black") for ln in B_TITLE]
           + [(ln, FS, "normal", "#444") for ln in B_SUB])
text_block("C_head", [(ln, FS_PANEL, "bold", "black") for ln in C_TITLE]
           + [(ln, FS, "normal", "#444") for ln in C_SUB])
if foot_lines:
    text_block("foot", [(ln, FS, "normal", "#444") for ln in foot_lines])


def despine(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


# --- A ---
sf = subs["A"]
ax = sf.subplots()
ys = list(range(len(IN)))[::-1]
xmax = max([int(cov[(s, "candidates")]["n_proteins"]) for s in IN] + [1])
for y, s in zip(ys, IN):
    c = cov[(s, "candidates")]
    n, na = int(c["n_proteins"]), int(c["n_any"])
    ax.barh(y, na, color="#2a78d6", height=0.65)
    ax.barh(y, n - na, left=na, color="#c9c9c9", height=0.65)
    u = cov[(s, "universe")]
    frac = f" ({na / n:.0%})" if n else ""
    ax.text(n + xmax * 0.015, y, f"{n:,}{frac}; proteome {float(u['frac_any']):.0%}",
            va="center", fontsize=FS)
ax.set_yticks(ys)
_label_count = {}
for s in IN:
    _label_count[LABEL[s]] = _label_count.get(LABEL[s], 0) + 1
A_LABELS = [f"{LABEL[s]} ({s})" if _label_count[LABEL[s]] > 1 else f"{LABEL[s]}" for s in IN]
ax.set_yticklabels(A_LABELS, style="italic")
ax.set_xlabel(a.candidate_label)
ax.set_xlim(0, xmax * 1.55)
ax.set_ylim(-0.6, len(IN) - 0.4)
despine(ax)
annot_label = (f"{a.annot_label} ({n_annot:,}, {n_annot / total:.0%})" if total
               else f"{a.annot_label} (0)")
sf.legend(handles=[Patch(color="#2a78d6", label=annot_label), Patch(color="#c9c9c9", label="no annotation")],
          loc="outside lower left", ncol=2, frameon=False)


def barpanel(sf, rows, labels, colorf, handles, test, what, whats):
    ax = sf.subplots()
    if not rows:
        n_fdr = 0 if not_available(test) else int(test["n_fdr_lt_0.05"])
        if not_available(test):
            msg = f"{whats[0].upper() + whats[1:]} not available: no {what} annotation in this study."
        elif n_fdr == 0:
            msg = f"No {what} passed FDR<0.05."
        else:
            msg = f"{n_fdr} {whats} at FDR<0.05, none with k>={a.min_k} candidate units."
        ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=FS + 2, color="#444",
                transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.6", facecolor="#f2f2f2", edgecolor="#bbbbbb"))
        ax.set_axis_off()
        return
    ys = list(range(len(rows)))[::-1]
    vals = [-math.log10(max(float(r["fdr_bh"]), 1e-300)) for r in rows]
    vmax = max(vals)
    for y, r, v in zip(ys, rows, vals):
        col, hatch = colorf(r)
        ax.barh(y, v, color=col, hatch=hatch, edgecolor="white" if not hatch else "#555", linewidth=0.5,
                height=0.7)
        ax.text(v + vmax * 0.015, y, f'{float(r["fold"]):.1f}x  k={r["k_cand_with_term"]}/{r["n_bg_with_term"]}',
                va="center", fontsize=FS)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, vmax * 1.4)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("-log10 FDR (BH within vocabulary)")
    despine(ax)
    sf.legend(handles=handles, loc="outside lower left", ncol=len(handles), frameon=False)


barpanel(subs["B"], pf, pf_labels, lambda r: (COL[kind(r)], None),
         [Patch(color=c, label=k_) for k_, c in COL.items()],
         tests.get(("whole_set", "family", "pfam")), "Pfam family", "Pfam families")
barpanel(subs["C"], go, go_labels, go_style,
         [Patch(facecolor="#b8b8b8", hatch="///", edgecolor="#555", label="IEA only"),
          Patch(color="#c7b3e0", label="IEA + IBA/ISS (computational)"),
          Patch(color="#7b5aa6", label="incl. curated/experimental code")],
         tests.get(("whole_set", "family", "go")), "GO term", "GO terms")

a.out.parent.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=200, bbox_inches="tight", pad_inches=0.15)
print(f"Pfam rows shown {len(pf)} of {len(pf_all)} (FDR<0.05, k>={a.min_k}); GO shown {len(go)} of {len(go_all)}")
