#!/usr/bin/env python3
"""Locus summary figures and tables from island_loci.json (frequency-stratified set).

Works for any pangenome study. One focal species is shown (--species); strains of other
species in the matrix are left out of every count and figure.

Outputs in --outdir:
  locus_table.tsv          one row per locus: frequency, haplotype states, locus type, Starship overlap
  locus_overview.png       frequency vs number of common states, and locus-type counts by frequency bin
  locus_strips/<key>.png   one figure per locus: strains (by Mash clade) x genes (flank | locus | flank)
  locus_strips_selected.png montage of the loci in --pick plus one locus per type and frequency bin
  starship_overlay.png     Af293 chromosomes with Starship spans and locus positions (needs --starships)
  idiomorph_candidates.tsv loci where the empty state has other genes between the flanks

Cell codes (lib/island_locus.py): 1 and 3 = gene in place (annotated or rescued); 6 = DNA at the
site but no annotated gene; 7 = DNA checked and absent; 0 = not DNA-checked; 2, 4, 5 = elsewhere,
rescued elsewhere, contig break.

Locus type, on the focal-species strains only:
  relocated          >= half of the strains have >= half of the locus genes elsewhere in the genome
                     (content is present but not at this site; a mobile-element signature)
  unresolved         >= half of the strains have >= half of the locus genes as contig break, elsewhere
                     or not DNA-checked (often a locus at a contig end); the pattern is not interpretable
  absent_in_focal    one state covers >= 90% of strains and it has no content (content is in other species)
  monomorphic        one state covers >= 90% of strains, with content
  rare_content       content (>= half of the locus columns present) in < 10% of strains, >= 3 strains
  swap_candidate two common states, one with content and one empty, and the empty state has
                     >= 1 other gene between the flanks (swap candidate; idiomorph only if checked)
  presence_absence   same two states, but no other gene between the flanks (insertion or deletion)
  content_variants   two or more common states that differ in which genes are present
"""
import argparse
import collections
import csv
import json
import os
import statistics
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

PRESENT = set("136")
ANNOT = set("13")
COMMON_FRAC = 0.10
CLASS_COL = {"nlr": "#c53030", "transporter": "#d69e2e", "secondary_metabolite": "#38a169",
             "other": "#718096", "unannotated": "#cbd5e0"}
TYPE_COL = {"absent_in_focal": "#f6e05e", "unresolved": "#e2e8f0", "relocated": "#2c7a7b", "monomorphic": "#a0aec0", "rare_content": "#805ad5", "swap_candidate": "#e53e3e",
            "presence_absence": "#3182ce", "content_variants": "#dd6b20"}
# strip cell colours: annotated, DNA only, DNA absent, not checked, other
CELL = {"1": 0, "3": 0, "6": 1, "7": 2, "0": 3, "2": 4, "4": 4, "5": 4}
CELL_COL = ["#2b6cb0", "#68d391", "#f7fafc", "#cbd5e0", "#ed8936"]
CELL_LAB = ["gene in place", "DNA at site, no gene", "DNA checked, absent", "not DNA-checked",
            "elsewhere / contig break"]
CLADE_COLS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999"]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def zstd_lines(path):
    if path.endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        yield from p.stdout
        p.wait()
    else:
        with open(path) as fh:
            yield from fh


def gap_between_flanks(rep, nl, nb):
    """Genes between the innermost flank columns that are not locus columns, or None."""
    d = rep.get("d") or []
    if len(d) < nl + nb + 1 or nl < 1:
        return None
    a, b = d[nl - 1], d[nl + nb]
    if not a or not b or len(a) < 2 or len(b) < 2 or a[0] != b[0]:
        return None
    lo, hi = sorted((a[1], b[1]))
    inside = sum(1 for x in d[nl:nl + nb] if x and len(x) > 2 and x[0] == a[0] and lo < x[1] < hi)
    return max(hi - lo - 1 - inside, 0)


def pattern(codes, nl, nb):
    tr = {"1": "P", "3": "P", "6": "P", "7": "A"}
    return "".join(tr.get(ch, "?") for ch in codes[nl:nl + nb])


def kind(pat):
    n = len(pat)
    if pat.count("P") >= n / 2:
        return "content"
    if pat.count("A") >= n / 2:
        return "empty"
    return "mixed"


def classify(L, focal):
    nl, nb = L["n_left"], L["n_locus"]
    states = collections.defaultdict(lambda: {"n": 0, "gaps": []})
    total = 0
    moved = 0
    unknown = 0
    for r in L["rows"]:
        strains = [s for s in r["strains"] if s in focal]
        if not strains:
            continue
        p = pattern(r["codes"], nl, nb)
        st = states[p]
        st["n"] += len(strains)
        if sum(ch in "24" for ch in r["codes"][nl:nl + nb]) >= nb / 2:
            moved += len(strains)
        if sum(ch in "0245" for ch in r["codes"][nl:nl + nb]) >= nb / 2:
            unknown += len(strains)
        g = gap_between_flanks(r.get("rep") or {}, nl, nb)
        if g is not None:
            st["gaps"].append((g, len(strains)))
        total += len(strains)
    if not total:
        return None
    ranked = sorted(states.items(), key=lambda kv: -kv[1]["n"])
    content_n = sum(v["n"] for p, v in states.items() if kind(p) == "content")
    out = {"n_strains": total, "content_n": content_n, "content_frac": content_n / total,
           "n_states": len(ranked), "top_share": ranked[0][1]["n"] / total,
           "n_common": sum(1 for _, v in ranked if v["n"] >= max(3, COMMON_FRAC * total)),
           "top_state": ranked[0][0], "pair": "", "elsewhere_frac": moved / total,
           "unknown_frac": unknown / total}
    common = [(p, v) for p, v in ranked if v["n"] >= max(3, COMMON_FRAC * total)]
    ctype = "content_variants"
    if out["elsewhere_frac"] >= 0.5:
        ctype = "relocated"
    elif out["unknown_frac"] >= 0.5:
        ctype = "unresolved"
    elif out["top_share"] >= 0.90:
        ctype = "monomorphic" if out["content_frac"] >= COMMON_FRAC else "absent_in_focal"
    elif out["content_frac"] < COMMON_FRAC and content_n >= 3:
        ctype = "rare_content"
    else:
        # pairs of common states: the lacking state has DNA-absent (A) where the other has a gene (P)
        best = None
        for x in range(len(common)):
            for y in range(x + 1, len(common)):
                px, vx = common[x]
                py, vy = common[y]
                lack_y = sum(1 for a_, b_ in zip(px, py) if a_ == "P" and b_ == "A")
                lack_x = sum(1 for a_, b_ in zip(px, py) if a_ == "A" and b_ == "P")
                if lack_y == 0 and lack_x == 0:
                    continue
                lacking = (py, vy, px, vx) if lack_y >= lack_x else (px, vx, py, vy)
                gaps = [g for g, n in lacking[1]["gaps"] for _ in range(n)]
                med = statistics.median(gaps) if gaps else None
                score = (med is not None and med >= 1, vx["n"] + vy["n"])
                if best is None or score > best[0]:
                    best = (score, lacking, med)
        if best:
            _, (pl, vl, ph, vh), med = best
            out["empty_gap"] = med
            out["pair"] = f"{ph}({vh['n']}) vs {pl}({vl['n']})"
            ctype = "swap_candidate" if med is not None and med >= 1 else "presence_absence"
    out["type"] = ctype
    return out


def af293_positions(args):
    """family -> (contig, rank) and rank -> (start, end) for the reference strain."""
    ref = args.ref_strain
    fam = {}
    for ln in zstd_lines(args.family_positions):
        f = ln.rstrip("\n").split("\t")
        if f[0] == ref:
            fam[f[1]] = (f[2], int(f[3]))
    genes = []
    for ln in zstd_lines(args.gene_positions):
        f = ln.rstrip("\n").split("\t")
        if f[0] == ref:
            genes.append((f[2], int(f[3]), int(f[4])))
    for r in read_tsv(args.rescue_positions):
        if r["Short"] == ref:
            genes.append((r["contig"], int(r["start"]), int(r["start"])))
    genes.sort()
    by_contig = collections.defaultdict(list)
    for c, s, e in genes:
        by_contig[c].append((s, e))
    # ranks run on across contigs (sorted contig order), so keep each contig's offset
    off, run = {}, 0
    for c in sorted(by_contig):
        off[c] = run
        run += len(by_contig[c])
    by_contig["_offset"] = off
    return fam, by_contig


SCATTERED = []
BAD = []  # loci whose ranks fall outside the reference strain's gene list


def starship_overlap(L, fam, by_contig, ships):
    nl, nb = L["n_left"], L["n_locus"]
    fams = L["families"]
    block = [fam[f] for f in fams[nl:nl + nb] if f in fam]
    flank = [fam[f] for f in (fams[nl - 1:nl] + fams[nl + nb:nl + nb + 1]) if f in fam]
    pts = block or flank
    if not pts:
        return None
    contig = collections.Counter(c for c, _ in pts).most_common(1)[0][0]
    g = by_contig.get(contig)
    if not g:
        return None
    off = by_contig["_offset"][contig]
    ranks = [r - off for c, r in pts if c == contig]
    lo, hi = min(ranks), max(ranks)
    if hi - lo > 3 * (nb + 2):
        SCATTERED.append(L["key"])  # columns spread over many genes: paralogs, no compact site
        return None
    if lo < 0 or hi >= len(g):
        BAD.append(L["key"])
        return None
    span = (g[lo][0], g[hi][1])
    hit = [s for s in ships if s["contig"] == contig and s["start"] < span[1] and s["end"] > span[0]]
    return {"contig": contig, "start": span[0], "end": span[1], "from_block": bool(block), "ships": hit}


def draw_strip(ax_main, ax_clade, L, rows_by_strain, order, clade_of, clade_idx, with_labels=True):
    nl, nb, nr = L["n_left"], L["n_locus"], L["n_right"]
    ncol = nl + nb + nr
    M = np.full((len(order), ncol), 3, dtype=int)
    for i, s in enumerate(order):
        codes = rows_by_strain[s]
        for j in range(min(ncol, len(codes))):
            M[i, j] = CELL.get(codes[j], 4)
    ax_main.imshow(M, aspect="auto", cmap=ListedColormap(CELL_COL), vmin=0, vmax=4,
                   interpolation="nearest")
    ax_main.axvline(nl - .5, color="k", lw=.8)
    ax_main.axvline(nl + nb - .5, color="k", lw=.8)
    ax_main.set_yticks([])
    ax_main.set_xticks(range(ncol))
    labs = []
    for j in range(ncol):
        dom = (L.get("family_domains") or [""] * ncol)[j] if j < len(L.get("family_domains") or []) else ""
        labs.append(dom.split(",")[0][:12] if dom and with_labels else str(j + 1))
    if ncol > 30:
        labs = [""] * ncol
    ax_main.set_xticklabels(labs, rotation=90, fontsize=5)
    cm = np.array([[clade_idx[clade_of.get(s, "NA")]] for s in order])
    ax_clade.imshow(cm, aspect="auto", cmap=ListedColormap(CLADE_COLS), vmin=0, vmax=len(CLADE_COLS) - 1,
                    interpolation="nearest")
    ax_clade.set_xticks([])
    ax_clade.set_yticks([])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--clades", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--species", required=True, help="focal species; others are dropped")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--pick", default="", help="comma list of locus keys to include in the montage")
    ap.add_argument("--starships", default="", help="Af293 Starship table (Table S7 TSV)")
    ap.add_argument("--labels", default="", help="TSV with columns key, label (known loci, e.g. MAT); keys are per run")
    ap.add_argument("--ref_strain", default="Asfu_Af293")
    ap.add_argument("--family_positions", default="")
    ap.add_argument("--gene_positions", default="")
    ap.add_argument("--rescue_positions", default="")
    a = ap.parse_args()

    d = json.load(open(a.loci))
    loci = d["loci"]
    cfg = list(csv.DictReader(open(a.config)))
    focal = {r["Short"] for r in cfg if r["Species"] == a.species}
    clade_of = {r["Short"]: r["mash_clade"] for r in read_tsv(a.clades)}
    clade_names = sorted({clade_of.get(s, "NA") for s in focal})
    clade_idx = {c: i for i, c in enumerate(clade_names)}
    os.makedirs(os.path.join(a.outdir, "locus_strips"), exist_ok=True)

    known = {r["key"]: r["label"] for r in read_tsv(a.labels)} if a.labels else {}
    ships, fam, by_contig = [], {}, {}
    if a.starships:
        for r in read_tsv(a.starships):
            ships.append({"id": r["starship_id"], "name": r["name"], "hap": r["haplotype"],
                          "contig": r["contig"], "start": int(r["start"]), "end": int(r["end"])})
        fam, by_contig = af293_positions(a)

    table = []
    for L in loci:
        c = classify(L, focal)
        if not c:
            continue
        c.update(key=L["key"], tier=L["tier"], size=L["size"], dom_class=L["dominant_class"], exemplar=L["exemplar"],
                 n_species_called=len(L["counts_by_species"]))
        ov = starship_overlap(L, fam, by_contig, ships) if ships else None
        c["starship"] = ";".join(f"{s['name']}-{s['hap']}" for s in ov["ships"]) if ov and ov["ships"] else ""
        c["af293_pos"] = f"{ov['contig']}:{ov['start']}-{ov['end']}" if ov else ""
        c["_ov"] = ov
        c["label"] = known.get(L["key"], "")
        table.append(c)
    cols = ["key", "label", "size", "dom_class", "n_strains", "content_n", "content_frac", "n_common", "n_states",
            "top_share", "elsewhere_frac", "unknown_frac", "tier", "type", "empty_gap", "pair", "starship", "af293_pos", "exemplar"]
    with open(os.path.join(a.outdir, "locus_table.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in table:
            fh.write("\t".join("" if c.get(k) is None else (f"{c[k]:.3f}" if isinstance(c.get(k), float)
                                                              else str(c[k])) for k in cols) + "\n")
    idio = [c for c in table if c["type"] == "swap_candidate"]
    with open(os.path.join(a.outdir, "idiomorph_candidates.tsv"), "w") as fh:
        fh.write("key\tpair\tempty_state_median_extra_genes\tsize\tdom_class\tstarship\n")
        for c in idio:
            fh.write(f"{c['key']}\t{c['pair']}\t{c.get('empty_gap')}\t{c['size']}\t{c['dom_class']}\t{c['starship']}\n")

    # ---- overview figure
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    for t, col in TYPE_COL.items():
        s = [c for c in table if c["type"] == t]
        ax[0].scatter([max(c["content_frac"], 0.002) for c in s], [c["n_common"] for c in s],
                      s=[18 + 8 * min(c["size"], 40) for c in s], c=col, alpha=.7, edgecolor="k",
                      linewidth=.3, label=f"{t} ({len(s)})")
    ss = [c for c in table if c["starship"]]
    ax[0].scatter([max(c["content_frac"], 0.002) for c in ss], [c["n_common"] for c in ss],
                  s=140, facecolors="none", edgecolors="k", linewidth=1.4, label=f"overlaps Af293 Starship ({len(ss)})")
    ax[0].set_xscale("log")
    ax[0].set_xlabel(f"fraction of {a.species} strains with the locus content (log)")
    ax[0].set_ylabel("states held by >=10% of strains")
    ax[0].set_title(f"{a.label or a.species}: {len(table)} loci")
    ax[0].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    edges = [0, 0.03, 0.10, 0.30, 0.70, 1.0001]
    bn = [f"{edges[i]:.0%}-{min(edges[i + 1], 1):.0%}" for i in range(len(edges) - 1)]
    bot = np.zeros(len(bn))
    for t, col in TYPE_COL.items():
        cnt = np.array([sum(1 for c in table if c["type"] == t and edges[i] <= c["content_frac"] < edges[i + 1])
                        for i in range(len(bn))])
        ax[1].bar(bn, cnt, bottom=bot, color=col, label=t)
        bot += cnt
    ax[1].set_xlabel("locus content frequency bin")
    ax[1].set_ylabel("loci")
    ax[1].set_title("Locus type by frequency")
    ax[1].legend(fontsize=7)
    cls = collections.Counter(c["dom_class"] for c in table)
    ax[2].bar(list(cls), [cls[k] for k in cls], color=[CLASS_COL.get(k, "#718096") for k in cls])
    ax[2].set_title("Dominant gene class")
    ax[2].tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "locus_overview.png"), dpi=130)
    plt.close(fig)

    # ---- per-locus strips
    byk = {L["key"]: L for L in loci}
    tab = {c["key"]: c for c in table}
    pdf = PdfPages(os.path.join(a.outdir, "locus_strips.pdf"))

    def one(L, c, order_cache={}):
        by_strain = {}
        for r in L["rows"]:
            for s in r["strains"]:
                if s in focal:
                    by_strain[s] = r["codes"]
        order = sorted(by_strain, key=lambda s: (clade_idx[clade_of.get(s, "NA")], pattern(by_strain[s], L["n_left"], L["n_locus"]), s))
        fig = plt.figure(figsize=(7.5, 7))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 30], wspace=.01)
        axc, axm = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        draw_strip(axm, axc, L, by_strain, order, clade_of, clade_idx)
        ttl = (f"{L['key']}  {L['size']} genes  content {c['content_frac']:.0%} of {c['n_strains']} strains  "
               f"type: {c['type']}" + (f"  Known: {c['label']}" if c["label"] else "") + (f"  Starship: {c['starship']}" if c["starship"] else ""))
        fig.suptitle(ttl, fontsize=8, x=0.02, ha="left")
        fig.legend([Patch(color=x) for x in CELL_COL] + [Patch(color=CLADE_COLS[i]) for i in range(len(clade_names))],
                   CELL_LAB + clade_names, ncol=5, fontsize=6, loc="lower center")
        fig.subplots_adjust(bottom=0.16, top=0.93)
        return fig

    for L in loci:
        c = tab.get(L["key"])
        if not c:
            continue
        fig = one(L, c)
        fig.savefig(os.path.join(a.outdir, "locus_strips", f"{L['key']}.png"), dpi=110)
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()

    # ---- montage of selected loci
    pick = [k for k in a.pick.split(",") if k in tab]
    for t in TYPE_COL:
        for lo, hi in ((0, 0.10), (0.10, 0.70), (0.70, 1.01)):
            cand = [c for c in table if c["type"] == t and lo <= c["content_frac"] < hi and c["key"] not in pick]
            cand.sort(key=lambda c: -c["size"])
            if cand:
                pick.append(cand[0]["key"])
    pick = pick[:12]
    n = len(pick)
    cols_n = 4
    rows_n = (n + cols_n - 1) // cols_n
    fig = plt.figure(figsize=(4.2 * cols_n, 4.2 * rows_n + 0.8))
    outer = fig.add_gridspec(rows_n, cols_n, hspace=.5, wspace=.15)
    for i, k in enumerate(pick):
        L, c = byk[k], tab[k]
        by_strain = {s: r["codes"] for r in L["rows"] for s in r["strains"] if s in focal}
        order = sorted(by_strain, key=lambda s: (clade_idx[clade_of.get(s, "NA")], pattern(by_strain[s], L["n_left"], L["n_locus"]), s))
        inner = outer[i // cols_n, i % cols_n].subgridspec(1, 2, width_ratios=[1, 25], wspace=.01)
        axc, axm = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        draw_strip(axm, axc, L, by_strain, order, clade_of, clade_idx, with_labels=False)
        axm.set_title(f"{k} {c['type']}" + (f" [{c['label']}]" if c["label"] else "") + f"\n{c['content_frac']:.0%} content" + (f", {c['starship']}" if c["starship"] else ""),
                      fontsize=7)
    fig.legend([Patch(color=x) for x in CELL_COL] + [Patch(color=CLADE_COLS[i]) for i in range(len(clade_names))],
               CELL_LAB + clade_names, ncol=6, fontsize=7, loc="lower center")
    fig.savefig(os.path.join(a.outdir, "locus_strips_selected.png"), dpi=110, bbox_inches="tight")
    plt.close(fig)

    # ---- Starship overlay
    if ships:
        contigs = sorted({s["contig"] for s in ships} | {c["_ov"]["contig"] for c in table if c["_ov"]})
        fig, axes = plt.subplots(len(contigs), 1, figsize=(12, 1.0 + 0.8 * len(contigs)), sharex=False)
        for ax_, ctg in zip(np.atleast_1d(axes), contigs):
            ln = max([g[1] for g in by_contig.get(ctg, [(0, 0)])] if ctg != '_offset' else [0])
            ax_.set_xlim(0, ln)
            ax_.set_ylim(0, 2)
            ax_.set_yticks([])
            ax_.set_ylabel(ctg, rotation=0, ha="right", va="center", fontsize=7)
            ax_.axhline(1, color="#cbd5e0", lw=.6)
            for s in ships:
                if s["contig"] == ctg:
                    ax_.add_patch(plt.Rectangle((s["start"], 0.15), s["end"] - s["start"], 0.7, color="#1a202c", alpha=.8))
                    ax_.text((s["start"] + s["end"]) / 2, 0.5, f"{s['name']}-{s['hap']}", color="w", ha="center",
                             va="center", fontsize=5)
            for c in table:
                ov = c["_ov"]
                if ov and ov["contig"] == ctg:
                    col = TYPE_COL[c["type"]]
                    ax_.plot([(ov["start"] + ov["end"]) / 2], [1.5], marker="v", color=col, ms=5 if not c["starship"] else 8,
                             markeredgecolor="k" if c["starship"] else col, markeredgewidth=1)
                    if c["starship"]:
                        ax_.text((ov["start"] + ov["end"]) / 2, 1.85, c["key"], fontsize=6, ha="center")
        fig.legend([plt.Line2D([], [], marker="v", ls="", color=col) for col in TYPE_COL.values()] +
                   [Patch(color="#1a202c")], list(TYPE_COL) + ["Af293 Starship (Table S7)"], ncol=6, fontsize=7,
                   loc="lower center")
        fig.suptitle(f"{a.ref_strain}: loci (triangles) and Starship spans (bars)", fontsize=9)
        fig.subplots_adjust(bottom=0.1, top=0.94, hspace=.6)
        fig.savefig(os.path.join(a.outdir, "starship_overlay.png"), dpi=120)
        plt.close(fig)

    if BAD:
        print(f"rank outside gene list (no Starship test): {len(BAD)} loci {BAD[:10]}", file=sys.stderr)
    print(f"no compact Af293 site (scattered columns): {len(SCATTERED)} loci", file=sys.stderr)
    tc = collections.Counter(c["type"] for c in table)
    print(f"{len(table)} loci in focal species; types {dict(tc)}; Starship overlaps "
          f"{sum(1 for c in table if c['starship'])}; swap candidates {len(idio)}", file=sys.stderr)


if __name__ == "__main__":
    main()
