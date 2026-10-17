#!/usr/bin/env python3
"""Shared-domain context, gene level (tasks 1, 2i, 2iii). Reads the stage 1 tables.

Writes (tables/):
  genes_per_genome_summary.tsv     per genus/species x label: proteomes with >= 1, median, min, max, total
  gene_family_class_split.tsv      per run x label: hit proteins and distinct families by family bin
  same_protein_combinations.tsv    per genus: label combinations on one protein (n proteins, n genomes)
  reference_window_pairs.tsv       Af293 / RS / UCSF_Cp_Silveira: gene pairs of two labels within
                                   5/10/20 genes and 50 kb on one contig, vs a rotation null
  allstrain_neighbour_rates.tsv    all strains: fraction of label-A genes with a label-B gene
                                   within 5 genes (descriptive, no null)

Rotation null: for each label and contig, the label's gene ranks are rotated (circular shift)
by one random offset; labels move independently; 1,000 permutations, seed 20261017.
A gene pair = two different genes; a gene with both labels does not pair with itself.
"""
from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from pathlib import Path
from statistics import median

import numpy as np

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
D = NII / "studies/fungi/coccidioides_pangenome/analysis/shared_domains_context"
T = D / "tables"
TIER1 = ["APH", "DUF3435", "DUF3723", "DUF3723_N"]
LABELS = TIER1 + ["PF11001", "HrmA", "NPHP3_N", "Condensation", "PP-binding", "AMP-binding", "SPOK2_N"]
BINS = ["core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "outgroup_only", "NA"]
WINDOWS = [5, 10, 20]
BP = 50000
NPERM = 1000
SEED = 20261017


def read_tsv(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def add_bh(rows, keys, pcol, qcol):
    """Benjamini-Hochberg within groups defined by keys, over rows with a numeric p."""
    groups = defaultdict(list)
    for r in rows:
        if isinstance(r.get(pcol), float):
            groups[tuple(r[k] for k in keys)].append(r)
    for rs in groups.values():
        rs.sort(key=lambda r: r[pcol])
        m = len(rs)
        q = 1.0
        for i in range(m - 1, -1, -1):
            q = min(q, rs[i][pcol] * m / (i + 1))
            rs[i][qcol] = round(q, 5)


def write_tsv(path, rows, header=None):
    if header is None:
        header = []
        for r in rows:
            for k in r:
                if k not in header:
                    header.append(k)
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def proteins(rows):
    """(short, pid) -> dict(labels=set, row=first row)"""
    out = {}
    for r in rows:
        k = (r["short"], r["protein_id"])
        if k not in out:
            out[k] = dict(labels=set(r["protein_labels"].split(",")), row=r)
    return out


def main():
    g_af = read_tsv(T / "genes.afum.tsv.gz")
    g_cp = read_tsv(T / "genes.cocci.tsv.gz")
    gc = read_tsv(T / "genome_counts.tsv.gz")
    P = {"afum": proteins(g_af), "cocci": proteins(g_cp)}

    # ---- task 1a: genes per genome ----
    rows = []
    groups = {
        "A. fumigatus (293)": lambda r: r["set"] == "afum" and r["species"] == "Aspergillus fumigatus",
        "A. fumigatus outgroups (A. fischeri, A. lentulus)": lambda r: r["set"] == "afum" and r["species"] != "Aspergillus fumigatus",
        "Coccidioides all (529)": lambda r: r["set"] == "cocci",
        "C. immitis (169)": lambda r: r["set"] == "cocci" and r["species"] == "Coccidioides immitis",
        "C. posadasii (360)": lambda r: r["set"] == "cocci" and r["species"] == "Coccidioides posadasii",
    }
    for gname, f in groups.items():
        for lab in LABELS:
            c = [int(r["n_proteins"]) for r in gc if f(r) and r["label"] == lab]
            rows.append(dict(group=gname, label=lab, n_proteomes=len(c), n_with_ge1=sum(x > 0 for x in c),
                             median=median(c), min=min(c), max=max(c), total=sum(c)))
    write_tsv(T / "genes_per_genome_summary.tsv", rows)

    # ---- task 1b: family class split ----
    rows = []
    for s, tags in (("afum", ["af"]), ("cocci", ["v5", "v2"])):
        for t in tags:
            for lab in LABELS:
                bins_p, fams = Counter(), {}
                n_not = 0
                for (short, pid), d in P[s].items():
                    if lab not in d["labels"]:
                        continue
                    r = d["row"]
                    fam = r.get(f"{t}_family", "")
                    if fam in ("not_in_run", "", None):
                        n_not += 1
                        continue
                    b = r.get(f"{t}_family_bin") or "NA"
                    bins_p[b] += 1
                    fams[fam] = b
                bins_f = Counter(fams.values())
                isl_f = {r["row"][f"{t}_family"] for r in P[s].values()
                         if lab in r["labels"] and r["row"].get(f"{t}_family_n_islands", "0") not in ("", "0")}
                row = dict(tag=t, label=lab, n_proteins_in_run=sum(bins_p.values()), n_proteins_not_in_run=n_not,
                           n_families=len(fams), n_families_in_ge1_island=len(isl_f))
                for b in BINS:
                    row[f"proteins_{b}"] = bins_p[b]
                for b in BINS:
                    row[f"families_{b}"] = bins_f[b]
                rows.append(row)
    write_tsv(T / "gene_family_class_split.tsv", rows)

    # ---- task 2i: same protein ----
    rows = []
    for s, sel in (("afum", lambda r: r["species"] == "Aspergillus fumigatus"),
                   ("cocci", lambda r: True)):
        comb = Counter()
        gen = defaultdict(set)
        tot = Counter()
        for (short, pid), d in P[s].items():
            if not sel(d["row"]):
                continue
            key = "+".join(l for l in LABELS if l in d["labels"])
            comb[key] += 1
            gen[key].add(short)
            for l in d["labels"]:
                tot[l] += 1
        for key, n in comb.most_common():
            labs = key.split("+")
            rows.append(dict(genus=("Aspergillus fumigatus (293)" if s == "afum" else "Coccidioides (529)"),
                             combination=key, n_labels=len(labs), n_proteins=n, n_genomes=len(gen[key]),
                             frac_of_first_label_proteins=round(n / tot[labs[0]], 4)))
    write_tsv(T / "same_protein_combinations.tsv", rows)

    # ---- task 2iii: reference windows ----
    order = read_tsv(T / "reference_gene_order.tsv.gz")
    rng = np.random.default_rng(SEED)
    out = []
    for ref, s in (("Asfu_Af293", "afum"), ("RS", "cocci"), ("UCSF_Cp_Silveira", "cocci")):
        contigs = defaultdict(list)
        for r in order:
            if r["short"] == ref:
                contigs[r["contig"]].append((int(r["rank"]), int(r["start"]), int(r["end"]), r["protein_id"]))
        cdata = {}
        lab_ranks = defaultdict(dict)   # label -> contig -> np.array of ranks
        for c, lst in contigs.items():
            lst.sort()
            st = np.array([x[1] for x in lst])
            en = np.array([x[2] for x in lst])
            cdata[c] = (st, en, len(lst))
            for lab in LABELS:
                rk = [x[0] for x in lst if (ref, x[3]) in P[s] and lab in P[s][(ref, x[3])]["labels"]]
                if rk:
                    lab_ranks[lab][c] = np.array(rk)
        n_lab_genes = {lab: sum(len(v) for v in lab_ranks[lab].values()) for lab in LABELS}
        n_lab_contigs = {lab: len(lab_ranks[lab]) for lab in LABELS}

        def count(lr, a, b):
            """counts for windows + bp, unordered gene pairs of two different genes"""
            res = np.zeros(len(WINDOWS) + 1, dtype=np.int64)
            for c in set(lr[a]) & set(lr[b]):
                st, en, n = cdata[c]
                ra, rb = np.unique(lr[a][c]), np.unique(lr[b][c])
                I = np.repeat(ra, len(rb))
                J = np.tile(rb, len(ra))
                keep = I != J
                lo, hi = np.minimum(I[keep], J[keep]), np.maximum(I[keep], J[keep])
                if len(lo) == 0:
                    continue
                key = np.unique(lo * n + hi)
                lo, hi = key // n, key % n
                dr = hi - lo
                for w, W in enumerate(WINDOWS):
                    res[w] += int((dr <= W).sum())
                gap = np.maximum(0, np.maximum(st[lo], st[hi]) - np.minimum(en[lo], en[hi]))
                res[-1] += int((gap <= BP).sum())
            return res

        obs = {}
        for a, b in combinations_with_replacement(LABELS, 2):
            obs[(a, b)] = count(lab_ranks, a, b)
        perm = {k: np.zeros((NPERM, len(WINDOWS) + 1), dtype=np.int64) for k in obs if k[0] != k[1]}
        for p in range(NPERM):
            lr = {}
            for lab in LABELS:
                lr[lab] = {c: (v + rng.integers(cdata[c][2])) % cdata[c][2] for c, v in lab_ranks[lab].items()}
            for k in perm:
                perm[k][p] = count(lr, *k)
        for (a, b), o in obs.items():
            for w, wname in enumerate([f"le{W}_genes" for W in WINDOWS] + ["le50kb"]):
                row = dict(reference=ref, label_a=a, label_b=b, window=wname,
                           n_genes_a=n_lab_genes[a], n_genes_b=n_lab_genes[b],
                           n_contigs_a=n_lab_contigs[a], n_contigs_b=n_lab_contigs[b],
                           observed_pairs=int(o[w]))
                if a != b:
                    pv = perm[(a, b)][:, w]
                    row.update(perm_mean=round(float(pv.mean()), 3),
                               perm_q97_5=float(np.percentile(pv, 97.5)),
                               obs_over_exp=(round(o[w] / pv.mean(), 3) if pv.mean() > 0 else "NA"),
                               p_upper=round((1 + int((pv >= o[w]).sum())) / (NPERM + 1), 5))
                else:
                    row.update(perm_mean="same_label(rotation preserves spacing)", p_upper="NA")
                out.append(row)
        print(f"{ref}: contigs={len(contigs)} genes={sum(v[2] for v in cdata.values())} "
              f"label genes={n_lab_genes}", file=sys.stderr)
    add_bh(out, ("reference", "window"), "p_upper", "q_upper_bh")
    write_tsv(T / "reference_window_pairs.tsv", out)

    # ---- all strains, descriptive neighbour rates (+/-5 genes) ----
    rows = []
    for s, sel in (("afum", lambda r: r["species"] == "Aspergillus fumigatus"), ("cocci", lambda r: True)):
        byc = defaultdict(list)
        for (short, pid), d in P[s].items():
            r = d["row"]
            if not sel(r) or r["contig"] == "NA" or r["gene_rank_on_contig"] == "":
                continue
            byc[(short, r["contig"])].append((int(r["gene_rank_on_contig"]), d["labels"]))
        for a in LABELS:
            n_a = 0
            near = Counter()
            for k, lst in byc.items():
                for (i, la) in lst:
                    if a not in la:
                        continue
                    n_a += 1
                    nb = set()
                    for (j, lb) in lst:
                        if j != i and abs(j - i) <= 5:
                            nb |= lb
                    for b in nb:
                        near[b] += 1
            for b in LABELS:
                rows.append(dict(genus=("Aspergillus fumigatus (293)" if s == "afum" else "Coccidioides (529)"),
                                 label_a=a, label_b=b, n_genes_a_placed=n_a, n_with_b_within_5=near[b],
                                 frac=round(near[b] / n_a, 4) if n_a else ""))
    write_tsv(T / "allstrain_neighbour_rates.tsv", rows)


if __name__ == "__main__":
    main()
