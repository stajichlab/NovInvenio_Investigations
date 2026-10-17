#!/usr/bin/env python3
"""Cross-genus best-hit homologous families: A. fumigatus full_v072 vs Coccidioides v5, v2.

Input: mmseqs easy-search tables (mmseqs/<a>_vs_<b>.m8.zst; run_mmseqs_search.sbatch:
-e 1e-5 -c 0.5 --cov-mode 0). Best hit per query = highest bits (ties: lower evalue,
then target id). Reciprocal best hit (RBH) = Af family A and Coccidioides family C with
best(A) = C and best(C) = A. v5 and v2 are handled separately (no v5-v2 mapping).

Per family attributes: bin and strain_count (frequency_table.tsv), island membership
(union of significant_islands.tsv member_families), gain/loss direction (rule_freq of
analysis/gain_loss/<run>/family_gain_loss.tsv.gz), Pfam names (family_pfam.tsv.gz),
Af Starship element / captain status (analysis/starship_island_overlap/full_v072/
annotated_families.tsv.gz).

Tests (Fisher exact, two-sided, plus one-sided greater; conditional OR with 95% CI):
  island: Af family in island (yes/no) x RBH partner in Coccidioides island (yes/no),
          over (a) all RBH pairs, (b) pairs with both families in bin shell or cloud.
  gain:   Af direction gain (yes/no) x partner direction gain (yes/no), same two sets.
  loss:   same with loss.
"""
from __future__ import annotations

import csv
import gzip
import io
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
OUT = CP / "analysis/cross_species_overrep"
RUNS = {
    "af": dict(name="full_v072", study=AF, pg=AF / "results/full_v072/config/pangenome"),
    "v5": dict(name="immitis_in_posadasii_out_v5", study=CP,
               pg=CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome"),
    "v2": dict(name="posadasii_in_immitis_out_v2", study=CP,
               pg=CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome"),
}
STARSHIP = AF / "analysis/starship_island_overlap/full_v072/annotated_families.tsv.gz"
SC = ("shell", "cloud")


def zst_lines(path):
    proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE)
    yield from io.TextIOWrapper(proc.stdout)
    proc.wait()


def best_hits(path):
    best, n_q_hits = {}, set()
    for line in zst_lines(path):
        p = line.rstrip("\n").split("\t")
        q, t, ev, bits = p[0], p[1], float(p[10]), float(p[11])
        n_q_hits.add(q)
        key = (-bits, ev, t)
        if q not in best or key < best[q][0]:
            best[q] = (key, t, float(p[2]), float(p[14]), float(p[15]), ev, bits)
    return {q: v[1:] for q, v in best.items()}


def load_run(k):
    cfg = RUNS[k]
    fam = {}
    with open(cfg["pg"] / "frequency_table.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam[r["family"]] = dict(bin=r["bin"], strain_count=r["strain_count"], island=False,
                                    direction="", pfam="")
    with open(cfg["pg"] / "significant_islands.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for f in r["member_families"].split(","):
                if f in fam:
                    fam[f]["island"] = True
    gl = cfg["study"] / "analysis/gain_loss" / cfg["name"] / "family_gain_loss.tsv.gz"
    with gzip.open(gl, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family"] in fam:
                fam[r["family"]]["direction"] = r["rule_freq"]
    fp = cfg["study"] / "analysis/gain_loss_enrichment" / cfg["name"] / "family_pfam.tsv.gz"
    with gzip.open(fp, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] in fam:
                fam[r["accession"]]["pfam"] = r["pfam_names"]
    return fam


def fisher_row(label, pairs, fa, fc, pred_a, pred_c):
    a = sum(pred_a(fa[x]) and pred_c(fc[y]) for x, y in pairs)
    b = sum(pred_a(fa[x]) and not pred_c(fc[y]) for x, y in pairs)
    c = sum((not pred_a(fa[x])) and pred_c(fc[y]) for x, y in pairs)
    d = sum((not pred_a(fa[x])) and not pred_c(fc[y]) for x, y in pairs)
    _, p2 = fisher_exact([[a, b], [c, d]])
    _, p1 = fisher_exact([[a, b], [c, d]], alternative="greater")
    res = odds_ratio([[a, b], [c, d]], kind="conditional")
    ci = res.confidence_interval(0.95)
    f1 = a / (a + b) if a + b else float("nan")
    f0 = c / (c + d) if c + d else float("nan")
    return dict(test=label, n_pairs=len(pairs), af_yes_partner_yes=a, af_yes_partner_no=b,
                af_no_partner_yes=c, af_no_partner_no=d,
                frac_partner_yes_given_af_yes=round(f1, 4),
                frac_partner_yes_given_af_no=round(f0, 4),
                ratio=round(f1 / f0, 3) if f0 else "inf",
                odds_ratio=round(float(res.statistic), 3), or_ci95_lo=round(float(ci.low), 3),
                or_ci95_hi=round(float(ci.high), 3) if np.isfinite(ci.high) else "inf",
                p_two_sided=f"{p2:.3e}", p_greater=f"{p1:.3e}")


def components(edges):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    comp = defaultdict(set)
    for x in list(parent):
        comp[find(x)].add(x)
    return list(comp.values())


def write_tsv(path, rows, header=None):
    header = header or list(rows[0].keys())
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def main():
    fams = {k: load_run(k) for k in RUNS}
    star = {}
    with gzip.open(STARSHIP, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            star[r["family"]] = (r["starship_elements"], r["captain"])
    md = OUT / "mmseqs"
    pair_rows, tests, cover, groups = [], [], [], []
    rbh = {}
    for c in ("v5", "v2"):
        bh_a = best_hits(md / f"af_vs_{c}.m8.zst")
        bh_c = best_hits(md / f"{c}_vs_af.m8.zst")
        fa, fc = fams["af"], fams[c]
        pairs = []
        for a, (t, fid, qcov, tcov, ev, bits) in bh_a.items():
            recip = bh_c.get(t, (None,))[0] == a
            if recip:
                pairs.append((a, t))
            pair_rows.append(dict(cocci_run=RUNS[c]["name"], af_family=a, cocci_family=t,
                                  reciprocal=recip, fident=fid, qcov=qcov, tcov=tcov, evalue=ev,
                                  bits=bits, af_bin=fa[a]["bin"], af_strains=fa[a]["strain_count"],
                                  af_island=fa[a]["island"], af_direction=fa[a]["direction"],
                                  cocci_bin=fc[t]["bin"], cocci_strains=fc[t]["strain_count"],
                                  cocci_island=fc[t]["island"], cocci_direction=fc[t]["direction"],
                                  af_pfam=fa[a]["pfam"], cocci_pfam=fc[t]["pfam"]))
        rbh[c] = dict(pairs)
        # coverage of partners by Af bin and island status
        for b in ["core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "outgroup_only",
                  "ALL", "island", "shell+cloud island", "shell+cloud non-island"]:
            if b == "ALL":
                sel = list(fa)
            elif b == "island":
                sel = [f for f in fa if fa[f]["island"]]
            elif b == "shell+cloud island":
                sel = [f for f in fa if fa[f]["island"] and fa[f]["bin"] in SC]
            elif b == "shell+cloud non-island":
                sel = [f for f in fa if not fa[f]["island"] and fa[f]["bin"] in SC]
            else:
                sel = [f for f in fa if fa[f]["bin"] == b]
            n = len(sel)
            nb = sum(f in bh_a for f in sel)
            nr = sum(f in rbh[c] for f in sel)
            cover.append(dict(cocci_run=RUNS[c]["name"], af_set=b, n_af_families=n,
                              n_with_best_hit=nb, frac_with_best_hit=round(nb / n, 4) if n else "",
                              n_with_rbh=nr, frac_with_rbh=round(nr / n, 4) if n else ""))
        # tests
        sc_pairs = [(a, t) for a, t in pairs if fa[a]["bin"] in SC and fc[t]["bin"] in SC]
        for lab, ps in (("all_rbh_pairs", pairs), ("both_shell_or_cloud", sc_pairs)):
            for what, pa, pc in (
                    ("island", lambda x: x["island"], lambda y: y["island"]),
                    ("gain", lambda x: x["direction"] == "gain", lambda y: y["direction"] == "gain"),
                    ("loss", lambda x: x["direction"] == "loss", lambda y: y["direction"] == "loss")):
                r = fisher_row(f"{what}|{lab}", ps, fa, fc, pa, pc)
                r = dict(cocci_run=RUNS[c]["name"], **r)
                tests.append(r)
        # paralog / group check: best-hit graph components (both directions, not only RBH)
        edges = [("A:" + a, "C:" + v[0]) for a, v in bh_a.items()] + \
                [("C:" + q, "A:" + v[0]) for q, v in bh_c.items()]
        comps = components(edges)
        comps.sort(key=len, reverse=True)
        node2comp = {}
        for i, cs in enumerate(comps):
            for x in cs:
                node2comp[x] = i
        rbh_per_comp = Counter(node2comp["A:" + a] for a, _ in pairs)
        indeg = Counter(v[0] for v in bh_a.values())
        for i, cs in enumerate(comps[:15]):
            na = sum(x.startswith("A:") for x in cs)
            afs = [x[2:] for x in cs if x.startswith("A:")]
            pf = Counter()
            for f in afs:
                for n in filter(None, fa[f]["pfam"].split("|")):
                    pf[n] += 1
            groups.append(dict(cocci_run=RUNS[c]["name"], rank=i + 1, n_nodes=len(cs), n_af=na,
                               n_cocci=len(cs) - na, n_rbh_pairs=rbh_per_comp.get(i, 0),
                               n_af_island=sum(fa[f]["island"] for f in afs),
                               top_af_pfam=";".join(f"{n}:{m}" for n, m in pf.most_common(4))))
        top_in = indeg.most_common(5)
        sizes = [len(cs) for cs in comps]
        print(f"[{c}] components with an RBH pair {len(rbh_per_comp)}, max RBH pairs per component {max(rbh_per_comp.values())}", file=sys.stderr)
        print(f"[{c}] best hits Af->C {len(bh_a)}, C->Af {len(bh_c)}, RBH pairs {len(pairs)} "
              f"(distinct Af {len(set(a for a,_ in pairs))}, distinct C {len(set(t for _,t in pairs))}); "
              f"components {len(comps)}, largest {sizes[:5]}; RBH pairs in the 5 largest comps "
              f"{sum(rbh_per_comp.get(i,0) for i in range(5))}; max in-degree of Af best hits {top_in}",
              file=sys.stderr)
        # Note: in the best-hit graph each node has one out-edge, so every component holds
        # at most one 2-cycle (= one RBH pair). RBH pairs are therefore 1:1 and a large
        # paralog group contributes at most one pair.
        # Sensitivity: Af family -> its (not necessarily reciprocal) best hit; Af families
        # are the units, so one Coccidioides family can be counted many times.
        bh_pairs = [(a, v[0]) for a, v in bh_a.items()]
        bh_sc = [(a, t) for a, t in bh_pairs if fa[a]["bin"] in SC and fc[t]["bin"] in SC]
        for lab, ps in (("besthit_af_to_cocci_all", bh_pairs), ("besthit_af_to_cocci_both_shell_or_cloud", bh_sc)):
            r = fisher_row(f"island|{lab}", ps, fa, fc, lambda x: x["island"], lambda y: y["island"])
            tests.append(dict(cocci_run=RUNS[c]["name"], **r))

    # cross-genus island families
    fa = fams["af"]
    xg = []
    for a in sorted(fa):
        if not fa[a]["island"]:
            continue
        t5, t2 = rbh["v5"].get(a), rbh["v2"].get(a)
        if t5 and t2 and fams["v5"][t5]["island"] and fams["v2"][t2]["island"]:
            # annotated_families.tsv.gz lists only families in an Af293 Table S7 element or
            # with a DUF3435 captain hit (380 rows); others are "not_listed".
            se, cap = star.get(a, ("not_listed", "not_listed"))
            xg.append(dict(af_family=a, af_bin=fa[a]["bin"], af_strains=fa[a]["strain_count"],
                           af_direction=fa[a]["direction"], af_pfam=fa[a]["pfam"],
                           af_starship_elements=se, af_captain=cap,
                           v5_family=t5, v5_bin=fams["v5"][t5]["bin"],
                           v5_strains=fams["v5"][t5]["strain_count"],
                           v5_direction=fams["v5"][t5]["direction"], v5_pfam=fams["v5"][t5]["pfam"],
                           v2_family=t2, v2_bin=fams["v2"][t2]["bin"],
                           v2_strains=fams["v2"][t2]["strain_count"],
                           v2_direction=fams["v2"][t2]["direction"], v2_pfam=fams["v2"][t2]["pfam"]))
    # three-way: Af families with an RBH partner in both runs; Af island x (both partners island)
    both = [a for a in fa if a in rbh["v5"] and a in rbh["v2"]]
    trip = [(a, a) for a in both]
    fboth = {a: dict(island=fams["v5"][rbh["v5"][a]]["island"] and fams["v2"][rbh["v2"][a]]["island"])
             for a in both}
    r = fisher_row("island|rbh_in_both_runs|partners_island_in_both", trip, fa, fboth,
                   lambda x: x["island"], lambda y: y["island"])
    tests.append(dict(cocci_run="v5_and_v2", **r))
    both_sc = [a for a in both if fa[a]["bin"] in SC and fams["v5"][rbh["v5"][a]]["bin"] in SC
               and fams["v2"][rbh["v2"][a]]["bin"] in SC]
    r = fisher_row("island|rbh_in_both_runs_all_shell_or_cloud|partners_island_in_both",
                   [(a, a) for a in both_sc], fa, fboth, lambda x: x["island"], lambda y: y["island"])
    tests.append(dict(cocci_run="v5_and_v2", **r))
    write_tsv(OUT / "family_pairs.tsv.gz", pair_rows)
    write_tsv(OUT / "family_tests.tsv", tests)
    write_tsv(OUT / "family_partner_coverage.tsv", cover)
    write_tsv(OUT / "family_besthit_components_top15.tsv", groups)
    write_tsv(OUT / "cross_genus_island_families.tsv", xg)
    print(f"cross-genus island families: {len(xg)}", file=sys.stderr)


if __name__ == "__main__":
    main()
