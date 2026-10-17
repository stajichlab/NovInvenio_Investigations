#!/usr/bin/env python3
"""Cross-species Pfam over-representation: A. fumigatus full_v072 vs C. immitis v5 and
C. posadasii v2.

Definitions (one set of tests per run; each run tested separately):
  A_island_pfam / A_island_clan
      Pipeline island test. Background = families with frequency_table bin shell or cloud.
      Family has a domain when one pfam.domtblout row (hmmscan --domE 1e-3 on the
      background representatives) has i-Evalue <= 1e-3. Island set = union of
      significant_islands.tsv member_families, intersected with background. One-sided
      Fisher (greater); BH over the terms present in >= 1 island family.
      Pfam level is recomputed here and checked against island_pfam_enrichment.tsv
      (counts must match exactly). Clan level is new: a family has a clan when one of
      its counted domains maps to that clan (Pfam-A.hmm.dat CL), same test.
  B_<scope>_<dir>_<vocab>
      Read from analysis/gain_loss_enrichment/<run>/<scope>/enrichment_domains.tsv.gz
      (bin/pangenome_gain_loss_enrichment.py: hypergeometric, BH per test set and
      vocabulary, K >= 5). scope shell_soft_core or all_accessory; dir gain or loss;
      vocab pfam or pfam_clan.
  C_accessory_vs_all / C_shellcloud_vs_core
      New here. Annotation = analysis/gain_loss_enrichment/<run>/family_pfam.tsv.gz
      (hmmsearch --cut_ga of all family representatives). Universe:
        accessory_vs_all: families with bin in core, soft_core, shell, cloud, singleton
          and >= 1 Pfam; test set = bin soft_core, shell, cloud, singleton.
        shellcloud_vs_core: families with bin core, shell, cloud and >= 1 Pfam;
          test set = shell + cloud.
      p = hypergeom.sf(k-1, M, K, N); term tested if K >= 5; BH over tested terms.

Cross-run summary per definition: terms tested in >= 2 runs; tiers; hypergeometric
overlap test of q < 0.05 terms; Spearman of -log10 p over shared terms.
"""
from __future__ import annotations

import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, fisher_exact, hypergeom, spearmanr

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
OUT = CP / "analysis/cross_species_overrep"
PFAM_DAT = Path("/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat")
PFAM_C = Path("/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-C")
RUNS = {
    "af": dict(name="full_v072", study=AF,
               pg=AF / "results/full_v072/config/pangenome"),
    "v5": dict(name="immitis_in_posadasii_out_v5", study=CP,
               pg=CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome"),
    "v2": dict(name="posadasii_in_immitis_out_v2", study=CP,
               pg=CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome"),
}
for r in RUNS.values():
    r["gle"] = r["study"] / "analysis/gain_loss_enrichment" / r["name"]
KEYS = ["af", "v5", "v2"]
ALPHA = 0.05


def acc_nov(a: str) -> str:
    return a.split(".")[0]


def read_dat():
    info, name2acc = {}, {}
    cur = {}
    with open(PFAM_DAT, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                cur = {"id": line[10:].strip()}
            elif line.startswith("#=GF AC"):
                cur["acc"] = acc_nov(line[10:].strip())
            elif line.startswith("#=GF DE"):
                cur["de"] = line[10:].strip()
            elif line.startswith("#=GF CL"):
                cur["clan"] = line[10:].strip()
            elif line.startswith("//") and cur.get("acc"):
                info[cur["acc"]] = cur
                name2acc[cur["id"]] = cur["acc"]
                cur = {}
    return info, name2acc


def read_clans():
    out, cur = {}, {}
    with open(PFAM_C, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                cur = {"id": line[10:].strip()}
            elif line.startswith("#=GF AC"):
                cur["acc"] = acc_nov(line[10:].strip())
            elif line.startswith("#=GF DE") and "acc" in cur and "de" not in cur:
                cur["de"] = line[10:].strip()
                out[cur["acc"]] = cur
    return out


CLANS: dict = {}


def read_freq(pg):
    with open(pg / "frequency_table.tsv", newline="") as fh:
        return {r["family"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_island_fams(pg):
    fams = set()
    with open(pg / "significant_islands.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fams.update(r["member_families"].split(","))
    return fams


def read_domtbl(pg, name2acc, max_ie=1e-3):
    fam = defaultdict(set)
    with open(pg / "pfam.domtblout") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if float(p[12]) > max_ie:
                continue
            fam[p[3]].add(acc_nov(p[1]) if p[1] != "-" else name2acc[p[0]])
    return fam


def read_family_pfam(path):
    fam = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["pfam_ids"]:
                fam[r["accession"]] = {acc_nov(x) for x in r["pfam_ids"].split("|")}
    return fam


def to_clans(fam_terms, info):
    out = {}
    for f, ts in fam_terms.items():
        cl = {info[t]["clan"] for t in ts if t in info and info[t].get("clan")}
        if cl:
            out[f] = cl
    return out


def bh(rows):
    if not rows:
        return
    q = false_discovery_control(np.array([r["p"] for r in rows]), method="bh")
    for r, qq in zip(rows, q):
        r["q"] = float(qq)


def island_test(island, background, fam_terms):
    island = island & background
    n_bg, n_is = len(background), len(island)
    tf = defaultdict(set)
    for f in background:
        for t in fam_terms.get(f, ()):
            tf[t].add(f)
    present = set()
    for f in island:
        present |= fam_terms.get(f, set())
    rows = []
    for t in sorted(present):
        a = len(tf[t] & island)
        K = len(tf[t])
        _, p = fisher_exact([[a, n_is - a], [K - a, n_bg - n_is - K + a]], alternative="greater")
        rows.append(dict(term=t, k=a, K=K, N=n_is, M=n_bg, p=float(p),
                         fold=(a / K) / (n_is / n_bg)))
    bh(rows)
    return {r["term"]: r for r in rows}


def hyper_test(set_fams, universe, fam_terms, min_k=5):
    N, M = len(set_fams & universe), len(universe)
    tf = defaultdict(set)
    for f in universe:
        for t in fam_terms[f]:
            tf[t].add(f)
    rows = []
    for t, fs in sorted(tf.items()):
        K = len(fs)
        if K < min_k:
            continue
        k = len(fs & set_fams)
        p = float(hypergeom.sf(k - 1, M, K, N)) if k > 0 else 1.0
        rows.append(dict(term=t, k=k, K=K, N=N, M=M, p=p, fold=k / (K * N / M)))
    bh(rows)
    return {r["term"]: r for r in rows}


def read_gle(path, test_set, vocab):
    out = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["test_set"] != test_set or r["vocab"] != vocab:
                continue
            t = acc_nov(r["term"]) if vocab == "pfam" else r["term"]
            out[t] = dict(term=t, k=int(r["k"]), K=int(r["K"]), N=int(r["n_set"]),
                          M=int(r["M"]), p=float(r["p"]), q=float(r["q"]), fold=float(r["fold"]))
    return out


def fmt(x, d=4):
    if x is None or x == "":
        return ""
    if isinstance(x, float):
        if x == 0 or (abs(x) >= 1e-3 and abs(x) < 1e4):
            return f"{x:.{d}g}" if abs(x) < 1 else f"{x:.{d}g}"
        return f"{x:.3e}"
    return str(x)


def cross_run(defn, res, info, fam_all, vocab):
    """res: key -> {term: row}. Returns (table rows, summary dict)."""
    terms = set().union(*[set(res[k]) for k in KEYS])
    rows = []
    for t in sorted(terms):
        n_tested = sum(t in res[k] for k in KEYS)
        if n_tested < 2:
            continue
        sig = {k: (t in res[k] and res[k][t]["q"] < ALPHA) for k in KEYS}
        nom = {k: (t in res[k] and res[k][t]["p"] < ALPHA) for k in KEYS}
        if sig["af"] and sig["v5"] and sig["v2"]:
            tier = "1"
        elif sig["af"] and (sig["v5"] or sig["v2"]):
            tier = "2"
        elif nom["af"] and nom["v5"] and nom["v2"]:
            tier = "3"
        else:
            tier = ""
        if vocab == "pfam":
            inf = info.get(t, {})
            name, de, clan = inf.get("id", ""), inf.get("de", ""), inf.get("clan", "")
        else:
            members = sorted(info[a]["id"] for a in info if info[a].get("clan") == t)
            cinf = CLANS.get(t, {})
            name, clan = cinf.get("id", t), t
            de = f"{cinf.get('de', '')} ({len(members)} Pfam-A members)"
        row = dict(definition=defn, term=t, name=name, description=de, clan=clan,
                   n_runs_tested=n_tested, tier=tier)
        for k in KEYS:
            r = res[k].get(t)
            row[f"{k}_k"] = r["k"] if r else ""
            row[f"{k}_K"] = r["K"] if r else ""
            row[f"{k}_N"] = r["N"] if r else ""
            row[f"{k}_M"] = r["M"] if r else ""
            row[f"{k}_fold"] = round(r["fold"], 3) if r else ""
            row[f"{k}_p"] = r["p"] if r else ""
            row[f"{k}_q"] = r["q"] if r else ""
            row[f"{k}_n_families_all"] = fam_all[k].get(t, 0)
        rows.append(row)

    s = dict(definition=defn)
    for k in KEYS:
        s[f"{k}_n_tested"] = len(res[k])
        s[f"{k}_n_q_lt_0.05"] = sum(r["q"] < ALPHA for r in res[k].values())
    s["n_tested_ge2_runs"] = len(rows)
    s["n_tested_all3"] = sum(r["n_runs_tested"] == 3 for r in rows)
    for tier in "123":
        s[f"n_tier{tier}"] = sum(r["tier"] == tier for r in rows)
    s["n_af_q_and_any_cocci_q_incl_tier1"] = s["n_tier1"] + s["n_tier2"]
    s["n_p_lt_0.05_all3_incl_tier1_2"] = sum(
        all(res[k].get(r["term"], {"p": 1})["p"] < ALPHA for k in KEYS) for r in rows)

    # overlap tests and Spearman
    def overlap(keys_c, label):
        U = [t for t in terms if t in res["af"] and all(t in res[k] for k in keys_c)]
        Un = len(U)
        a = sum(res["af"][t]["q"] < ALPHA for t in U)
        b = sum(all(res[k][t]["q"] < ALPHA for k in keys_c) for t in U)
        o = sum(res["af"][t]["q"] < ALPHA and all(res[k][t]["q"] < ALPHA for k in keys_c) for t in U)
        exp = a * b / Un if Un else float("nan")
        p = float(hypergeom.sf(o - 1, Un, a, b)) if (Un and o > 0) else (1.0 if Un else float("nan"))
        s[f"{label}_universe"] = Un
        s[f"{label}_n_af_sig"] = a
        s[f"{label}_n_cocci_sig"] = b
        s[f"{label}_overlap"] = o
        s[f"{label}_expected"] = round(exp, 3) if Un else ""
        s[f"{label}_hypergeom_p"] = p
    overlap(["v5"], "af_x_v5")
    overlap(["v2"], "af_x_v2")
    overlap(["v5", "v2"], "af_x_v5andv2")
    for c in ["v5", "v2"]:
        U = [t for t in terms if t in res["af"] and t in res[c]]
        if len(U) >= 3:
            x = [-math.log10(max(res["af"][t]["p"], 1e-300)) for t in U]
            y = [-math.log10(max(res[c][t]["p"], 1e-300)) for t in U]
            rho, p = spearmanr(x, y)
            s[f"spearman_af_{c}_rho"] = round(float(rho), 4)
            s[f"spearman_af_{c}_p"] = float(p)
            s[f"spearman_af_{c}_n"] = len(U)
    return rows, s


def write_tsv(path, rows, header=None):
    if not rows:
        return
    header = header or list(rows[0].keys())
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(fmt(r.get(h, "")) if isinstance(r.get(h), float) else str(r.get(h, ""))
                               for h in header) + "\n")


def main():
    info, name2acc = read_dat()
    CLANS.update(read_clans())
    data = {}
    check = []
    for k in KEYS:
        cfg = RUNS[k]
        freq = read_freq(cfg["pg"])
        isl = read_island_fams(cfg["pg"])
        scan = read_domtbl(cfg["pg"], name2acc)
        fp = read_family_pfam(cfg["gle"] / "family_pfam.tsv.gz")
        data[k] = dict(freq=freq, isl=isl, scan=scan, fp=fp)
        # copy-number context: families (all bins) with term, cut_ga annotation
        ctx_p, ctx_c = defaultdict(int), defaultdict(int)
        for f, ts in fp.items():
            for t in ts:
                ctx_p[t] += 1
            for c in {info[t]["clan"] for t in ts if t in info and info[t].get("clan")}:
                ctx_c[c] += 1
        data[k]["ctx_pfam"], data[k]["ctx_clan"] = ctx_p, ctx_c

    results = {}
    # A
    resA, resAc = {}, {}
    for k in KEYS:
        d = data[k]
        bg = {f for f, r in d["freq"].items() if r["bin"] in ("shell", "cloud")}
        resA[k] = island_test(d["isl"], bg, d["scan"])
        resAc[k] = island_test(d["isl"], bg, to_clans(d["scan"], info))
        # check vs published
        pub = {}
        with open(RUNS[k]["pg"] / "island_pfam_enrichment.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pub[acc_nov(r["pfam_accession"])] = r
        mism = sum((int(pub[t]["n_with_domain_in_islands"]), int(pub[t]["n_with_domain_in_background"]))
                   != (r["k"], r["K"]) for t, r in resA[k].items() if t in pub)
        missing = len(set(pub) ^ set(resA[k]))
        maxrel = max(abs(r["p"] - float(pub[t]["fisher_p"])) / float(pub[t]["fisher_p"])
                     for t, r in resA[k].items() if t in pub and float(pub[t]["fisher_p"]) > 0)
        check.append(dict(run=RUNS[k]["name"], n_published=len(pub), n_recomputed=len(resA[k]),
                          n_terms_not_in_both=missing, n_count_mismatch=mism,
                          max_rel_diff_p=f"{maxrel:.2e}"))
    results["A_island_pfam"] = (resA, "pfam")
    results["A_island_clan"] = (resAc, "pfam_clan")
    # B
    for scope in ("shell_soft_core", "all_accessory"):
        for dr in ("gain", "loss"):
            for vocab in ("pfam", "pfam_clan"):
                res = {k: read_gle(RUNS[k]["gle"] / scope / "enrichment_domains.tsv.gz", dr, vocab)
                       for k in KEYS}
                v = "pfam" if vocab == "pfam" else "clan"
                results[f"B_{scope}_{dr}_{v}"] = (res, vocab)
    # C
    resC1, resC2 = {}, {}
    for k in KEYS:
        d = data[k]
        bins = {f: r["bin"] for f, r in d["freq"].items()}
        u1 = {f for f, b in bins.items() if b in ("core", "soft_core", "shell", "cloud", "singleton")
              and f in d["fp"]}
        s1 = {f for f in u1 if bins[f] != "core"}
        resC1[k] = hyper_test(s1, u1, d["fp"])
        u2 = {f for f, b in bins.items() if b in ("core", "shell", "cloud") and f in d["fp"]}
        s2 = {f for f in u2 if bins[f] in ("shell", "cloud")}
        resC2[k] = hyper_test(s2, u2, d["fp"])
    results["C_accessory_vs_all_pfam"] = (resC1, "pfam")
    results["C_shellcloud_vs_core_pfam"] = (resC2, "pfam")

    summaries, tiers = [], []
    for defn, (res, vocab) in results.items():
        ctx = {k: data[k]["ctx_pfam" if vocab == "pfam" else "ctx_clan"] for k in KEYS}
        rows, s = cross_run(defn, res, info, ctx, vocab)
        rows.sort(key=lambda r: (r["tier"] or "9", r["af_p"] if r["af_p"] != "" else 1, r["term"]))
        write_tsv(OUT / f"domains_{defn}.tsv", rows)
        summaries.append(s)
        tiers += [r for r in rows if r["tier"]]
        # full per-run results (all tested terms) for PF11001 lookup etc.
    write_tsv(OUT / "tiers.tsv", tiers)
    hdr = []
    for s in summaries:
        for h in s:
            if h not in hdr:
                hdr.append(h)
    write_tsv(OUT / "summary_domains.tsv", summaries, hdr)
    write_tsv(OUT / "logs" / "check_island_recompute.tsv", check)
    for c in check:
        print(c, file=sys.stderr)


if __name__ == "__main__":
    main()
