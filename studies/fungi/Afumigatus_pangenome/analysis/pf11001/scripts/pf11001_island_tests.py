#!/usr/bin/env python3
"""PF11001 island-enrichment check, stage 1.

For each of the three pangenome runs:
  * Re-implement the pipeline's island Pfam Fisher test (one-sided, BH-FDR)
    from the raw tables (pfam.domtblout, frequency_table.tsv,
    significant_islands.tsv) and compare with island_pfam_enrichment.tsv.
  * Repeat the same test with the --cut_ga hmmsearch family annotation
    (analysis/gain_loss_enrichment/<run>/family_pfam.tsv.gz), same background.
  * Write one row per PF11001 family (both annotation sources) with class,
    carrier count, island membership, gain/loss direction and member counts.
  * Cross-run selection-effect numbers and Fisher's combined p.

Outputs go under analysis/pf11001/tables/ of each study.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import combine_pvalues, false_discovery_control, fisher_exact
from scipy.stats.contingency import odds_ratio

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
RUNS = {
    "full_v072": dict(
        study=AF, pg=AF / "results/full_v072/config/pangenome",
        fam_pfam=AF / "analysis/gain_loss_enrichment/full_v072/family_pfam.tsv.gz",
        gl=AF / "analysis/gain_loss/full_v072/family_gain_loss.tsv.gz",
        out=AF / "analysis/pf11001/tables"),
    "immitis_in_posadasii_out_v5": dict(
        study=CP, pg=CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome",
        fam_pfam=CP / "analysis/gain_loss_enrichment/immitis_in_posadasii_out_v5/family_pfam.tsv.gz",
        gl=CP / "analysis/gain_loss/immitis_in_posadasii_out_v5/family_gain_loss.tsv.gz",
        out=CP / "analysis/pf11001/tables"),
    "posadasii_in_immitis_out_v2": dict(
        study=CP, pg=CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome",
        fam_pfam=CP / "analysis/gain_loss_enrichment/posadasii_in_immitis_out_v2/family_pfam.tsv.gz",
        gl=CP / "analysis/gain_loss/posadasii_in_immitis_out_v2/family_gain_loss.tsv.gz",
        out=CP / "analysis/pf11001/tables"),
}
DOM_NAME = "AFUB_07903_YDR124W_hel"
DOM_ACC = "PF11001"


def read_domtbl_scan(path, max_ie=1e-3):
    """hmmscan domtblout: target=Pfam name (col0), acc (col1), query=family (col3),
    i-Evalue col12, domain score col13."""
    fam_doms = defaultdict(set)
    hits = []  # PF11001 rows
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            ie = float(p[12])
            if p[0] == DOM_NAME:
                hits.append(dict(family=p[3], qlen=int(p[5]), full_E=p[6], full_score=p[7],
                                 i_Evalue=p[12], dom_score=p[13], ali_from=p[17], ali_to=p[18],
                                 env_from=p[19], env_to=p[20], passes_ie=ie <= max_ie,
                                 passes_ga=float(p[13]) >= 31.3))
            if ie > max_ie:
                continue
            fam_doms[p[3]].add(p[0])
    return fam_doms, hits


def read_family_pfam(path):
    fam_doms = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["pfam_ids"]:
                fam_doms[r["accession"]] = set(r["pfam_ids"].split("|"))
    return fam_doms


def read_freq(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["family"]] = r
    return out


def read_islands(path):
    isl = []
    with open(path, newline="") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            isl.append(dict(idx=i, n_strains=int(r["n_strains"]),
                            fams=r["member_families"].split(",")))
    return isl


def enrichment(island_fams, background, fam_doms):
    island = island_fams & background
    n_bg, n_is = len(background), len(island)
    doms = set()
    for f in island:
        doms |= fam_doms.get(f, set())
    dom_fams = defaultdict(set)
    for f in background:
        for d in fam_doms.get(f, ()):
            dom_fams[d].add(f)
    rows = []
    for d in sorted(doms):
        fd = dom_fams[d]
        a = len(fd & island)
        b = n_is - a
        c = len(fd) - a
        e = n_bg - n_is - c
        _, p = fisher_exact([[a, b], [c, e]], alternative="greater")
        rows.append(dict(domain=d, a=a, K=len(fd), n_island=n_is, n_bg=n_bg, p=p,
                         table=(a, b, c, e)))
    q = false_discovery_control(np.array([r["p"] for r in rows]), method="bh")
    for r, qq in zip(rows, q):
        r["q"] = float(qq)
    rows.sort(key=lambda r: (r["p"], r["domain"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def published(path):
    with open(path, newline="") as fh:
        return {r["domain"]: r for r in csv.DictReader(fh, delimiter="\t")}


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def effect(row):
    a, b, c, e = row["table"]
    res = odds_ratio([[a, b], [c, e]], kind="conditional")
    ci = res.confidence_interval(0.95)
    frac = a / row["K"] if row["K"] else float("nan")
    exp = row["n_island"] / row["n_bg"]
    return dict(frac_dom_fams_in_islands=round(frac, 4), expected_frac=round(exp, 4),
                fold=round(frac / exp, 3), expected_count=round(row["K"] * exp, 2),
                odds_ratio_cond=round(float(res.statistic), 3),
                or_ci95_lo=round(float(ci.low), 3), or_ci95_hi=(round(float(ci.high), 3)
                                                                 if np.isfinite(ci.high) else "inf"))


def main():
    summary = []
    per_run_rows = {}
    for run, cfg in RUNS.items():
        pg = cfg["pg"]
        out = cfg["out"]
        fam_doms, pf_hits = read_domtbl_scan(pg / "pfam.domtblout")
        freq = read_freq(pg / "frequency_table.tsv")
        background = {f for f, r in freq.items() if r["bin"] in ("shell", "cloud")}
        islands = read_islands(pg / "significant_islands.tsv")
        island_fams = set()
        for i in islands:
            island_fams.update(i["fams"])
        rows = enrichment(island_fams, background, fam_doms)
        pub = published(pg / "island_pfam_enrichment.tsv")
        # compare with published
        n_diff_counts = 0
        max_rel_p = 0.0
        max_rel_q = 0.0
        for r in rows:
            pr = pub.get(r["domain"])
            if pr is None:
                n_diff_counts += 1
                continue
            if (int(pr["n_with_domain_in_islands"]), int(pr["n_with_domain_in_background"]),
                    int(pr["n_island_families"]), int(pr["n_background_families"])) != (
                    r["a"], r["K"], r["n_island"], r["n_bg"]):
                n_diff_counts += 1
            pp, pq = float(pr["fisher_p"]), float(pr["fdr_q"])
            if pp > 0:
                max_rel_p = max(max_rel_p, abs(r["p"] - pp) / pp)
            if pq > 0:
                max_rel_q = max(max_rel_q, abs(r["q"] - pq) / pq)
        missing_in_mine = len(set(pub) - {r["domain"] for r in rows})
        write_tsv(out / f"island_fisher_reproduced.{run}.tsv.gz",
                  ["rank", "domain", "a", "K", "n_island", "n_bg", "p", "q"], rows)
        me = next(r for r in rows if r["domain"] == DOM_NAME)
        pubme = pub[DOM_NAME]
        n_sig = sum(r["q"] < 0.05 for r in rows)
        # neighbours by rank
        neigh = [f"{r['rank']}:{r['domain']}({r['a']}/{r['K']},p={r['p']:.2e},q={r['q']:.3f})"
                 for r in rows if abs(r["rank"] - me["rank"]) <= 3]

        # cut_ga sensitivity
        fp = read_family_pfam(cfg["fam_pfam"])
        rows_ga = enrichment(island_fams, background, fp)
        me_ga = next((r for r in rows_ga if r["domain"] == DOM_ACC), None)
        write_tsv(out / f"island_fisher_cutga.{run}.tsv.gz",
                  ["rank", "domain", "a", "K", "n_island", "n_bg", "p", "q"], rows_ga)

        s = dict(run=run, annotation="pipeline hmmscan --domE 1e-3, i-Evalue<=1e-3, shell+cloud reps",
                 n_domains_tested=len(rows), n_q_lt_0p05=n_sig, pf11001_rank_by_p=me["rank"],
                 a=me["a"], K=me["K"], n_island=me["n_island"], n_bg=me["n_bg"],
                 p=f"{me['p']:.4e}", q=f"{me['q']:.4e}", q_lt_0p05=me["q"] < 0.05,
                 published_counts=f"{pubme['n_with_domain_in_islands']}/{pubme['n_with_domain_in_background']}"
                                  f" of {pubme['n_island_families']}/{pubme['n_background_families']}",
                 published_p=pubme["fisher_p"], published_q=pubme["fdr_q"],
                 n_domains_count_mismatch=n_diff_counts, n_published_not_reproduced=missing_in_mine,
                 max_rel_diff_p=f"{max_rel_p:.2e}", max_rel_diff_q=f"{max_rel_q:.2e}",
                 neighbours=" | ".join(neigh), **effect(me))
        summary.append(s)
        if me_ga:
            sg = dict(run=run, annotation="our hmmsearch --cut_ga on all family reps (family_pfam.tsv.gz); same shell+cloud background",
                      n_domains_tested=len(rows_ga), n_q_lt_0p05=sum(r["q"] < 0.05 for r in rows_ga),
                      pf11001_rank_by_p=me_ga["rank"], a=me_ga["a"], K=me_ga["K"],
                      n_island=me_ga["n_island"], n_bg=me_ga["n_bg"], p=f"{me_ga['p']:.4e}",
                      q=f"{me_ga['q']:.4e}", q_lt_0p05=me_ga["q"] < 0.05, **effect(me_ga))
            summary.append(sg)
        per_run_rows[run] = dict(pipe={r["domain"]: r for r in rows},
                                 ga={r["domain"]: r for r in rows_ga})

        # ---- per-family table ----
        fam_isl = defaultdict(list)
        for i in islands:
            for f in i["fams"]:
                fam_isl[f].append(i)
        members = defaultdict(list)
        with open(pg / "cluster/tier1_cluster.tsv") as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")
                members[rep].append(mem)
        gl = {}
        with gzip.open(cfg["gl"], "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gl[r["family"]] = r
        pipe_fams = {h["family"]: h for h in pf_hits}
        ga_fams = {f for f, d in fp.items() if DOM_ACC in d}
        allf = sorted(set(pipe_fams) | ga_fams)
        frows = []
        for f in allf:
            fr = freq.get(f, {})
            isl = fam_isl.get(f, [])
            mem = members.get(f, [])
            strains = {m.split("|", 1)[0] for m in mem}
            h = pipe_fams.get(f)
            g = gl.get(f, {})
            frows.append(dict(
                run=run, family=f, bin=fr.get("bin", "NA"), strain_count=fr.get("strain_count", "NA"),
                frequency=fr.get("frequency", "NA"), bin_out=fr.get("bin_out", "NA"),
                in_pipeline_background=f in background,
                pipeline_hmmscan_hit=h is not None,
                pipeline_i_Evalue=h["i_Evalue"] if h else "", pipeline_dom_score=h["dom_score"] if h else "",
                pipeline_counted=(h is not None and h["passes_ie"] and f in background),
                cutga_hmmsearch_hit=f in ga_fams,
                n_islands=len(isl),
                n_islands_single_strain=sum(i["n_strains"] == 1 for i in isl),
                n_islands_multi_strain=sum(i["n_strains"] >= 2 for i in isl),
                max_island_n_strains=max((i["n_strains"] for i in isl), default=0),
                island_idx=",".join(str(i["idx"]) for i in isl[:50]) + (",..." if len(isl) > 50 else ""),
                gl_class=g.get("class", "NA"), direction=g.get("direction", "NA"),
                outgroup_count=g.get("outgroup_count", "NA"),
                n_member_proteins=len(mem), n_member_strains=len(strains),
            ))
        per_run_rows[run]["families"] = frows
        per_run_rows[run]["islands"] = islands
        per_run_rows[run]["fam_isl"] = fam_isl
        hdr = list(frows[0].keys())
        write_tsv(out / f"pf11001_families.{run}.tsv", hdr, frows)

        # independent islands: islands with >=1 PF11001 family (pipeline-counted or cut_ga)
        pf_set = {r["family"] for r in frows}
        isl_pf = [i for i in islands if pf_set & set(i["fams"])]
        # union-find on islands sharing any member family
        parent = {i["idx"]: i["idx"] for i in isl_pf}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        by_fam = defaultdict(list)
        for i in isl_pf:
            for f in i["fams"]:
                by_fam[f].append(i["idx"])
        for f, idxs in by_fam.items():
            for j in idxs[1:]:
                ra, rb = find(idxs[0]), find(j)
                if ra != rb:
                    parent[ra] = rb
        comps = defaultdict(list)
        for i in isl_pf:
            comps[find(i["idx"])].append(i)
        n_pf_per_isl = [len(pf_set & set(i["fams"])) for i in isl_pf]
        isl_rows = []
        for i in isl_pf:
            isl_rows.append(dict(run=run, island_idx=i["idx"], n_strains=i["n_strains"],
                                 n_member_families=len(i["fams"]),
                                 pf11001_families=",".join(sorted(pf_set & set(i["fams"]))),
                                 component=find(i["idx"])))
        write_tsv(out / f"pf11001_islands.{run}.tsv", list(isl_rows[0].keys()), isl_rows)
        s_isl = dict(run=run, n_pf11001_families_any_source=len(pf_set),
                     n_pf11001_families_in_any_island=sum(1 for f in pf_set if fam_isl.get(f)),
                     n_islands_with_pf11001=len(isl_pf),
                     n_islands_single_strain=sum(i["n_strains"] == 1 for i in isl_pf),
                     n_islands_multi_strain=sum(i["n_strains"] >= 2 for i in isl_pf),
                     n_islands_with_ge2_pf11001_families=sum(n >= 2 for n in n_pf_per_isl),
                     n_components_islands_sharing_any_family=len(comps),
                     largest_component_n_islands=max(len(c) for c in comps.values()))
        per_run_rows[run]["island_summary"] = s_isl

    # ---- cross-run selection effect ----
    cross = []
    af = per_run_rows["full_v072"]["pipe"]
    for crun in ("immitis_in_posadasii_out_v5", "posadasii_in_immitis_out_v2"):
        co = per_run_rows[crun]["pipe"]
        common = set(af) & set(co)
        paf, pco = af[DOM_NAME]["p"], co[DOM_NAME]["p"]
        n_both_005 = sum(af[d]["p"] <= 0.05 and co[d]["p"] <= 0.05 for d in common)
        n_af_005 = sum(af[d]["p"] <= 0.05 for d in common)
        n_co_005 = sum(co[d]["p"] <= 0.05 for d in common)
        n_both_le = sum(af[d]["p"] <= paf and co[d]["p"] <= pco for d in common)
        both_le = sorted(d for d in common if af[d]["p"] <= paf and co[d]["p"] <= pco)
        n_both_q = sum(af[d]["q"] < 0.05 and co[d]["q"] < 0.05 for d in common)
        both_q = sorted(d for d in common if af[d]["q"] < 0.05 and co[d]["q"] < 0.05)
        stat, pcomb = combine_pvalues([paf, pco], method="fisher")
        # Fisher combined p for every common domain, rank of PF11001
        comb = sorted(((combine_pvalues([af[d]["p"], co[d]["p"]], method="fisher")[1], d) for d in common))
        rank = [d for _, d in comb].index(DOM_NAME) + 1
        top = "; ".join(f"{d}:{p:.2e}" for p, d in comb[:8])
        cross.append(dict(pair=f"full_v072 x {crun}", n_domains_af=len(af), n_domains_cocci=len(co),
                          n_domains_tested_in_both=len(common),
                          n_af_p_le_0p05=n_af_005, n_cocci_p_le_0p05=n_co_005,
                          n_both_p_le_0p05=n_both_005,
                          expected_both_p_le_0p05_if_independent=round(n_af_005 * n_co_005 / len(common), 2),
                          n_both_p_le_pf11001_p=n_both_le, domains_both_p_le_pf11001_p=",".join(both_le),
                          n_both_q_lt_0p05=n_both_q, domains_both_q_lt_0p05=",".join(both_q),
                          fisher_combined_p_pf11001=f"{pcomb:.3e}", fisher_chi2=round(float(stat), 3),
                          pf11001_rank_by_combined_p=rank, top_combined=top))
    # all three (v5 and v2 share genomes: not independent; shown only for reference)
    p3 = [per_run_rows[r]["pipe"][DOM_NAME]["p"] for r in RUNS]
    cross.append(dict(pair="all three (NOT independent: v5 and v2 share genomes)",
                      fisher_combined_p_pf11001=f"{combine_pvalues(p3, method='fisher')[1]:.3e}"))

    out = AF / "analysis/pf11001/tables"
    hdr = []
    for s in summary:
        for k in s:
            if k not in hdr:
                hdr.append(k)
    write_tsv(out / "island_enrichment_summary.tsv", hdr, summary)
    hdr = []
    for s in cross:
        for k in s:
            if k not in hdr:
                hdr.append(k)
    write_tsv(out / "cross_run_selection.tsv", hdr, cross)
    isum = [per_run_rows[r]["island_summary"] for r in RUNS]
    write_tsv(out / "island_membership_summary.tsv", list(isum[0].keys()), isum)
    for s in summary:
        print(json.dumps(s, default=str))
    for s in cross:
        print(json.dumps(s, default=str))
    for s in isum:
        print(json.dumps(s, default=str))


if __name__ == "__main__":
    main()
