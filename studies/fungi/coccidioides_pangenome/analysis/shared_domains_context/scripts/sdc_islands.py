#!/usr/bin/env python3
"""Shared-domain context, family/island level (tasks 2ii, 3, 4).

Family labels: family representative hmmsearch --cut_ga Pfam (analysis/gain_loss_enrichment/
<run>/family_pfam.tsv.gz). HrmA label = any of PF28515 / PF28646 / PF28647.
Islands: significant_islands.tsv (row index = island_idx, as in starship_island_overlap).

Writes (tables/):
  island_pair_cooccurrence.tsv         per run x label pair: observed, permutation mean, p.
                                       count_type any = islands with >= 1 family with a and >= 1 with b;
                                       distinct_families = same, but a and b on different families
                                       (diagonal: islands with >= 2 families with a);
                                       family_pairs = ordered (a-family, b-family) pairs, f != g, sharing
                                       >= 1 island (diagonal: unordered pairs). Column observed_islands
                                       holds the family-pair count for that type.
                                       Null: label sets of island-resident families shuffled among
                                       island-resident families (1,000 permutations, seed 20261017).
  island_content.<tag>.tsv.gz          every island with >= 1 tier-1 family
  island_tier1_count_summary.tsv       islands by number of tier-1 labels (1-4), per run
  island_largest_multidomain.tsv       10 largest islands with >= 2 tier-1 labels, per run
  island_composition_summary.tsv       per run: class mix and label mix over tier-1 islands
  island_other_pfam_top.tsv            per run: top 25 other Pfams over tier-1 islands
  family_label_counts.tsv              per run x label: families, bins, island-resident
  crossgenus_pairs_by_domain.tsv       task 4: Af families with label vs Coccidioides best-hit partner
  crossgenus_island_families_by_domain.tsv  labels in cross_genus_island_families.tsv
"""
from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
D = CP / "analysis/shared_domains_context"
T = D / "tables"
XS = CP / "analysis/cross_species_overrep"
RUNS = {
    "af": ("full_v072", AF / "results/full_v072/config/pangenome",
           AF / "analysis/gain_loss_enrichment/full_v072/family_pfam.tsv.gz",
           AF / "analysis/starship_island_overlap/full_v072/with_unplaced_af293_coords/per_island.tsv.gz"),
    "v5": ("immitis_in_posadasii_out_v5",
           CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome",
           CP / "analysis/gain_loss_enrichment/immitis_in_posadasii_out_v5/family_pfam.tsv.gz",
           CP / "analysis/starship_island_overlap/immitis_in_posadasii_out_v5/per_island.tsv.gz"),
    "v2": ("posadasii_in_immitis_out_v2",
           CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome",
           CP / "analysis/gain_loss_enrichment/posadasii_in_immitis_out_v2/family_pfam.tsv.gz",
           CP / "analysis/starship_island_overlap/posadasii_in_immitis_out_v2/per_island.tsv.gz"),
}
ACC2LAB = {"PF01636": "APH", "PF11917": "DUF3435", "PF12520": "DUF3723", "PF29418": "DUF3723_N",
           "PF11001": "PF11001", "PF28515": "HrmA", "PF28646": "HrmA", "PF28647": "HrmA",
           "PF24883": "NPHP3_N", "PF00668": "Condensation", "PF00550": "PP-binding",
           "PF00501": "AMP-binding", "PF27653": "SPOK2_N"}
NAME2LAB = {"APH": "APH", "DUF3435": "DUF3435", "DUF3723": "DUF3723", "DUF3723_N": "DUF3723_N",
            "AFUB_07903_YDR124W_hel": "PF11001", "HrmA": "HrmA", "HrmA_middle": "HrmA", "HrmA_N": "HrmA",
            "NPHP3_N": "NPHP3_N", "Condensation": "Condensation", "PP-binding": "PP-binding",
            "AMP-binding": "AMP-binding", "SPOK2_N": "SPOK2_N"}
TIER1 = ["APH", "DUF3435", "DUF3723", "DUF3723_N"]
LABELS = TIER1 + ["PF11001", "HrmA", "NPHP3_N", "Condensation", "PP-binding", "AMP-binding", "SPOK2_N"]
BINS = ["core", "soft_core", "shell", "cloud", "singleton", "nonrep_only", "outgroup_only"]
NPERM = 1000
SEED = 20261017


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
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def load_run(tag):
    run, pg, fp, pi = RUNS[tag]
    fam_pf = {}       # family -> list of Pfam names (all, cut_ga)
    fam_lab = {}      # family -> set of labels
    with gzip.open(fp, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["pfam_ids"]:
                continue
            accs = r["pfam_ids"].split("|")
            names = r["pfam_names"].split("|")
            fam_pf[r["accession"]] = names
            labs = {ACC2LAB[a.split(".")[0]] for a in accs if a.split(".")[0] in ACC2LAB}
            if labs:
                fam_lab[r["accession"]] = labs
    freq = {}
    with open(pg / "frequency_table.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            freq[r["family"]] = r["bin"]
    islands = []
    with open(pg / "significant_islands.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            islands.append(dict(n_strains=int(r["n_strains"]), island_size=int(r["island_size"]),
                                example_strain=r["example_strain"],
                                members=r["member_families"].split(",")))
    pisl = {}
    with gzip.open(pi, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pisl[int(r["island_idx"])] = r
    return dict(run=run, fam_pf=fam_pf, fam_lab=fam_lab, freq=freq, islands=islands, pisl=pisl)


def pair_counts(M, L, nlab):
    """M: islands x families (csr, 0/1). L: families x labels (0/1 int).
    Returns obs_any[nlab,nlab], obs_distinct[nlab,nlab]; diagonal of distinct = islands
    with >= 2 different families carrying the label."""
    C = (M @ L)                              # islands x labels: families with label
    C = np.asarray(C)
    pres = C > 0
    any_ = pres.T.astype(np.int64) @ pres.astype(np.int64)
    dist = np.zeros((nlab, nlab), dtype=np.int64)
    for a in range(nlab):
        dist[a, a] = int((C[:, a] >= 2).sum())
        for b in range(a + 1, nlab):
            both = np.asarray(M @ (L[:, a] & L[:, b])).ravel()
            ok = pres[:, a] & pres[:, b] & ~((C[:, a] == 1) & (C[:, b] == 1) & (both == 1))
            dist[a, b] = dist[b, a] = int(ok.sum())
    return any_, dist


def family_pair_counts(G, L):
    """G: families x families (csr, 0/1, diagonal = 1), co-membership in >= 1 island.
    Returns number of distinct family pairs (f != g) with f carrying label a and g label b
    that share >= 1 island; diagonal = unordered pairs of two different families with label a."""
    GL = np.asarray(G @ L)
    X = L.T @ GL                      # sum_f,g La_f G_fg Lb_g (includes f == g)
    same = L.T @ L                    # f == g terms (G_ff = 1)
    out = X - same
    d = np.diag(out).copy() // 2
    np.fill_diagonal(out, d)
    return out


def main():
    rng = np.random.default_rng(SEED)
    co_rows, content_summary, largest, comp_rows, other_rows, famcount_rows = [], [], [], [], [], []
    data = {}
    for tag in RUNS:
        R = load_run(tag)
        data[tag] = R
        isl = R["islands"]
        fams = sorted({f for i in isl for f in i["members"]})
        fidx = {f: k for k, f in enumerate(fams)}
        rows, cols = [], []
        for i, it in enumerate(isl):
            for f in set(it["members"]):
                rows.append(i)
                cols.append(fidx[f])
        M = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(len(isl), len(fams)))
        nlab = len(LABELS)
        L = np.zeros((len(fams), nlab), dtype=np.int64)
        for f, k in fidx.items():
            for lab in R["fam_lab"].get(f, ()):
                L[k, LABELS.index(lab)] = 1
        G = (M.T @ M)
        G.data[:] = 1
        G = G.tocsr()
        obs_any, obs_dist = pair_counts(M, L, nlab)
        obs_fp = family_pair_counts(G, L)
        perm_any = np.zeros((NPERM, nlab, nlab), dtype=np.int64)
        perm_dist = np.zeros((NPERM, nlab, nlab), dtype=np.int64)
        perm_fp = np.zeros((NPERM, nlab, nlab), dtype=np.int64)
        for p in range(NPERM):
            Lp = L[rng.permutation(len(fams))]
            perm_any[p], perm_dist[p] = pair_counts(M, Lp, nlab)
            perm_fp[p] = family_pair_counts(G, Lp)
        nfam_lab = L.sum(axis=0)
        for a in range(nlab):
            for b in range(a, nlab):
                for kind, obs, perm in (("any", obs_any, perm_any), ("distinct_families", obs_dist, perm_dist),
                                        ("family_pairs", obs_fp, perm_fp)):
                    if kind == "any" and a == b:
                        continue
                    o = int(obs[a, b])
                    pv = perm[:, a, b]
                    e = float(pv.mean())
                    co_rows.append({"run": R["run"], "tag": tag, "label_a": LABELS[a], "label_b": LABELS[b],
                                    "count_type": kind,
                                    "n_island_families_a": int(nfam_lab[a]), "n_island_families_b": int(nfam_lab[b]),
                                    "observed_islands": o, "perm_mean": round(e, 3),
                                    "perm_q2.5": float(np.percentile(pv, 2.5)),
                                    "perm_q97.5": float(np.percentile(pv, 97.5)),
                                    "obs_over_exp": (round(o / e, 3) if e > 0 else "NA"),
                                    "p_upper": round((1 + int((pv >= o).sum())) / (NPERM + 1), 5),
                                    "p_lower": round((1 + int((pv <= o).sum())) / (NPERM + 1), 5)})
        print(f"{tag}: islands={len(isl)} island_families={len(fams)} labelled={int((L.sum(1) > 0).sum())}",
              file=sys.stderr)

        # ---- task 3: island content ----
        fam_lab, fam_pf, freq, pisl = R["fam_lab"], R["fam_pf"], R["freq"], R["pisl"]
        crows = []
        ntier = Counter()
        agg_bins, agg_lab, agg_other = Counter(), Counter(), Counter()
        for i, it in enumerate(isl):
            mem = list(dict.fromkeys(it["members"]))
            labs = Counter(l for f in mem for l in fam_lab.get(f, ()))
            t1 = [l for l in TIER1 if labs[l]]
            if not t1:
                continue
            ntier[len(t1)] += 1
            bins = Counter(freq.get(f, "NA") for f in mem)
            other = Counter(n for f in mem for n in dict.fromkeys(fam_pf.get(f, []))
                            if NAME2LAB.get(n) is None)
            pr = pisl.get(i, {})
            row = dict(island_idx=i, n_strains=it["n_strains"], example_strain=it["example_strain"],
                       island_size=it["island_size"], n_member_families=len(mem),
                       n_tier1_labels=len(t1), tier1_labels=",".join(t1),
                       labels_families=",".join(f"{l}:{labs[l]}" for l in LABELS if labs[l]),
                       n_unannotated_families=sum(1 for f in mem if f not in fam_pf),
                       classes=",".join(f"{b}:{bins[b]}" for b in BINS + ["NA"] if bins[b]),
                       top10_other_pfam=",".join(f"{n}:{c}" for n, c in other.most_common(10)),
                       duf3435_rep_family=("yes" if labs["DUF3435"] else "no"),
                       n_captain_families_SIO=pr.get("n_captain_families", "NA"))
            if tag == "af":
                row.update(n_starship_families_S7=pr.get("n_starship_families", "NA"),
                           starship_elements_family_level=pr.get("starship_elements_family_level", "") or "none",
                           starship_elements_coordinate=pr.get("starship_elements_coordinate", "") or "none")
            crows.append(row)
            for b in BINS + ["NA"]:
                agg_bins[b] += bins[b]
            for l in LABELS:
                if labs[l]:
                    agg_lab[l] += 1
            for n in other:
                agg_other[n] += 1
        write_tsv(T / f"island_content.{tag}.tsv.gz", crows)
        n_t1 = len(crows)
        for k in (1, 2, 3, 4):
            content_summary.append(dict(run=R["run"], tag=tag, n_tier1_labels=k, n_islands=ntier[k],
                                        n_islands_total=len(isl), n_islands_with_tier1=n_t1))
        multi = [("any_strains", r) for r in sorted((r for r in crows if r["n_tier1_labels"] >= 2),
                                                     key=lambda r: (-r["n_member_families"], -r["n_strains"]))[:10]]
        multi += [("ge2_strains", r) for r in sorted((r for r in crows if r["n_tier1_labels"] >= 2 and r["n_strains"] >= 2),
                                                     key=lambda r: (-r["n_member_families"], -r["n_strains"]))[:10]]
        for sel, r in multi:
            largest.append(dict(run=R["run"], tag=tag, selection=sel, **{k: r[k] for k in
                                ("island_idx", "n_strains", "n_member_families", "tier1_labels", "labels_families",
                                 "classes", "n_captain_families_SIO", "top10_other_pfam")},
                                **({"starship_elements_family_level": r["starship_elements_family_level"],
                                    "starship_elements_coordinate": r["starship_elements_coordinate"]}
                                   if tag == "af" else {})))
        nfam_tot = sum(agg_bins.values())
        comp = dict(run=R["run"], tag=tag, n_tier1_islands=n_t1,
                    single_strain=sum(r["n_strains"] == 1 for r in crows),
                    median_strains=float(np.median([r["n_strains"] for r in crows])) if crows else "",
                    max_strains=max((r["n_strains"] for r in crows), default=""),
                    median_member_families=float(np.median([r["n_member_families"] for r in crows])) if crows else "",
                    max_member_families=max((r["n_member_families"] for r in crows), default=""),
                    with_captain_SIO=sum(1 for r in crows if str(r["n_captain_families_SIO"]) not in ("0", "NA")),
                    with_duf3435_rep=sum(r["duf3435_rep_family"] == "yes" for r in crows),
                    family_slots=nfam_tot)
        for b in BINS + ["NA"]:
            comp[f"frac_family_slots_{b}"] = round(agg_bins[b] / nfam_tot, 4) if nfam_tot else ""
        for l in LABELS:
            comp[f"frac_islands_with_{l}"] = round(agg_lab[l] / n_t1, 4) if n_t1 else ""
        if tag == "af":
            comp["with_starship_family_S7"] = sum(1 for r in crows if str(r["n_starship_families_S7"]) not in ("0", "NA"))
            comp["with_starship_coordinate_S7"] = sum(1 for r in crows if r["starship_elements_coordinate"] != "none")
        comp_rows.append(comp)
        for n, c in agg_other.most_common(25):
            other_rows.append(dict(run=R["run"], tag=tag, pfam=n, n_tier1_islands_with=c,
                                   frac=round(c / n_t1, 4)))
        # family counts per label
        isl_fams = set(fams)
        for l in LABELS:
            fl = [f for f, s in fam_lab.items() if l in s]
            bc = Counter(freq.get(f, "NA") for f in fl)
            famcount_rows.append(dict(run=R["run"], tag=tag, label=l, n_families=len(fl),
                                      n_island_resident=sum(f in isl_fams for f in fl),
                                      **{b: bc[b] for b in BINS}, NA=bc["NA"]))
    add_bh(co_rows, ("tag", "count_type"), "p_upper", "q_upper_bh")
    write_tsv(T / "island_pair_cooccurrence.tsv", co_rows)
    write_tsv(T / "island_tier1_count_summary.tsv", content_summary)
    write_tsv(T / "island_largest_multidomain.tsv", largest)
    write_tsv(T / "island_composition_summary.tsv", comp_rows)
    write_tsv(T / "island_other_pfam_top.tsv", other_rows)
    write_tsv(T / "family_label_counts.tsv", famcount_rows)

    # ---- task 4: cross-genus best-hit pairs ----
    af_lab = data["af"]["fam_lab"]
    pairs = defaultdict(list)
    with gzip.open(XS / "family_pairs.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pairs[r["cocci_run"]].append(r)
    xrows = []
    for tag in ("v5", "v2"):
        run = RUNS[tag][0]
        c_lab = data[tag]["fam_lab"]
        for l in LABELS:
            af_f = {f for f, s in af_lab.items() if l in s}
            for subset in ("best_hit", "reciprocal"):
                ps = [r for r in pairs[run] if r["af_family"] in af_f and (subset == "best_hit" or r["reciprocal"] == "True")]
                shared = [r for r in ps if l in c_lab.get(r["cocci_family"], ())]
                def rng_(k, rs=shared):
                    v = [float(r[k]) for r in rs]
                    return f"{min(v):.3f}-{max(v):.3f} (median {np.median(v):.3f})" if v else ""
                xrows.append(dict(cocci_run=run, tag=tag, label=l, pair_type=subset, n_af_families_with_label=len(af_f),
                                  n_af_with_partner=len(ps), n_partner_has_label=len(shared),
                                  n_distinct_cocci_partners_with_label=len({r["cocci_family"] for r in shared}),
                                  fident_range=rng_("fident"), qcov_range=rng_("qcov"), tcov_range=rng_("tcov"),
                                  af_bins=",".join(f"{k}:{v}" for k, v in Counter(r["af_bin"] for r in shared).most_common()),
                                  cocci_bins=",".join(f"{k}:{v}" for k, v in Counter(r["cocci_bin"] for r in shared).most_common()),
                                  af_island_true=sum(r["af_island"] == "True" for r in shared),
                                  cocci_island_true=sum(r["cocci_island"] == "True" for r in shared),
                                  both_island=sum(r["af_island"] == "True" and r["cocci_island"] == "True" for r in shared)))
    write_tsv(T / "crossgenus_pairs_by_domain.tsv", xrows)
    cg = []
    with open(XS / "cross_genus_island_families.tsv", newline="") as fh:
        cgr = list(csv.DictReader(fh, delimiter="\t"))
    for l in LABELS:
        sel = [r for r in cgr if l in af_lab.get(r["af_family"], ())]
        cg.append(dict(label=l, n_cross_genus_island_af_families=len(sel), n_total_rows=len(cgr),
                       v5_partner_has_label=sum(l in data["v5"]["fam_lab"].get(r["v5_family"], ()) for r in sel),
                       v2_partner_has_label=sum(l in data["v2"]["fam_lab"].get(r["v2_family"], ()) for r in sel),
                       af_families=";".join(r["af_family"] for r in sel),
                       af_bins=";".join(r["af_bin"] for r in sel),
                       af_starship_elements=";".join(r["af_starship_elements"] or "-" for r in sel),
                       af_captain=";".join(r["af_captain"] for r in sel)))
    write_tsv(T / "crossgenus_island_families_by_domain.tsv", cg)


if __name__ == "__main__":
    main()
