#!/usr/bin/env python3
"""Step 5: sensitivity re-scoring of loci. Occupancy counts only genes of the
locus's own families (step 2 counted any HET/NLR gene). L001 is split into two
physical sites (L001a NB-ARC fragments; L001b DUF7104-NPHP3_N), because step 2
joined them through shared markers ~20 genes apart. Same W / marker rules."""
import io, subprocess, sys
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/het_nlr_2026-09-25"
W, MIN_FRAC, NMARK, MIN_RES = 12, 0.5, 6, 10

census = pd.read_csv(A / "het_family_census.tsv", sep="\t")
lab = dict(zip(census.family, census.label)); rl = {v: k for k, v in lab.items()}
freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")
nb = pd.read_csv(A / "het_neighbourhoods.tsv.gz", sep="\t")
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
sp = dict(zip(ss.Short, ss.Species.map({"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"})))
lt = pd.read_csv(A / "het_loci.tsv", sep="\t")

sites = {}
for r in lt[lt.n_markers > 0].itertuples():
    if r.locus == "L001":
        f = r.families.split(",")
        cls = census.set_index("label").family_class
        sites["L001a"] = [x for x in f if cls[x] == "NLR_NOD"]
        sites["L001b"] = [x for x in f if cls[x] != "NLR_NOD"]
    else:
        sites[r.locus] = r.families.split(",")

nanch = nb[nb.offset == 0].groupby("anchor_family").anchor.nunique()
nbo = nb[nb.offset != 0]
def markers_for(labels):
    mk = Counter()
    for L in labels:
        f = rl[L]
        if f not in nanch:
            continue
        c = nbo[nbo.anchor_family == f].groupby("family").anchor.nunique() / nanch[f]
        for m, v in c[c >= MIN_FRAC].items():
            if freq.bin.get(m, "") in ("core", "soft_core"):
                mk[m] = max(mk[m], v)
    return [m for m, _ in mk.most_common(NMARK)]

gp = pd.read_csv(io.BytesIO(subprocess.run(["zstd", "-dc", str(P / "gene_positions.tsv.zst")], capture_output=True, check=True).stdout), sep="\t")
gp["member"] = gp.Short + "|" + gp.protein_id
m2f = {}
with open(P / "cluster/tier1_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t"); m2f[mem] = rep
gp["family"] = gp.member.map(m2f)
gp = gp.sort_values(["Short", "contig", "start"]).reset_index(drop=True)
gp["idx"] = gp.groupby(["Short", "contig"]).cumcount()
gp["n_on"] = gp.groupby(["Short", "contig"]).idx.transform("size")
het = set(pd.read_csv(A / "het_gene_locations.tsv.gz", sep="\t").member)
fam_arr, mem_arr = gp.family.values, gp.member.values

rows, arows, occ = [], [], []
for S, labels in sites.items():
    own = {rl[x] for x in labels}
    mlist = markers_for(labels)
    if not mlist:
        continue
    sub = gp[gp.family.isin(mlist)]
    by = {k: v for k, v in sub.groupby("Short")}
    st = []
    for s in ss.Short:
        ms = by.get(s); state, found = "no_marker", set()
        if ms is not None:
            state, full = "unresolved", False
            for r in ms.itertuples():
                lo, hi = max(0, r.idx - W), min(r.n_on - 1, r.idx + W)
                i0, i1 = r.Index - (r.idx - lo), r.Index + (hi - r.idx) + 1
                for k in range(i0, i1):
                    if fam_arr[k] in own and mem_arr[k] in het:
                        found.add(lab[fam_arr[k]])
                full = full or (r.idx - W >= 0 and r.idx + W <= r.n_on - 1)
            state = "occupied" if found else ("empty_site" if full else "unresolved")
        st.append(dict(site=S, Short=s, species=sp[s], state=state, families=",".join(sorted(found))))
    o = pd.DataFrame(st); occ.append(o)
    d = dict(site=S, families=",".join(labels), n_markers=len(mlist), markers=",".join(mlist))
    for spc in ("Ci", "Cp"):
        c = Counter(o[o.species == spc].state)
        for k in ("occupied", "empty_site", "unresolved", "no_marker"):
            d[f"{spc}_{k}"] = c.get(k, 0)
        d[f"{spc}_resolved"] = d[f"{spc}_occupied"] + d[f"{spc}_empty_site"]
        d[f"{spc}_occ_frac"] = round(d[f"{spc}_occupied"] / d[f"{spc}_resolved"], 3) if d[f"{spc}_resolved"] else None
    a, b, c_, e = d["Ci_occupied"], d["Ci_empty_site"], d["Cp_occupied"], d["Cp_empty_site"]
    d["fisher_p"] = fisher_exact([[a, b], [c_, e]])[1] if (a + b) and (c_ + e) else None
    rows.append(d)
    res = o[o.state.isin(["occupied", "empty_site"])]
    n = res.groupby("species").size()
    cnt = Counter((f, r.species) for r in res[res.state == "occupied"].itertuples() for f in r.families.split(","))
    for f in sorted({k[0] for k in cnt}):
        x, y = cnt[(f, "Ci")], cnt[(f, "Cp")]
        nci, ncp = int(n.get("Ci", 0)), int(n.get("Cp", 0))
        arows.append(dict(site=S, family_label=f, Ci_carriers=x, Ci_resolved=nci, Ci_frac=round(x / nci, 4) if nci else None,
                          Cp_carriers=y, Cp_resolved=ncp, Cp_frac=round(y / ncp, 4) if ncp else None,
                          fisher_p=fisher_exact([[x, nci - x], [y, ncp - y]])[1] if nci and ncp else None))
t = pd.DataFrame(rows)
ok = t.fisher_p.notna(); t.loc[ok, "bh_q"] = multipletests(t.loc[ok, "fisher_p"], method="fdr_bh")[1]
inr = lambda fr, nn: (nn >= MIN_RES) & (fr >= 0.05) & (fr <= 0.95)
t["poly_within"] = inr(t.Ci_occ_frac.fillna(-1), t.Ci_resolved) | inr(t.Cp_occ_frac.fillna(-1), t.Cp_resolved)
t["poly_between"] = t.bh_q < 0.05
t.to_csv(A / "het_sites_ownfamily.tsv", sep="\t", index=False)
al = pd.DataFrame(arows)
ok = al.fisher_p.notna(); al.loc[ok, "bh_q"] = multipletests(al.loc[ok, "fisher_p"], method="fdr_bh")[1]
al["poly_within"] = inr(al.Ci_frac.fillna(-1), al.Ci_resolved) | inr(al.Cp_frac.fillna(-1), al.Cp_resolved)
al["poly_between"] = al.bh_q < 0.05
al.to_csv(A / "het_sites_ownfamily_alleles.tsv", sep="\t", index=False)
pd.concat(occ).to_csv(A / "het_sites_ownfamily_occupancy.tsv.gz", sep="\t", index=False, compression="gzip")
print(t.drop(columns=["markers"]).to_string(), file=sys.stderr)
print(al[(al.Ci_carriers + al.Cp_carriers) >= 5].to_string(), file=sys.stderr)
