#!/usr/bin/env python3
"""Step 2: group HET/NLR families into genomic loci by shared conserved
(core/soft_core) flanking families; score per-strain locus occupancy with an
empty-site check; Ci vs Cp Fisher tests with BH. Same rules as
analysis/hrmA_2026-09-24/hrmA_step3_loci.py (W=12, MIN_FRAC=0.5, MIN_SHARED=2,
up to 6 markers). Read-only on pipeline outputs."""
import io, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/het_nlr_2026-09-25"
W, MIN_FRAC, MIN_SHARED, NMARK = 12, 0.5, 2, 6
MIN_RESOLVED = 10  # a species frequency is only called when >= this many strains are resolved


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


census = pd.read_csv(A / "het_family_census.tsv", sep="\t")
lab = dict(zip(census.family, census.label))
freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")
nb = pd.read_csv(A / "het_neighbourhoods.tsv.gz", sep="\t")
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
sp = dict(zip(ss.Short, ss.Species.map({"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"})))

nb_off = nb[nb.offset != 0]
nanch = nb[nb.offset == 0].groupby("anchor_family").anchor.nunique()
markers = {}
for fam, d in nb_off.groupby("anchor_family"):
    c = d.groupby("family").anchor.nunique() / nanch[fam]
    c = c[c >= MIN_FRAC]
    c = c[[freq.bin.get(f, "") in ("core", "soft_core") for f in c.index]]
    markers[fam] = c.sort_values(ascending=False)

fams = list(census.family)
parent = {f: f for f in fams}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


mk_sets = {f: set(markers[f].index) for f in fams if f in markers and len(markers[f])}
for i, a in enumerate(fams):
    if a not in mk_sets:
        continue
    for b in fams[i + 1:]:
        if b in mk_sets and len(mk_sets[a] & mk_sets[b]) >= min(MIN_SHARED, len(mk_sets[a]), len(mk_sets[b])):
            parent[find(a)] = find(b)
groups = defaultdict(list)
for f in fams:
    groups[find(f)].append(f)
nmem = census.set_index("family").n_members
loci = sorted(groups.values(), key=lambda g: -nmem.loc[g].sum())
locus_of = {f: f"L{i+1:03d}" for i, g in enumerate(loci) for f in g}
census["locus"] = census.family.map(locus_of)
census["n_markers_family"] = census.family.map(lambda f: len(mk_sets.get(f, ())))
census.to_csv(A / "het_family_census.tsv", sep="\t", index=False)

print("loading genes", file=sys.stderr)
gp = zread(P / "gene_positions.tsv.zst")
gp["member"] = gp.Short + "|" + gp.protein_id
m2f = {}
with open(P / "cluster/tier1_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        m2f[mem] = rep
gp["family"] = gp.member.map(m2f)
gp = gp.sort_values(["Short", "contig", "start"]).reset_index(drop=True)
gp["idx"] = gp.groupby(["Short", "contig"]).cumcount()
gp["n_on"] = gp.groupby(["Short", "contig"]).idx.transform("size")
hits = pd.read_csv(A / "het_gene_locations.tsv.gz", sep="\t")
gp["is_het"] = gp.member.isin(set(hits.member))
gp["lab"] = gp.family.map(lab)
is_het = gp.is_het.values
labs = gp.lab.values
fam_arr = gp.family.values
rp = pd.read_csv(P / "rescue_positions.tsv", sep="\t")
rp = rp[rp.family.isin(set(fams))]
rp_by = {k: v for k, v in rp.groupby(["Short", "contig"])}

all_markers = set()
locus_markers = {}
for li, g in enumerate(loci):
    mk = Counter()
    for f in g:
        if f in markers:
            for m, v in markers[f].items():
                mk[m] = max(mk[m], v)
    locus_markers[f"L{li+1:03d}"] = [m for m, _ in mk.most_common(NMARK)]
    all_markers |= set(locus_markers[f"L{li+1:03d}"])
msub = gp[gp.family.isin(all_markers)]
by_fam = {k: v for k, v in msub.groupby("family")}

rows, occ = [], []
for li, g in enumerate(loci):
    L = f"L{li+1:03d}"
    mlist = locus_markers[L]
    gset = set(g)
    if not mlist:
        rows.append(dict(locus=L, families=",".join(lab[f] for f in g), n_markers=0, markers=""))
        continue
    sub = pd.concat([by_fam[m] for m in mlist if m in by_fam])
    sub_by = {k: v for k, v in sub.groupby("Short")}
    for s in ss.Short:
        ms = sub_by.get(s)
        state, found, own, genes = "no_marker", "", False, []
        if ms is not None:
            state = "unresolved"
            fset, full = set(), False
            genes = []
            for r in ms.itertuples():
                lo, hi = max(0, r.idx - W), min(r.n_on - 1, r.idx + W)
                i0, i1 = r.Index - (r.idx - lo), r.Index + (hi - r.idx) + 1
                w = np.where(is_het[i0:i1])[0]
                for k in w:
                    fset.add(str(labs[i0 + k]))
                    genes.append(gp.member.iat[i0 + k])
                    own = own or fam_arr[i0 + k] in gset
                if r.idx - W >= 0 and r.idx + W <= r.n_on - 1:
                    full = True
            if fset:
                state = "occupied"
                found = ",".join(sorted(fset))
            elif full:
                state = "empty_site"
            if state != "occupied":
                for r in ms.itertuples():
                    x = rp_by.get((s, r.contig))
                    if x is None:
                        continue
                    lo, hi = max(0, r.idx - W), min(r.n_on - 1, r.idx + W)
                    i0, i1 = r.Index - (r.idx - lo), r.Index + (hi - r.idx)
                    x = x[(x.start >= gp.start.iat[i0]) & (x.start <= gp.end.iat[i1])]
                    if len(x):
                        state += "+rescue"
                        found = ",".join(sorted(set(x.family.map(lab))))
                        break
        occ.append(dict(locus=L, Short=s, species=sp[s], state=state, families=found, own_family=own,
                        genes=",".join(sorted(set(genes))) if ms is not None else ""))
    o = pd.DataFrame([x for x in occ if x["locus"] == L])
    d = dict(locus=L, families=",".join(lab[f] for f in g), n_markers=len(mlist), markers=",".join(mlist))
    for spc in ("Ci", "Cp"):
        oo = o[o.species == spc]
        c = Counter(oo.state)
        for k in ("occupied", "empty_site", "unresolved", "no_marker", "empty_site+rescue", "unresolved+rescue", "no_marker+rescue"):
            d[f"{spc}_{k}"] = c.get(k, 0)
        d[f"{spc}_occupied_by_other_family_only"] = int(((oo.state == "occupied") & ~oo.own_family).sum())
    a, b, c_, e = d["Ci_occupied"], d["Ci_empty_site"], d["Cp_occupied"], d["Cp_empty_site"]
    d["Ci_resolved"], d["Cp_resolved"] = a + b, c_ + e
    d["Ci_occ_frac_resolved"] = round(a / (a + b), 3) if a + b else None
    d["Cp_occ_frac_resolved"] = round(c_ / (c_ + e), 3) if c_ + e else None
    d["fisher_p_resolved"] = fisher_exact([[a, b], [c_, e]])[1] if (a + b) and (c_ + e) else None
    rows.append(d)
lt = pd.DataFrame(rows)
ok = lt.fisher_p_resolved.notna()
lt.loc[ok, "bh_q_resolved"] = multipletests(lt.loc[ok, "fisher_p_resolved"], method="fdr_bh")[1]


def within(fr, n):
    return (n >= MIN_RESOLVED) & (fr >= 0.05) & (fr <= 0.95)


lt["poly_within"] = within(lt.Ci_occ_frac_resolved.fillna(-1), lt.Ci_resolved.fillna(0)) | \
    within(lt.Cp_occ_frac_resolved.fillna(-1), lt.Cp_resolved.fillna(0))
lt["poly_between"] = lt.bh_q_resolved < 0.05
# family composition summary
cz = census.set_index("family")
lt["family_classes"] = [";".join(f"{k}:{v}" for k, v in Counter(cz.family_class[f] for f in g).most_common()) for g in loci]
lt["architectures"] = [";".join(sorted({cz.modal_architecture[f] for f in g})) for g in loci]
lt["n_members_total"] = [int(cz.n_members[g].sum()) for g in loci]
lt.to_csv(A / "het_loci.tsv", sep="\t", index=False)
pd.DataFrame(occ).to_csv(A / "het_locus_occupancy.tsv.gz", sep="\t", index=False, compression="gzip")
print(lt[lt.n_markers > 0].drop(columns=["markers"]).head(60).to_string(), file=sys.stderr)

# ---------------- allele (family-at-locus) table
occ = pd.DataFrame(occ)
arows = []
for L, d in occ.groupby("locus"):
    res = d[d.state.isin(["occupied", "empty_site"])]
    nres = res.groupby("species").size()
    cnt = Counter()
    for r in res[res.state == "occupied"].itertuples():
        for f in r.families.split(","):
            cnt[(f, r.species)] += 1
    for f in sorted({k[0] for k in cnt}):
        a, c = cnt[(f, "Ci")], cnt[(f, "Cp")]
        nci, ncp = int(nres.get("Ci", 0)), int(nres.get("Cp", 0))
        pv = fisher_exact([[a, nci - a], [c, ncp - c]])[1] if nci and ncp else None
        arows.append(dict(locus=L, family_label=f, Ci_carriers=a, Ci_resolved=nci,
                          Ci_frac=round(a / nci, 4) if nci else None, Cp_carriers=c, Cp_resolved=ncp,
                          Cp_frac=round(c / ncp, 4) if ncp else None, fisher_p=pv))
al = pd.DataFrame(arows)
ok = al.fisher_p.notna()
al.loc[ok, "bh_q"] = multipletests(al.loc[ok, "fisher_p"], method="fdr_bh")[1]
al["poly_within"] = within(al.Ci_frac.fillna(-1), al.Ci_resolved) | within(al.Cp_frac.fillna(-1), al.Cp_resolved)
al["poly_between"] = al.bh_q < 0.05
cz2 = census.set_index("label")
al["family_class"] = al.family_label.map(cz2.family_class)
al["modal_architecture"] = al.family_label.map(cz2.modal_architecture)
al["home_locus"] = al.family_label.map(cz2.locus)
al.to_csv(A / "het_locus_alleles.tsv", sep="\t", index=False)

# allele diversity per locus and species: distinct occupant combinations, Simpson 1-sum(p^2)
dv = []
for (L, spc), d in occ[occ.state == "occupied"].groupby(["locus", "species"]):
    c = d.families.value_counts()
    p = c / c.sum()
    dv.append(dict(locus=L, species=spc, n_occupied=int(c.sum()), n_combinations=len(c),
                   top_combination=c.index[0], top_frac=round(p.iloc[0], 3), simpson=round(1 - (p ** 2).sum(), 3)))
dv = pd.DataFrame(dv)
dv.to_csv(A / "het_locus_allele_diversity.tsv", sep="\t", index=False)
# Ci vs Cp difference in combination spectrum (chi-square on combos with >=5 strains total)
from scipy.stats import chi2_contingency
cs = []
for L, d in occ[occ.state == "occupied"].groupby("locus"):
    t = pd.crosstab(d.families, d.species)
    t = t[t.sum(axis=1) >= 5]
    if t.shape[0] >= 2 and t.shape[1] == 2:
        chi, p, dof, _ = chi2_contingency(t)
        cs.append(dict(locus=L, n_combinations_tested=t.shape[0], chi2=round(chi, 2), dof=dof, p=p))
cs = pd.DataFrame(cs)
if len(cs):
    cs["bh_q"] = multipletests(cs.p, method="fdr_bh")[1]
cs.to_csv(A / "het_locus_allele_spectrum_tests.tsv", sep="\t", index=False)
print(al.head(60).to_string(), file=sys.stderr)
