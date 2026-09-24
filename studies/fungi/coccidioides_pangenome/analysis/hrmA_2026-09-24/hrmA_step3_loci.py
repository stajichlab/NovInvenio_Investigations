#!/usr/bin/env python3
"""Step 3: group HrmA families into genomic loci by shared conserved (core/soft_core)
flanking families; score per-strain locus occupancy with an empty-site check;
species tests at the locus level. Read-only on pipeline outputs."""
import io, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/hrmA_2026-09-24"
W = 12          # gene window around marker genes
MIN_FRAC = 0.5  # marker must flank >= this fraction of a family's anchors
MIN_SHARED = 2  # families sharing >= this many markers are one locus


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


census = pd.read_csv(A / "hrmA_family_census.tsv", sep="\t")
lab = dict(zip(census.family, census.label))
freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")
nb = pd.read_csv(A / "hrmA_neighbourhoods.tsv.gz", sep="\t")
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
sp = dict(zip(ss.Short, ss.Species.map({"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"})))
nsp = Counter(sp.values())

# markers per family
nb_off = nb[nb.offset != 0]
nanch = nb[nb.offset == 0].groupby("anchor_family").anchor.nunique()
markers = {}
for fam, d in nb_off.groupby("anchor_family"):
    c = d.groupby("family").anchor.nunique() / nanch[fam]
    c = c[c >= MIN_FRAC]
    c = c[[freq.bin.get(f, "") in ("core", "soft_core") for f in c.index]]
    markers[fam] = c.sort_values(ascending=False)

# union-find on shared markers
fams = list(census.family)
parent = {f: f for f in fams}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for i, a in enumerate(fams):
    for b in fams[i + 1:]:
        if a in markers and b in markers and len(markers[a]) and len(markers[b]) and \
                len(set(markers[a].index) & set(markers[b].index)) >= min(MIN_SHARED, len(markers[a]), len(markers[b])):
            parent[find(a)] = find(b)
groups = defaultdict(list)
for f in fams:
    groups[find(f)].append(f)
loci = sorted(groups.values(), key=lambda g: -census.set_index("family").loc[g, "n_members"].sum())
locus_of = {}
for i, g in enumerate(loci):
    for f in g:
        locus_of[f] = f"L{i+1:02d}"
census["locus"] = census.family.map(locus_of)
census.to_csv(A / "hrmA_family_census.tsv", sep="\t", index=False)

# --- gene table with families for all strains
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
hits = pd.read_csv(A / "hrmA_gene_locations.tsv", sep="\t")
hitset = set(hits.member)
gp["is_hrmA"] = gp.member.isin(hitset)
rp = pd.read_csv(P / "rescue_positions.tsv", sep="\t")
rp = rp[rp.family.isin(set(fams))]

rows, occ = [], []
for li, g in enumerate(loci):
    L = f"L{li+1:02d}"
    mk = Counter()
    for f in g:
        if f in markers:
            for m, v in markers[f].items():
                mk[m] = max(mk[m], v)
    mlist = [m for m, _ in mk.most_common(6)]
    sub = gp[gp.family.isin(mlist)]
    for s in ss.Short:
        ms = sub[sub.Short == s]
        state, fam_found = "no_marker", ""
        if len(ms):
            state = "unresolved"
            for r in ms.itertuples():
                lo, hi = max(0, r.idx - W), min(r.n_on - 1, r.idx + W)
                win = gp.iloc[r.Index - (r.idx - lo): r.Index + (hi - r.idx) + 1]
                h = win[win.is_hrmA]
                if len(h):
                    state = "occupied"
                    fam_found = ",".join(sorted(set(h.family.map(lambda x: lab.get(x, x)))))
                    break
                if r.idx - W >= 0 and r.idx + W <= r.n_on - 1:
                    state = "empty_site"
            if state != "occupied":
                # rescue-only call inside the marker window span?
                rr = rp[(rp.Short == s)]
                for r in ms.itertuples():
                    lo, hi = max(0, r.idx - W), min(r.n_on - 1, r.idx + W)
                    win = gp.iloc[r.Index - (r.idx - lo): r.Index + (hi - r.idx) + 1]
                    x = rr[(rr.contig == r.contig) & (rr.start >= win.start.min()) & (rr.start <= win.end.max())]
                    if len(x):
                        state = state + "+rescue"
                        fam_found = ",".join(sorted(set(x.family.map(lambda y: lab.get(y, y)))))
                        break
        occ.append(dict(locus=L, Short=s, species=sp[s], state=state, families=fam_found))
    o = pd.DataFrame([x for x in occ if x["locus"] == L])
    d = dict(locus=L, families=",".join(lab[f] for f in g), n_markers=len(mlist),
             markers=",".join(mlist))
    for spc in ("Ci", "Cp"):
        c = Counter(o[o.species == spc].state)
        for k in ("occupied", "empty_site", "unresolved", "no_marker", "empty_site+rescue", "unresolved+rescue", "no_marker+rescue"):
            d[f"{spc}_{k}"] = c.get(k, 0)
    # Fisher on resolved strains only: occupied vs empty_site (pure)
    a, b = d["Ci_occupied"], d["Ci_empty_site"]
    c_, e = d["Cp_occupied"], d["Cp_empty_site"]
    d["Ci_occ_frac_resolved"] = round(a / (a + b), 3) if a + b else None
    d["Cp_occ_frac_resolved"] = round(c_ / (c_ + e), 3) if c_ + e else None
    d["fisher_p_resolved"] = fisher_exact([[a, b], [c_, e]])[1] if (a + b) and (c_ + e) else None
    rows.append(d)
lt = pd.DataFrame(rows)
ok = lt.fisher_p_resolved.notna()
lt.loc[ok, "bh_q_resolved"] = multipletests(lt.loc[ok, "fisher_p_resolved"], method="fdr_bh")[1]
lt.to_csv(A / "hrmA_loci.tsv", sep="\t", index=False)
pd.DataFrame(occ).to_csv(A / "hrmA_locus_occupancy.tsv.gz", sep="\t", index=False, compression="gzip")
print(lt.drop(columns=["markers"]).to_string(), file=sys.stderr)
