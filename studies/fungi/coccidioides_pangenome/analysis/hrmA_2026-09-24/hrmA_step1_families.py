#!/usr/bin/env python3
"""Step 1: HrmA-domain family census, presence/absence per species, locations,
neighbourhoods. Read-only on pipeline outputs; writes TSVs into this directory.
Run with the NII pixi python."""
import gzip, io, subprocess, sys
from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/hrmA_2026-09-24"
HRMA_MODELS = ["HrmA", "HrmA_N", "HrmA_middle"]
FLANK = 10


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


def parse_domtbl(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        rows.append(dict(target=f[0], tlen=int(f[2]), model=f[3], full_E=float(f[6]),
                         full_score=float(f[7]), dom_iE=float(f[12]), dom_score=float(f[13]),
                         hmm_from=int(f[15]), hmm_to=int(f[16]), env_from=int(f[19]), env_to=int(f[20])))
    return pd.DataFrame(rows)


# --- samples
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
species = dict(zip(ss.Short, ss.Species))
strains = list(ss.Short)
sp_short = {"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"}
n_sp = Counter(species[s] for s in strains)

# --- hmm hits
reps_dom = parse_domtbl(A / "reps.domtblout")
all_dom = parse_domtbl(A / "allprot.domtblout")
all_dom["Short"] = all_dom.target.str.split("|").str[0]
all_dom["protein_id"] = all_dom.target.str.split("|", n=1).str[1]

# --- cluster map (member -> rep)
print("loading cluster map", file=sys.stderr)
m2f = {}
fam_members = defaultdict(list)
with open(P / "cluster/tier1_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        m2f[mem] = rep
        fam_members[rep].append(mem)
print("families:", len(fam_members), "members:", len(m2f), file=sys.stderr)

all_dom["family"] = all_dom.target.map(m2f)
unmapped = all_dom[all_dom.family.isna()].target.unique()
print("hit proteins not in cluster map:", len(unmapped), file=sys.stderr)
(A / "hits_not_in_cluster.txt").write_text("\n".join(unmapped) + "\n")

hrma_all = all_dom[all_dom.model.isin(HRMA_MODELS)]
hrma_reps = reps_dom[reps_dom.model.isin(HRMA_MODELS)]
fam_set = set(hrma_reps.target) | set(hrma_all.family.dropna())
print("HrmA families (rep hit or member hit):", len(fam_set), file=sys.stderr)

freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")

# per-protein model strings
prot_models = hrma_all.groupby("target").model.apply(lambda s: "+".join(sorted(set(s))))
prot_pf11001 = set(all_dom[all_dom.model == "AFUB_07903_YDR124W_hel"].target)

census = []
for fam in sorted(fam_set):
    mem = fam_members[fam]
    rep_models = "+".join(sorted(set(hrma_reps[hrma_reps.target == fam].model))) or "-"
    mem_hits = [m for m in mem if m in prot_models.index]
    mc = Counter(prot_models[m] for m in mem_hits)
    per_strain = Counter(m.split("|")[0] for m in mem)
    census.append(dict(family=fam, rep_models=rep_models,
                       rep_len=None,
                       n_members=len(mem), n_members_hrma_hit=len(mem_hits),
                       member_model_combos=";".join(f"{k}:{v}" for k, v in mc.most_common()),
                       n_members_PF11001=sum(1 for m in mem if m in prot_pf11001),
                       n_strains_with_member=len(per_strain),
                       max_copies_per_strain=max(per_strain.values()),
                       n_strains_multicopy=sum(1 for v in per_strain.values() if v > 1),
                       bin=freq.bin.get(fam, "NA"), strain_count=freq.strain_count.get(fam, None),
                       frequency=freq.frequency.get(fam, None)))
census = pd.DataFrame(census)

# rep length from fasta
lens = {}
cur = None
with open(P / "cluster/tier1_rep_seq.fasta") as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].split()[0]
            if cur in fam_set:
                lens[cur] = 0
            else:
                cur = None
        elif cur:
            lens[cur] += len(line.strip())
census["rep_len"] = census.family.map(lens)

# HrmA-hit proteins outside the HrmA family set? (all are in by construction) ;
# members of HrmA families that are not hit by any HrmA model:
census["frac_members_hit"] = (census.n_members_hrma_hit / census.n_members).round(3)

# --- presence matrix rows
print("reading presence matrix", file=sys.stderr)
pm_rows = {}
with open(P / "presence_matrix.rescued.tsv") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        fam = line[:line.index("\t")]
        if fam in fam_set:
            pm_rows[fam] = line.rstrip("\n").split("\t")[1:]
cols = header[1:]
pm = pd.DataFrame.from_dict(pm_rows, orient="index", columns=cols)
assert set(cols) == set(strains), "presence matrix strains differ from samplesheet"

pres = []
for fam in census.family:
    row = pm.loc[fam]
    d = dict(family=fam)
    for sp, lab in sp_short.items():
        sub = row[[s for s in cols if species[s] == sp]]
        c = Counter(sub)
        d[f"{lab}_present"] = c.get("present", 0)
        d[f"{lab}_genome_only"] = c.get("genome_only", 0)
        d[f"{lab}_absent"] = c.get("absent", 0)
        d[f"{lab}_n"] = len(sub)
    for mode, keys in (("any", ("present", "genome_only")), ("annot", ("present",))):
        ci = sum(d[f"Ci_{k}"] for k in keys)
        cp = sum(d[f"Cp_{k}"] for k in keys)
        _, p = fisher_exact([[ci, d["Ci_n"] - ci], [cp, d["Cp_n"] - cp]])
        d[f"Ci_frac_{mode}"] = round(ci / d["Ci_n"], 4)
        d[f"Cp_frac_{mode}"] = round(cp / d["Cp_n"], 4)
        d[f"fisher_p_{mode}"] = p
    pres.append(d)
pres = pd.DataFrame(pres)
for mode in ("any", "annot"):
    pres[f"bh_q_{mode}"] = multipletests(pres[f"fisher_p_{mode}"], method="fdr_bh")[1]
census = census.merge(pres, on="family")

# copy number per strain by species (annotated members)
cn = []
for fam in census.family:
    ps = Counter(m.split("|")[0] for m in fam_members[fam])
    for lab_sp, lab in sp_short.items():
        v = [ps[s] for s in strains if species[s] == lab_sp and ps[s] > 0]
        cn.append((fam, lab, len(v), (sum(v) / len(v)) if v else 0, max(v) if v else 0))
cn = pd.DataFrame(cn, columns=["family", "sp", "n", "mean", "max"])
for lab in ("Ci", "Cp"):
    s = cn[cn.sp == lab].set_index("family")
    census[f"{lab}_mean_copies"] = census.family.map(s["mean"]).round(2)
    census[f"{lab}_max_copies"] = census.family.map(s["max"])

census = census.sort_values(["rep_models", "n_members"], ascending=[True, False])
census.to_csv(A / "hrmA_family_census.tsv", sep="\t", index=False)
pm.loc[census.family].to_csv(A / "hrmA_presence_matrix.tsv.gz", sep="\t", compression="gzip")

# --- positions
print("loading gene positions", file=sys.stderr)
gp = zread(P / "gene_positions.tsv.zst")
gp["member"] = gp.Short + "|" + gp.protein_id
gp["family"] = gp.member.map(m2f)
clen = pd.read_csv(A / "contig_lengths.tsv.gz", sep="\t", header=None, names=["Short", "contig", "contig_len"])
gp = gp.merge(clen, on=["Short", "contig"], how="left")
print("genes without contig length:", gp.contig_len.isna().sum(), file=sys.stderr)
gp["dist_end"] = pd.concat([gp.start - 1, gp.contig_len - gp.end], axis=1).min(axis=1)
gp = gp.sort_values(["Short", "contig", "start"]).reset_index(drop=True)
gp["idx"] = gp.groupby(["Short", "contig"]).cumcount()
gp["n_on_contig"] = gp.groupby(["Short", "contig"]).idx.transform("size")

# background: all genes, and genes by bin
gp["bin"] = gp.family.map(freq.bin)
bg = []
for label, sub in [("all_genes", gp)] + [(f"bin_{b}", gp[gp.bin == b]) for b in ["core", "soft_core", "shell", "cloud", "singleton"]]:
    bg.append(dict(set=label, n=len(sub), median_dist_end=sub.dist_end.median(),
                   frac_le20kb=round((sub.dist_end <= 20000).mean(), 4),
                   frac_le50kb=round((sub.dist_end <= 50000).mean(), 4),
                   median_contig_len=sub.contig_len.median(),
                   frac_contig_lt100kb=round((sub.contig_len < 100000).mean(), 4)))

hit_members = set(hrma_all.target)
hg = gp[gp.member.isin(hit_members)].copy()
hg["hrmA_models"] = hg.member.map(prot_models)
hg["species"] = hg.Short.map(species).map(sp_short)
hg["pf11001_same_protein"] = hg.member.isin(prot_pf11001)
hg.to_csv(A / "hrmA_gene_locations.tsv", sep="\t", index=False)
print("HrmA-hit proteins:", len(hit_members), "with positions:", len(hg), file=sys.stderr)
for label, sub in [("HrmA_hit_genes_all", hg)] + [(f"HrmA_hit_{m}", hg[hg.hrmA_models.str.contains(m + r"(?:\+|$)", regex=True)]) for m in HRMA_MODELS]:
    bg.append(dict(set=label, n=len(sub), median_dist_end=sub.dist_end.median(),
                   frac_le20kb=round((sub.dist_end <= 20000).mean(), 4),
                   frac_le50kb=round((sub.dist_end <= 50000).mean(), 4),
                   median_contig_len=sub.contig_len.median(),
                   frac_contig_lt100kb=round((sub.contig_len < 100000).mean(), 4)))
pd.DataFrame(bg).to_csv(A / "hrmA_location_vs_background.tsv", sep="\t", index=False)

# --- neighbourhoods
print("neighbourhoods", file=sys.stderr)
gpi = gp.set_index(["Short", "contig", "idx"])
neigh = []
for r in hg.itertuples():
    for off in range(-FLANK, FLANK + 1):
        j = r.idx + off
        if j < 0 or j >= r.n_on_contig:
            continue
        g = gpi.loc[(r.Short, r.contig, j)]
        neigh.append(dict(anchor=r.member, anchor_family=r.family, Short=r.Short, species=r.species,
                          contig=r.contig, offset=off, member=g.member, family=g.family,
                          bin=g.bin, start=g.start, end=g.end,
                          hrmA_models=prot_models.get(g.member, ""),
                          PF11001=g.member in prot_pf11001))
neigh = pd.DataFrame(neigh)
neigh.to_csv(A / "hrmA_neighbourhoods.tsv.gz", sep="\t", index=False, compression="gzip")
print("done", file=sys.stderr)
