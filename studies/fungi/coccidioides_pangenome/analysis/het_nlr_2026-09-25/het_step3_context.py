#!/usr/bin/env python3
"""Step 3: location vs contig ends, per-strain copy number tests, HrmA/APH
proximity, islands, pairs, rescue-on-paralog overlap. Reads step 1/2 outputs.
Read-only on pipeline outputs."""
import gzip, io, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/het_nlr_2026-09-25"
H = STUDY / "analysis/hrmA_2026-09-24"
NLR = ["NLR_NOD", "NLR_like_AAA"]


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


census = pd.read_csv(A / "het_family_census.tsv", sep="\t")
lab = dict(zip(census.family, census.label))
locus = dict(zip(census.family, census.locus))
fams = set(census.family)
freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
sp = dict(zip(ss.Short, ss.Species.map({"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"})))
inv = pd.read_csv(P / "strain_inventory.tsv", sep="\t").set_index("Short")
hg = pd.read_csv(A / "het_gene_locations.tsv.gz", sep="\t")
summary = []

# ---------------- all genes
gp = zread(P / "gene_positions.tsv.zst")
gp["member"] = gp.Short + "|" + gp.protein_id
clen = pd.read_csv(H / "contig_lengths.tsv.gz", sep="\t", header=None, names=["Short", "contig", "contig_len"])
gp = gp.merge(clen, on=["Short", "contig"], how="left")
gp["dist_end"] = pd.concat([gp.start - 1, gp.contig_len - gp.end], axis=1).min(axis=1)
gp = gp.sort_values(["Short", "contig", "start"]).reset_index(drop=True)
g = gp.groupby(["Short", "contig"])
gp["idx"] = g.cumcount()
gp["het"] = gp.member.isin(set(hg.member))
hrm = pd.read_csv(H / "hrmA_gene_locations.tsv", sep="\t")
hcen = pd.read_csv(H / "hrmA_family_census.tsv", sep="\t")
hlocus = dict(zip(hcen.label, hcen.locus))
gp["hrmA"] = gp.member.isin(set(hrm.member))
hrm_lab = dict(zip(hrm.member, hrm.label))
dom = [l.split() for l in gzip.open(H / "aph_hac6_allprot.domtblout.gz", "rt") if not l.startswith("#")]
aph = {f[0] for f in dom if f[3] == "APH"}
gp["aph"] = gp.member.isin(aph)
g = gp.groupby(["Short", "contig"])
gp["adj_aph"] = g.aph.shift(1, fill_value=False) | g.aph.shift(-1, fill_value=False)

# --------------- 1. contig-end location tests
rest = gp[~gp.het]
tests = []
sets = [("all_HET_NLR", hg)] + [(c, hg[hg.het_class == c]) for c in sorted(hg.het_class.unique())] + \
       [("NLR_NOD+NLR_like_AAA", hg[hg.het_class.isin(NLR)])]
for name, sub in sets:
    for tname, thr in (("within_20kb", 20000), ("within_50kb", 50000)):
        a = int((sub.dist_end <= thr).sum()); b = len(sub) - a
        c = int((rest.dist_end <= thr).sum()); d = len(rest) - c
        orr, p = fisher_exact([[a, b], [c, d]])
        tests.append(dict(set=name, test=tname, n=len(sub), yes=a, frac=round(a / len(sub), 4),
                          rest_frac=round(c / len(rest), 4), odds_ratio=round(orr, 3), fisher_p=p,
                          median_contig_len=int(sub.contig_len.median()), median_dist_end=int(sub.dist_end.median())))
    tests.append(dict(set=name, test="nearest_end_telomeric", n=len(sub), yes=int((sub.nearest_end_telo == 1).sum()),
                      frac=round((sub.nearest_end_telo == 1).mean(), 4)))
pd.DataFrame(tests).to_csv(A / "het_location_tests.tsv", sep="\t", index=False)

# per family location
fl = []
for fam, d in hg.groupby("family"):
    fl.append(dict(label=lab[fam], family=fam, n_genes=len(d), median_dist_end=int(d.dist_end.median()),
                   frac_le20kb=round((d.dist_end <= 20000).mean(), 3), frac_le50kb=round((d.dist_end <= 50000).mean(), 3),
                   median_contig_len=int(d.contig_len.median()), median_genes_to_end=d.genes_to_end.median(),
                   frac_genes_to_end_0=round((d.genes_to_end == 0).mean(), 3)))
pd.DataFrame(fl).sort_values("label").to_csv(A / "het_family_locations.tsv", sep="\t", index=False)

# --------------- 2. HrmA / APH proximity
near = []
hr_rows = gp[gp.hrmA]
hr_by = {k: v for k, v in hr_rows.groupby(["Short", "contig"])}
gp_het = gp[gp.het]
for r in gp_het.itertuples():
    x = hr_by.get((r.Short, r.contig))
    if x is None:
        near.append((r.member, None, "", "", None)); continue
    dd = (x.idx - r.idx).abs()
    j = dd.idxmin()
    lb = hrm_lab[x.member[j]]
    near.append((r.member, int(dd[j]), lb, hlocus.get(lb, ""), abs(int(x.start[j]) - int(r.start))))
near = pd.DataFrame(near, columns=["member", "genes_to_nearest_hrmA", "hrmA_label", "hrmA_locus", "bp_to_nearest_hrmA"])
near = near.merge(hg[["member", "label", "het_class", "species"]], on="member")
near.to_csv(A / "het_to_hrmA_distance.tsv.gz", sep="\t", index=False, compression="gzip")
# background: fraction of all non-HET genes within 10 genes of an HrmA gene
hr_idx = defaultdict(list)
for (s, c), x in hr_by.items():
    hr_idx[(s, c)] = x.idx.values
bg_within = 0
nonhet = gp[~gp.het & ~gp.hrmA]
for (s, c), x in nonhet.groupby(["Short", "contig"]):
    if (s, c) in hr_idx:
        hv = hr_idx[(s, c)]
        bg_within += int((np.abs(x.idx.values[:, None] - hv[None, :]).min(axis=1) <= 10).sum())
a = int((near.genes_to_nearest_hrmA <= 10).sum())
orr, p = fisher_exact([[a, len(near) - a], [bg_within, len(nonhet) - bg_within]])
summary += [("HET_NLR_genes", len(near)), ("HET_NLR_within10genes_of_HrmA", a),
            ("nonHET_nonHrmA_within10genes_of_HrmA", f"{bg_within}/{len(nonhet)}"),
            ("fisher_OR_within10_HrmA", round(orr, 3)), ("fisher_p_within10_HrmA", p)]
a = gp_het.adj_aph.sum(); b = nonhet.adj_aph.sum()
orr, p = fisher_exact([[int(a), len(gp_het) - int(a)], [int(b), len(nonhet) - int(b)]])
summary += [("HET_NLR_with_adjacent_APH", f"{int(a)}/{len(gp_het)}"), ("nonHET_with_adjacent_APH", f"{int(b)}/{len(nonhet)}"),
            ("fisher_OR_adjacent_APH", round(orr, 3)), ("fisher_p_adjacent_APH", p)]
del nonhet, rest

# --------------- 3. per-strain copy number
d = ss[["Short", "Species"]].merge(inv.reset_index(), on="Short")
d["species"] = d.Short.map(sp)
for c in ["all_HET_NLR"] + sorted(hg.het_class.unique()) + ["NLR_total"]:
    if c == "all_HET_NLR":
        s = hg.groupby("Short").size()
    elif c == "NLR_total":
        s = hg[hg.het_class.isin(NLR)].groupby("Short").size()
    else:
        s = hg[hg.het_class == c].groupby("Short").size()
    d[f"n_{c}"] = d.Short.map(s).fillna(0).astype(int)
d.to_csv(A / "het_genes_per_strain.tsv", sep="\t", index=False)
cn = []
for c in [x for x in d.columns if x.startswith("n_") and x != "n_contigs" and x != "n50"]:
    for subset, dd in (("all_strains", d), ("representatives_only", d[d.is_representative == 1])):
        ci, cp = dd[dd.species == "Ci"][c], dd[dd.species == "Cp"][c]
        u, p = mannwhitneyu(ci, cp)
        row = dict(count=c, subset=subset, Ci_n=len(ci), Cp_n=len(cp), Ci_median=ci.median(), Ci_mean=round(ci.mean(), 2),
                   Ci_range=f"{ci.min()}-{ci.max()}", Cp_median=cp.median(), Cp_mean=round(cp.mean(), 2),
                   Cp_range=f"{cp.min()}-{cp.max()}", mannwhitney_U=u, mannwhitney_p=p)
        for lb, x in (("Ci", dd[dd.species == "Ci"]), ("Cp", dd[dd.species == "Cp"])):
            rho, pr = spearmanr(x[c], x.n_contigs)
            row[f"{lb}_spearman_rho_vs_n_contigs"] = round(rho, 3)
            row[f"{lb}_spearman_p"] = pr
            rho, pr = spearmanr(x[c], x.total_length)
            row[f"{lb}_spearman_rho_vs_assembly_len"] = round(rho, 3)
            row[f"{lb}_spearman_p_len"] = pr
        cn.append(row)
pd.DataFrame(cn).to_csv(A / "het_copy_number_tests.tsv", sep="\t", index=False)

# --------------- 4. islands
isl = pd.read_csv(P / "report_tables/islands_with_domains.tsv", sep="\t")
irows = []
for r in isl.itertuples():
    for fam in set(r.member_families.split(",")) & fams:
        irows.append(dict(label=lab[fam], family=fam, locus=locus[fam], island_n_strains=r.n_strains,
                          example_strain=r.example_strain, island_size=r.island_size,
                          n_supporting_pairs=r.n_supporting_pairs, classifications=r.classifications,
                          locus_id=r.locus_id, locus_contig=r.locus_contig, locus_start=r.locus_start,
                          locus_end=r.locus_end, pfam_domains=r.pfam_domains))
irows = pd.DataFrame(irows)
irows.to_csv(A / "het_island_membership.tsv", sep="\t", index=False)
summary += [("islands_total", len(isl)), ("islands_with_HET_NLR_family", irows.locus_id.nunique() if len(irows) else 0),
            ("HET_NLR_families_in_islands", irows.family.nunique() if len(irows) else 0)]

# --------------- 5. pairs
pc = zread(P / "pair_classification.tsv.zst")
pp = pc[pc.family_a.isin(fams) | pc.family_b.isin(fams)].copy()
pp["het_family"] = pp.family_a.where(pp.family_a.isin(fams), pp.family_b)
pp["partner"] = pp.family_b.where(pp.family_a.isin(fams), pp.family_a)
pp["label"] = pp.het_family.map(lab)
pp["partner_label"] = pp.partner.map(lab).fillna("")
pp["partner_bin"] = pp.partner.map(freq.bin)
pp.to_csv(A / "het_pairs.tsv.gz", sep="\t", index=False, compression="gzip")
tested = set(pc.family_a) | set(pc.family_b)
census["in_pair_table"] = census.family.isin(tested)
summary += [("HET_NLR_families_in_pair_table", int(census.in_pair_table.sum())), ("pairs_with_HET_NLR_family", len(pp))]
del pc

# --------------- 6. rescue overlap (genome_only cells)
rp = pd.read_csv(P / "rescue_positions.tsv", sep="\t")
rp = rp[rp.family.isin(fams)]
m2f = {}
need = set(rp.Short)
with open(P / "cluster/tier1_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        if mem.split("|")[0] in need:
            m2f[mem] = rep
gpn = gp[gp.Short.isin(need)].copy()
gpn["family"] = gpn.member.map(m2f)
gidx = {k: v for k, v in gpn.groupby(["Short", "contig"])}
hset = set(hg.member)
orow = []
for r in rp.itertuples():
    x = gidx.get((r.Short, r.contig))
    ov = x[(x.start - 500 <= r.start) & (x.end + 500 >= r.start)] if x is not None else pd.DataFrame()
    orow.append(dict(label=lab[r.family], family=r.family, Short=r.Short, species=sp[r.Short], contig=r.contig,
                     start=r.start,
                     overlapping_gene=";".join(ov.member) if len(ov) else "",
                     overlapping_family=";".join(ov.family.map(lambda y: lab.get(y, y))) if len(ov) else "",
                     overlapping_is_het_nlr=bool(ov.member.isin(hset).any()) if len(ov) else False))
ro = pd.DataFrame(orow)
ro.to_csv(A / "het_rescue_overlap.tsv", sep="\t", index=False)
if len(ro):
    summary += [("rescue_positions_HET_families", len(ro)), ("rescue_overlapping_any_gene", int((ro.overlapping_gene != "").sum())),
                ("rescue_overlapping_HET_NLR_gene", int(ro.overlapping_is_het_nlr.sum()))]
census.to_csv(A / "het_family_census.tsv", sep="\t", index=False)
pd.DataFrame(summary, columns=["metric", "value"]).to_csv(A / "het_context_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(summary).to_string(), file=sys.stderr)
