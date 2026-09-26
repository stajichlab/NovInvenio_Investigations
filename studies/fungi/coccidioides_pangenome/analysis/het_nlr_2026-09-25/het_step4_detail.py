#!/usr/bin/env python3
"""Step 4: neighbour Pfam, top neighbours, allele-rep pairwise identity per locus,
relaxed NOD check, APH adjacency per family, exemplar locus views, candidate
ranking. Reads step 1-3 outputs and run_small_scans.sh outputs."""
import gzip, io, subprocess, sys, itertools
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/het_nlr_2026-09-25"
H = STUDY / "analysis/hrmA_2026-09-24"


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


census = pd.read_csv(A / "het_family_census.tsv", sep="\t")
lab = dict(zip(census.family, census.label))
cz = census.set_index("label")
freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")
inv = pd.read_csv(P / "strain_inventory.tsv", sep="\t").set_index("Short")
lt = pd.read_csv(A / "het_loci.tsv", sep="\t")
al = pd.read_csv(A / "het_locus_alleles.tsv", sep="\t")
occ = pd.read_csv(A / "het_locus_occupancy.tsv.gz", sep="\t")
hg = pd.read_csv(A / "het_gene_locations.tsv.gz", sep="\t")
nb = pd.read_csv(A / "het_neighbourhoods.tsv.gz", sep="\t")

# ---- Pfam per rep (hmmscan: col0 = model, col3 = query)
pf = defaultdict(set)
for path, op in ((H / "neighbour_reps_pfam.domtblout", open), (A / "reps_pfam.domtblout.gz", lambda p: gzip.open(p, "rt"))):
    for line in op(path):
        if not line.startswith("#"):
            f = line.split()
            pf[f[3]].add(f[0])
fam_pfam = {k: ",".join(sorted(v)) for k, v in pf.items()}
census["rep_pfam_full"] = census.family.map(fam_pfam).fillna("")

# ---- top neighbours per HET family
nb["nb_label"] = nb.family.map(lab).fillna("")
nb["pfam"] = nb.family.map(fam_pfam).fillna("")
top = []
for fam, d in nb[nb.offset != 0].groupby("anchor_family"):
    na = nb[(nb.anchor_family == fam) & (nb.offset == 0)].anchor.nunique()
    cnt = d.groupby("family").anchor.nunique().sort_values(ascending=False).head(8)
    for f, k in cnt.items():
        top.append(dict(label=lab[fam], anchor_family=fam, neighbour_family=f, neighbour_label=lab.get(f, ""),
                        frac_anchors=round(k / na, 3), n_anchors=int(k), median_offset=d[d.family == f].offset.median(),
                        bin=freq.bin.get(f, "NA"), pfam=fam_pfam.get(f, "")))
pd.DataFrame(top).to_csv(A / "het_top_neighbours.tsv", sep="\t", index=False)

# ---- rep pairwise identity, per locus
pw = pd.read_csv(A / "het_rep_pairwise.tsv", sep="\t", header=None,
                 names="q s pident length qlen slen qstart qend sstart send evalue bits".split())
pw = pw[pw.q != pw.s].sort_values("bits", ascending=False).drop_duplicates(["q", "s"])
pw["qcov"] = ((pw.qend - pw.qstart + 1) / pw.qlen).round(3)
pw["scov"] = ((pw.send - pw.sstart + 1) / pw.slen).round(3)
pw["q_label"], pw["s_label"] = pw.q.map(lab), pw.s.map(lab)
floc = dict(zip(census.family, census.locus))
pw["q_locus"], pw["s_locus"] = pw.q.map(floc), pw.s.map(floc)
pw.to_csv(A / "het_rep_pairwise_annotated.tsv", sep="\t", index=False)
best = {(r.q_label, r.s_label): r for r in pw.itertuples()}
rows = []
for L, d in al.groupby("locus"):
    fl = [f for f in d.family_label if (d.set_index("family_label").loc[f, ["Ci_carriers", "Cp_carriers"]].sum() >= 5)]
    for a, b in itertools.combinations(sorted(fl), 2):
        r = best.get((a, b)) or best.get((b, a))
        rows.append(dict(locus=L, fam_a=a, fam_b=b, arch_a=cz.modal_architecture[a], arch_b=cz.modal_architecture[b],
                         pident=r.pident if r else None, aln_len=r.length if r else None,
                         cov_a=(r.qcov if r.q_label == a else r.scov) if r else None,
                         cov_b=(r.scov if r.q_label == a else r.qcov) if r else None))
lp = pd.DataFrame(rows)
lp.to_csv(A / "het_locus_allele_pairwise.tsv", sep="\t", index=False)

# ---- relaxed NOD check
rel = defaultdict(set)
for line in gzip.open(A / "het_prots_nod_relaxed.domtblout.gz", "rt"):
    if not line.startswith("#"):
        f = line.split()
        rel[f[0]].add(f[3])
hg["relaxed_NOD_models"] = hg.member.map(lambda m: ",".join(sorted(rel.get(m, ())))).fillna("")
hg["relaxed_NOD_any"] = hg.relaxed_NOD_models != ""
rn = hg.groupby("label").agg(n=("member", "size"), relaxed_NOD=("relaxed_NOD_any", "sum"),
                             relaxed_NACHT_or_NBARC=("relaxed_NOD_models", lambda s: int(s.str.contains("NACHT|NB-ARC").sum())),
                             median_len=("tlen", "median")).reset_index()
rc = hg.groupby("het_class").agg(n=("member", "size"), relaxed_NOD=("relaxed_NOD_any", "sum"),
                                 relaxed_NACHT_or_NBARC=("relaxed_NOD_models", lambda s: int(s.str.contains("NACHT|NB-ARC").sum())),
                                 median_len=("tlen", "median")).reset_index()
rn.to_csv(A / "het_relaxed_nod_by_family.tsv", sep="\t", index=False)
rc.to_csv(A / "het_relaxed_nod_by_class.tsv", sep="\t", index=False)

# ---- APH adjacency per family (from step-1 neighbourhoods + hrmA APH domtblout)
dom = [l.split() for l in gzip.open(H / "aph_hac6_allprot.domtblout.gz", "rt") if not l.startswith("#")]
aph = {f[0] for f in dom if f[3] == "APH"}
adj = nb[nb.offset.abs() == 1].assign(aph=lambda d: d.member.isin(aph)).groupby("anchor").aph.any()
hg["APH_adjacent"] = hg.member.map(adj).fillna(False)
census["frac_APH_adjacent"] = census.label.map(hg.groupby("label").APH_adjacent.mean().round(3))
hg.to_csv(A / "het_gene_locations.tsv.gz", sep="\t", index=False, compression="gzip")

# ---- per-locus exemplar views (best-assembled occupied strain per species)
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
hcls = dict(zip(hg.member, hg.architecture))
ex, exs = [], []
occ["n_contigs"] = occ.Short.map(inv.n_contigs)
mkr = dict(zip(lt.locus, lt.markers.fillna("")))
for L in lt[lt.n_markers > 0].locus:
    mset = set(mkr[L].split(","))
    for spc in ("Ci", "Cp"):
        d = occ[(occ.locus == L) & (occ.species == spc) & (occ.state == "occupied")].sort_values("n_contigs")
        # one exemplar per distinct top allele combination (up to 3), best assembled
        for comb, dd in sorted(d.groupby("families"), key=lambda x: -len(x[1]))[:3]:
            r = dd.iloc[0]
            genes = r.genes.split(",")
            g = gp[gp.member.isin(genes)]
            for contig, gg in g.groupby("contig"):
                lo, hi = gg.idx.min() - 3, gg.idx.max() + 3
                win = gp[(gp.Short == r.Short) & (gp.contig == contig) & (gp.idx >= lo) & (gp.idx <= hi)]
                exs.append(dict(locus=L, species=spc, combination=comb, n_strains_with_combination=len(dd),
                                exemplar=f"{r.Short}:{contig}:{int(gg.start.min())}-{int(gg.end.max())}",
                                n_contigs=int(r.n_contigs), het_genes=",".join(gg.member)))
                for w in win.itertuples():
                    ex.append(dict(locus=L, species=spc, combination=comb, Short=r.Short, contig=contig, idx=w.idx,
                                   start=w.start, end=w.end, member=w.member, family_label=lab.get(w.family, ""),
                                   is_marker=w.family in mset, bin=freq.bin.get(w.family, "NA"),
                                   het_architecture=hcls.get(w.member, ""), rep_pfam=fam_pfam.get(w.family, "")))
pd.DataFrame(ex).to_csv(A / "het_exemplar_loci.tsv", sep="\t", index=False)
pd.DataFrame(exs).to_csv(A / "het_exemplar_summary.tsv", sep="\t", index=False)
census.to_csv(A / "het_family_census.tsv", sep="\t", index=False)
print("done", file=sys.stderr)
