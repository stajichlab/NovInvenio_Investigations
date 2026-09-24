#!/usr/bin/env python3
"""Step 2: neighbourhood conservation, telomere ends, islands, pairs, rescue
overlap, A. fumigatus HAC homologs. Reads step-1 outputs. Read-only on pipeline."""
import io, subprocess, sys, itertools, random
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/hrmA_2026-09-24"
random.seed(1)


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


def parse_domtbl(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        rows.append((f[3], f[0], float(f[12])))  # query(protein) , model name, i-Evalue (hmmscan: target=model, query=protein)
    return rows


census = pd.read_csv(A / "hrmA_family_census.tsv", sep="\t")
census = census.drop(columns=[c for c in ("label",) if c in census.columns])
hf = list(census.family)
fam_label = {f: f"F{i+1:02d}" for i, f in enumerate(hf)}
census.insert(0, "label", census.family.map(fam_label))
census.to_csv(A / "hrmA_family_census.tsv", sep="\t", index=False)
freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")
loc = pd.read_csv(A / "hrmA_gene_locations.tsv", sep="\t")
nb = pd.read_csv(A / "hrmA_neighbourhoods.tsv.gz", sep="\t")
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
species = dict(zip(ss.Short, ss.Species.map({"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"})))
inv = pd.read_csv(P / "strain_inventory.tsv", sep="\t").set_index("Short")

# --- Pfam on neighbour reps (hmmscan: col0 = model name, col3 = query protein)
pf = defaultdict(set)
for line in open(A / "neighbour_reps_pfam.domtblout"):
    if line.startswith("#"):
        continue
    f = line.split()
    pf[f[3]].add(f[0])
fam_pfam = {k: ",".join(sorted(v)) for k, v in pf.items()}

# --- 1. per-family location + telomere
tel = pd.read_csv(A / "contig_telomere_ends.tsv.gz", sep="\t")
loc = loc.merge(tel[["Short", "contig", "left_telo", "right_telo"]], on=["Short", "contig"], how="left")
loc["near_left"] = (loc.start - 1) <= (loc.contig_len - loc.end)
loc["nearest_end_telo"] = loc.left_telo.where(loc.near_left, loc.right_telo)
loc["genes_to_end"] = pd.concat([loc.idx, loc.n_on_contig - 1 - loc.idx], axis=1).min(axis=1)
loc["label"] = loc.family.map(fam_label)
loc.to_csv(A / "hrmA_gene_locations.tsv", sep="\t", index=False)

rows = []
for fam, d in loc.groupby("family"):
    rows.append(dict(label=fam_label[fam], family=fam, n_genes=len(d), n_strains=d.Short.nunique(),
                     median_dist_end_bp=int(d.dist_end.median()),
                     frac_le20kb=round((d.dist_end <= 20000).mean(), 3),
                     frac_le50kb=round((d.dist_end <= 50000).mean(), 3),
                     median_contig_len=int(d.contig_len.median()),
                     median_genes_to_contig_end=d.genes_to_end.median(),
                     frac_nearest_end_telomeric=round(d.nearest_end_telo.mean(), 3),
                     frac_on_contig_with_any_telomere=round(((d.left_telo + d.right_telo) > 0).mean(), 3)))
famloc = pd.DataFrame(rows).sort_values("label")
famloc.to_csv(A / "hrmA_family_locations.tsv", sep="\t", index=False)

# genome-wide background for telomere-at-nearest-end (all genes)
gp = zread(P / "gene_positions.tsv.zst")
cl = tel.rename(columns={"len": "contig_len"})
gp = gp.merge(cl, on=["Short", "contig"], how="left")
gp["dl"] = gp.start - 1
gp["dr"] = gp.contig_len - gp.end
gp["dist_end"] = gp[["dl", "dr"]].min(axis=1)
gp["nearest_end_telo"] = gp.left_telo.where(gp.dl <= gp.dr, gp.right_telo)
bgrows = []
for label, sub in [("all_genes", gp), ("HrmA_hit_genes", loc)]:
    bgrows.append(dict(set=label, n=len(sub), frac_le20kb=round((sub.dist_end <= 20000).mean(), 4),
                       frac_le50kb=round((sub.dist_end <= 50000).mean(), 4),
                       frac_nearest_end_telomeric=round(sub.nearest_end_telo.mean(), 4),
                       frac_le50kb_and_telomeric=round(((sub.dist_end <= 50000) & (sub.nearest_end_telo == 1)).mean(), 4)))
# Fisher tests HrmA vs rest of genome
hit = set(loc.member)
gp["member"] = gp.Short + "|" + gp.protein_id
rest = gp[~gp.member.isin(hit)]
tests = []
for name, fn in [("within_20kb", lambda d: d.dist_end <= 20000), ("within_50kb", lambda d: d.dist_end <= 50000),
                 ("nearest_end_telomeric", lambda d: d.nearest_end_telo == 1),
                 ("within_50kb_and_telomeric", lambda d: (d.dist_end <= 50000) & (d.nearest_end_telo == 1))]:
    a = int(fn(loc).sum()); b = len(loc) - a; c = int(fn(rest).sum()); dd = len(rest) - c
    orr, p = fisher_exact([[a, b], [c, dd]])
    tests.append(dict(test=name, hrmA_yes=a, hrmA_no=b, rest_yes=c, rest_no=dd, odds_ratio=round(orr, 3), fisher_p=p))
pd.DataFrame(bgrows).to_csv(A / "hrmA_location_telomere_background.tsv", sep="\t", index=False)
pd.DataFrame(tests).to_csv(A / "hrmA_location_tests.tsv", sep="\t", index=False)
tel_summary = dict(n_contigs=len(tel), frac_contigs_any_telo=round(((tel.left_telo + tel.right_telo) > 0).mean(), 4),
                   n_telomeric_ends=int(tel.left_telo.sum() + tel.right_telo.sum()))
print("telomere summary", tel_summary, file=sys.stderr)
pd.DataFrame([tel_summary]).to_csv(A / "telomere_end_summary.tsv", sep="\t", index=False)
per_strain_tel = tel.groupby("Short").apply(lambda d: int(d.left_telo.sum() + d.right_telo.sum()))
per_strain_tel.rename("n_telomeric_ends").to_csv(A / "telomeric_ends_per_strain.tsv", sep="\t")
del gp, rest

# --- 2. neighbourhood conservation
nb["label"] = nb.anchor_family.map(fam_label)
nb["nb_bin"] = nb.family.map(freq.bin)
nb["pfam"] = nb.family.map(fam_pfam).fillna("")
nb.to_csv(A / "hrmA_neighbourhoods.tsv.gz", sep="\t", index=False, compression="gzip")
cons, topn = [], []
for fam, d in nb.groupby("anchor_family"):
    anchors = d.anchor.unique()
    sets = {a: set(g[g.offset != 0].family) for a, g in d.groupby("anchor")}
    full = sum(1 for a, g in d.groupby("anchor") if g.offset.min() == -10 and g.offset.max() == 10)
    # ordered +-2 signature (orientation-canonical)
    sigs = []
    for a, g in d.groupby("anchor"):
        o = dict(zip(g.offset, g.family))
        if all(k in o for k in (-2, -1, 1, 2)):
            t = (o[-2], o[-1], o[1], o[2])
            sigs.append(min(t, t[::-1]))
    sc = Counter(sigs)
    modal_frac = round(sc.most_common(1)[0][1] / len(sigs), 3) if sigs else None
    al = list(sets)
    pairs = list(itertools.combinations(al, 2))
    if len(pairs) > 5000:
        pairs = random.sample(pairs, 5000)
    jac = [len(sets[x] & sets[y]) / len(sets[x] | sets[y]) for x, y in pairs if sets[x] | sets[y]]
    cons.append(dict(label=fam_label[fam], family=fam, n_anchors=len(anchors),
                     n_full_window=full, n_with_pm2=len(sigs), modal_pm2_order_frac=modal_frac,
                     n_distinct_pm2_orders=len(sc),
                     mean_pairwise_jaccard_pm10=round(sum(jac) / len(jac), 3) if jac else None))
    cnt = Counter(f for a in anchors for f in sets[a])
    for f, k in cnt.most_common(8):
        offs = d[d.family == f].offset
        topn.append(dict(label=fam_label[fam], anchor_family=fam, neighbour_family=f,
                         frac_anchors=round(k / len(anchors), 3), n_anchors=k,
                         median_offset=offs.median(), bin=freq.bin.get(f, "NA"),
                         strain_count_Ci=freq.strain_count.get(f, None), pfam=fam_pfam.get(f, "")))
pd.DataFrame(cons).sort_values("label").to_csv(A / "hrmA_neighbourhood_conservation.tsv", sep="\t", index=False)
pd.DataFrame(topn).to_csv(A / "hrmA_top_neighbours.tsv", sep="\t", index=False)

# --- 3. exemplar neighbourhoods (2 best-assembled carriers per species per family, families >=10 anchors)
def gff_products(short, members):
    want = {m.split("|", 1)[1] for m in members}
    out = {}
    for line in open(STUDY / f"data_dir/gff3/{short}.gff3"):
        if "\tmRNA\t" not in line:
            continue
        attr = line.rstrip("\n").split("\t")[8]
        kv = dict(x.split("=", 1) for x in attr.strip(";").split(";") if "=" in x)
        if kv.get("ID") in want:
            out[f"{short}|{kv['ID']}"] = kv.get("product", "")
    return out

ex_rows = []
big = census[census.n_members >= 10].family
for fam in big:
    d = nb[nb.anchor_family == fam]
    anc = d[d.offset == 0][["anchor", "Short", "species"]].drop_duplicates()
    anc["n_contigs"] = anc.Short.map(inv.n_contigs)
    for sp in ("Ci", "Cp"):
        for a in anc[anc.species == sp].sort_values("n_contigs").anchor.head(2):
            g = d[d.anchor == a].sort_values("offset")
            prod = gff_products(g.Short.iloc[0], g.member)
            for r in g.itertuples():
                ex_rows.append(dict(label=fam_label[fam], anchor=a, species=sp, n_contigs=inv.n_contigs[r.Short],
                                    offset=r.offset, member=r.member, family=r.family,
                                    family_label=fam_label.get(r.family, ""), bin=r.nb_bin, pfam=r.pfam,
                                    hrmA_models=r.hrmA_models if isinstance(r.hrmA_models, str) else "",
                                    product=prod.get(r.member, "")))
pd.DataFrame(ex_rows).to_csv(A / "hrmA_exemplar_neighbourhoods.tsv", sep="\t", index=False)

# --- 4. islands
isl = pd.read_csv(P / "report_tables/islands_with_domains.tsv", sep="\t")
irows = []
for r in isl.itertuples():
    mems = set(r.member_families.split(","))
    for fam in mems & set(hf):
        irows.append(dict(label=fam_label[fam], family=fam, island_n_strains=r.n_strains, example_strain=r.example_strain,
                          island_size=r.island_size, n_supporting_pairs=r.n_supporting_pairs,
                          classifications=r.classifications, locus_id=r.locus_id, locus_contig=r.locus_contig,
                          locus_start=r.locus_start, locus_end=r.locus_end,
                          pfam_domains=r.pfam_domains))
pd.DataFrame(irows).to_csv(A / "hrmA_island_membership.tsv", sep="\t", index=False)
print("island rows:", len(irows), file=sys.stderr)

# --- 5. pairs
pc = zread(P / "pair_classification.tsv.zst")
s = set(hf)
pp = pc[pc.family_a.isin(s) | pc.family_b.isin(s)].copy()
pp["hrmA_family"] = pp.family_a.where(pp.family_a.isin(s), pp.family_b)
pp["partner"] = pp.family_b.where(pp.family_a.isin(s), pp.family_a)
pp["label"] = pp.hrmA_family.map(fam_label)
pp["partner_label"] = pp.partner.map(fam_label).fillna("")
pp["partner_bin"] = pp.partner.map(freq.bin)
pp["partner_pfam"] = pp.partner.map(fam_pfam).fillna("")
nbfr = nb[nb.offset != 0].groupby(["anchor_family", "family"]).anchor.nunique()
nanch = nb[nb.offset == 0].groupby("anchor_family").anchor.nunique()
pp["partner_frac_anchors_within10"] = [round(nbfr.get((h, q), 0) / nanch[h], 3) for h, q in zip(pp.hrmA_family, pp.partner)]
pp.to_csv(A / "hrmA_pairs.tsv", sep="\t", index=False)
print("pairs involving HrmA families:", len(pp), file=sys.stderr)
tested = set(pc.family_a) | set(pc.family_b)
pd.DataFrame(dict(label=[fam_label[f] for f in hf], family=hf, in_pair_table=[f in tested for f in hf])).to_csv(
    A / "hrmA_in_pair_table.tsv", sep="\t", index=False)

# --- 6. rescue overlap for genome_only cells
rp = pd.read_csv(P / "rescue_positions.tsv", sep="\t")
rp = rp[rp.family.isin(s)]
gp2 = zread(P / "gene_positions.tsv.zst")
m2f = {}
need = set(rp.Short)
gp2 = gp2[gp2.Short.isin(need)]
members_needed = set(gp2.Short + "|" + gp2.protein_id)
with open(P / "cluster/tier1_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        if mem in members_needed:
            m2f[mem] = rep
gp2["family"] = (gp2.Short + "|" + gp2.protein_id).map(m2f)
gidx = {k: v for k, v in gp2.groupby(["Short", "contig"])}
hitset = set(loc.member)
orow = []
for r in rp.itertuples():
    g = gidx.get((r.Short, r.contig))
    ov = g[(g.start - 500 <= r.start) & (g.end + 500 >= r.start)] if g is not None else pd.DataFrame()
    orow.append(dict(label=fam_label[r.family], family=r.family, Short=r.Short, contig=r.contig, start=r.start,
                     overlapping_gene=";".join(ov.Short + "|" + ov.protein_id) if len(ov) else "",
                     overlapping_family=";".join(ov.family.map(lambda x: fam_label.get(x, x))) if len(ov) else "",
                     overlapping_is_hrmA_hit=any((ov.Short + "|" + ov.protein_id).isin(hitset)) if len(ov) else False))
ro = pd.DataFrame(orow)
ro.to_csv(A / "hrmA_rescue_overlap.tsv", sep="\t", index=False)

# --- 7. A. fumigatus HAC-region homologs
dm = pd.read_csv(A / "af293_hac_region_vs_reps.tsv", sep="\t", header=None,
                 names="q s pident length qlen slen qstart qend sstart send evalue bits".split())
dm["qcov"] = ((dm.qend - dm.qstart + 1) / dm.qlen).round(2)
nbfam_anchor = nb[nb.offset != 0].groupby("family").apply(lambda d: d.groupby("label").anchor.nunique().to_dict())
dm["hrmA_family_label"] = dm.s.map(fam_label).fillna("")
dm["in_hrmA_neighbourhoods"] = dm.s.map(lambda f: nbfam_anchor.get(f, {}))
dm["bin"] = dm.s.map(freq.bin)
dm.to_csv(A / "af293_hac_region_hits_annotated.tsv", sep="\t", index=False)
print("done", file=sys.stderr)
