#!/usr/bin/env python3
"""Step 4: genome-wide APH-adjacency baseline, HAC6 proximity, per-strain HrmA counts."""
import io, subprocess
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/hrmA_2026-09-24"
out = []
dom = [l.split() for l in __import__("gzip").open(A / "aph_hac6_allprot.domtblout.gz", "rt") if not l.startswith("#")]
aph = {f[0] for f in dom if f[3] == "APH"}
hac6 = {f[0] for f in dom if f[3].startswith("HAC6")}
gp = pd.read_csv(io.BytesIO(subprocess.run(["zstd", "-dc", str(P / "gene_positions.tsv.zst")], capture_output=True, check=True).stdout), sep="\t")
gp["member"] = gp.Short + "|" + gp.protein_id
gp = gp.sort_values(["Short", "contig", "start"]).reset_index(drop=True)
g = gp.groupby(["Short", "contig"])
gp["idx"] = g.cumcount()
hits = pd.read_csv(A / "hrmA_gene_locations.tsv", sep="\t")
lab = dict(zip(hits.member, hits.label))
gp["aph"] = gp.member.isin(aph); gp["hac6"] = gp.member.isin(hac6); gp["hrmA"] = gp.member.isin(set(hits.member))
gp["adj_aph"] = g.aph.shift(1, fill_value=False) | g.aph.shift(-1, fill_value=False)
gp["adj_hrmA"] = g.hrmA.shift(1, fill_value=False) | g.hrmA.shift(-1, fill_value=False)
a = gp[gp.hrmA]; b = gp[~gp.hrmA & ~gp.aph]
orr, p = fisher_exact([[int(a.adj_aph.sum()), int((~a.adj_aph).sum())], [int(b.adj_aph.sum()), int((~b.adj_aph).sum())]])
out += [("APH_genes_domain_GA", int(gp.aph.sum())), ("APH_genes_per_strain_median", gp.groupby("Short").aph.sum().median()),
        ("HrmA_genes_with_adjacent_APH", f"{int(a.adj_aph.sum())}/{len(a)}"),
        ("other_nonAPH_genes_with_adjacent_APH", f"{int(b.adj_aph.sum())}/{len(b)}"),
        ("fisher_OR", round(orr, 1)), ("fisher_p", p),
        ("APH_genes_adjacent_to_HrmA", f"{int(gp[gp.aph].adj_hrmA.sum())}/{int(gp.aph.sum())}"),
        ("HAC6_genes", int(gp.hac6.sum()))]
rows = []
for (s, c), d in gp[gp.hac6 | gp.hrmA].groupby(["Short", "contig"]):
    h = d[d.hrmA]
    for r in d[d.hac6].itertuples():
        if len(h):
            j = (h.idx - r.idx).abs().idxmin()
            rows.append((s, r.member, int(abs(h.idx[j] - r.idx)), lab[h.member[j]], abs(int(h.start[j]) - int(r.start))))
hd = pd.DataFrame(rows, columns=["Short", "hac6_gene", "genes_to_nearest_hrmA", "hrmA_label", "bp_to_nearest_hrmA"])
hd.to_csv(A / "hac6_to_hrmA_distance.tsv", sep="\t", index=False)
out += [("HAC6_genes_on_contig_with_HrmA", len(hd)), ("HAC6_within_10_genes", int((hd.genes_to_nearest_hrmA <= 10).sum())),
        ("HAC6_within_20_genes", int((hd.genes_to_nearest_hrmA <= 20).sum()))]
ss = pd.read_csv(P / "samplesheet.with_clades.csv"); inv = pd.read_csv(P / "strain_inventory.tsv", sep="\t")
d = ss[["Short", "Species"]].merge(inv, on="Short")
d["n_HrmA_PF28515"] = d.Short.map(hits[hits.hrmA_models == "HrmA"].groupby("Short").size()).fillna(0).astype(int)
d["n_HrmA_N_or_middle"] = d.Short.map(hits[hits.hrmA_models != "HrmA"].groupby("Short").size()).fillna(0).astype(int)
d.to_csv(A / "hrmA_genes_per_strain.tsv", sep="\t", index=False)
for sp_, x in d.groupby("Species"):
    out += [(f"{sp_}_PF28515_per_strain_median", x.n_HrmA_PF28515.median()), (f"{sp_}_PF28515_range", f"{x.n_HrmA_PF28515.min()}-{x.n_HrmA_PF28515.max()}"),
            (f"{sp_}_spearman_vs_n_contigs", spearmanr(x.n_HrmA_PF28515, x.n_contigs))]
out.append(("mannwhitney_Ci_vs_Cp", mannwhitneyu(d[d.Species.str.contains("immitis")].n_HrmA_PF28515, d[d.Species.str.contains("posadasii")].n_HrmA_PF28515)))
pd.DataFrame(out, columns=["metric", "value"]).to_csv(A / "hrmA_aph_hac6_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(out).to_string())
