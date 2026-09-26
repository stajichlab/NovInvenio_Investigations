#!/usr/bin/env python3
"""Step 1: HET/NLR domain census over all 529 proteomes.

Protein-level domain architecture and class, mapping to tier-1 families,
family presence/absence per species (Fisher, BH), gene locations and +/-10
neighbourhoods. Pattern follows analysis/hrmA_2026-09-24/hrmA_step1_families.py.
Read-only on pipeline outputs; writes into this directory.
Run with the NII pixi python."""
import gzip, io, subprocess, sys
from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/het_nlr_2026-09-25"
H = STUDY / "analysis/hrmA_2026-09-24"
FLANK = 10

# ---- domain groups (see report, Methods)
NOD = {"NACHT", "NB-ARC"}
NLR_ACC = {"NACHT_N", "NACHT_sigma", "NPHP3_N", "NPHP3_hel", "TPR_NPHP3", "WHD_GPIID", "WHD_NACHT_C",
           "WHD_NWD1", "WHD_APAF1", "WHD_CED4", "Beta-prop_NWD2_C"}
WHD_AMBIG = {"WHD_AAA_fung", "WHD_Fungal_DR"}
AAA = {"AAA_16", "AAA_22"}
HET_N = {"HET", "HeLo", "HET-S", "HET-s_218-289", "SesA", "Goodbye", "Het-6_barrel"}
EFFECTOR = {"PNP_UDP_1", "Patatin", "CHAT", "Peptidase_C14", "TIR", "TIR_2"}
HETC = {"Het-C"}


def rep_group(m):
    if m.startswith("Ank"):
        return "Ank"
    if m.startswith("WD40"):
        return "WD40"
    if m.startswith("TPR_") and m != "TPR_NPHP3":
        return "TPR"
    return m


SHORT = {"NACHT": "NACHT", "NB-ARC": "NB-ARC", "AAA_16": "AAA", "AAA_22": "AAA", "HET": "HET", "HeLo": "HeLo",
         "HET-S": "HET-S", "HET-s_218-289": "HETs_PFD", "SesA": "SesA", "Goodbye": "Goodbye",
         "Het-6_barrel": "Het6_barrel", "PNP_UDP_1": "PNP_UDP", "Patatin": "Patatin", "CHAT": "CHAT",
         "Peptidase_C14": "Caspase", "TIR": "TIR", "TIR_2": "TIR", "Het-C": "Het-C",
         "NACHT_N": "NACHT_N", "NACHT_sigma": "sigma", "NPHP3_N": "NPHP3_N", "NPHP3_hel": "NPHP3_hel",
         "TPR_NPHP3": "TPR_NPHP3", "Beta-prop_NWD2_C": "NWD2_prop"}


def label(m):
    g = rep_group(m)
    if g in ("Ank", "WD40", "TPR"):
        return g
    if m.startswith("WHD_"):
        return "WHD" if m in NLR_ACC else "WHD_fung"
    return SHORT.get(m, m)


def classify(models):
    s = set(models)
    reps = {rep_group(m) for m in s} & {"Ank", "WD40", "TPR"}
    if s & NOD:
        return "NLR_NOD"
    if s & AAA and (s & NLR_ACC or s & HET_N):
        return "NLR_like_AAA"
    if s & NLR_ACC:
        return "NLR_accessory_only"
    if "HET" in s:
        return "HET_domain"
    if s & (HET_N - {"HET"}):
        return "HeLo_SesA_Goodbye_HETs"
    if s & HETC:
        return "Het-C"
    if s & AAA and reps:
        return "AAA_repeat_lowconf"
    return "other"


def zread(path, **kw):
    out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(out), sep="\t", **kw)


def main():
    rows = []
    with gzip.open(A / "allprot.domtblout.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            rows.append((f[0], int(f[2]), f[3], float(f[7]), float(f[12]), float(f[13]), int(f[19]), int(f[20])))
    dom = pd.DataFrame(rows, columns=["target", "tlen", "model", "full_score", "dom_iE", "dom_score", "env_from", "env_to"])
    print("domain rows:", len(dom), "proteins:", dom.target.nunique(), file=sys.stderr)

    # architecture: ordered labels by env_from, consecutive duplicates collapsed;
    # AAA dropped from the string where it overlaps a NOD hit on the same protein.
    arch, cls, mods = {}, {}, {}
    for t, d in dom.groupby("target"):
        d = d.sort_values("env_from")
        ms = list(d.model)
        nod = d[d.model.isin(NOD)]
        labs = []
        for r in d.itertuples():
            if r.model in AAA and len(nod) and ((nod.env_from <= r.env_to) & (nod.env_to >= r.env_from)).any():
                continue
            lb = label(r.model)
            if not labs or labs[-1] != lb:
                labs.append(lb)
        # merge NACHT/NB-ARC overlap into one label
        arch[t] = "-".join(labs)
        cls[t] = classify(ms)
        mods[t] = ",".join(sorted(set(ms)))
    prot = pd.DataFrame(dict(target=list(arch), architecture=list(arch.values()),
                             het_class=[cls[t] for t in arch], models=[mods[t] for t in arch]))
    prot["tlen"] = prot.target.map(dom.groupby("target").tlen.first())
    prot["Short"] = prot.target.str.split("|").str[0]
    prot["has_Ank"] = prot.architecture.str.contains(r"(?:^|-)Ank(?:-|$)")
    prot["has_WD40"] = prot.architecture.str.contains(r"(?:^|-)WD40(?:-|$)")
    prot["has_TPR"] = prot.architecture.str.contains(r"(?:^|-)TPR(?:-|$)")

    # --- cluster map
    m2f, fam_members = {}, defaultdict(list)
    with open(P / "cluster/tier1_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            m2f[mem] = rep
            fam_members[rep].append(mem)
    prot["family"] = prot.target.map(m2f)
    print("classified proteins not in cluster map:", prot.family.isna().sum(), file=sys.stderr)
    prot.to_csv(A / "het_protein_domains.tsv.gz", sep="\t", index=False, compression="gzip")

    HETCLS = ["NLR_NOD", "NLR_like_AAA", "NLR_accessory_only", "HET_domain", "HeLo_SesA_Goodbye_HETs", "Het-C"]
    hp = prot[prot.het_class.isin(HETCLS)].copy()
    print(hp.het_class.value_counts().to_string(), file=sys.stderr)

    ss = pd.read_csv(P / "samplesheet.with_clades.csv")
    spmap = {"Coccidioides immitis": "Ci", "Coccidioides posadasii": "Cp"}
    species = dict(zip(ss.Short, ss.Species.map(spmap)))
    strains = list(ss.Short)
    freq = pd.read_csv(P / "frequency_table.tsv", sep="\t").set_index("family")

    fam_set = sorted(set(hp.family.dropna()))
    cpos = {"NLR_NOD": 0, "NLR_like_AAA": 1, "NLR_accessory_only": 2, "HET_domain": 3, "HeLo_SesA_Goodbye_HETs": 4, "Het-C": 5}
    pc = prot.set_index("target")
    census = []
    for fam in fam_set:
        mem = fam_members[fam]
        cl = [pc.het_class[m] for m in mem if m in pc.index and pc.het_class[m] in cpos]
        ar = [pc.architecture[m] for m in mem if m in pc.index and pc.het_class[m] in cpos]
        cc = Counter(cl)
        # family class = highest-priority class carried by >= 20% of classified members, else modal
        fam_cls = sorted(cc, key=lambda k: (cc[k] / len(cl) < 0.2, cpos[k]))[0]
        per_strain = Counter(m.split("|")[0] for m in mem)
        census.append(dict(family=fam, family_class=fam_cls, n_members=len(mem), n_members_classified=len(cl),
                           frac_members_classified=round(len(cl) / len(mem), 3),
                           class_counts=";".join(f"{k}:{v}" for k, v in cc.most_common()),
                           modal_architecture=Counter(ar).most_common(1)[0][0],
                           architectures=";".join(f"{k}:{v}" for k, v in Counter(ar).most_common(4)),
                           rep_architecture=pc.architecture.get(fam, "-") if fam in pc.index else "-",
                           median_len=int(pd.Series([pc.tlen[m] for m in mem if m in pc.index]).median()),
                           n_strains_with_member=len(per_strain), max_copies_per_strain=max(per_strain.values()),
                           n_strains_multicopy=sum(1 for v in per_strain.values() if v > 1),
                           bin=freq.bin.get(fam, "NA"), strain_count_Ci_rep=freq.strain_count.get(fam, None)))
    census = pd.DataFrame(census)

    # --- presence matrix rows
    pm_rows = {}
    fs = set(fam_set)
    with open(P / "presence_matrix.rescued.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fam = line[:line.index("\t")]
            if fam in fs:
                pm_rows[fam] = line.rstrip("\n").split("\t")[1:]
    cols = header[1:]
    pm = pd.DataFrame.from_dict(pm_rows, orient="index", columns=cols)
    assert set(cols) == set(strains)
    pres = []
    for fam in census.family:
        row = pm.loc[fam]
        d = dict(family=fam)
        for lab in ("Ci", "Cp"):
            sub = row[[s for s in cols if species[s] == lab]]
            c = Counter(sub)
            d[f"{lab}_present"] = c.get("present", 0)
            d[f"{lab}_genome_only"] = c.get("genome_only", 0)
            d[f"{lab}_n"] = len(sub)
        for mode, keys in (("any", ("present", "genome_only")), ("annot", ("present",))):
            ci = sum(d[f"Ci_{k}"] for k in keys)
            cp = sum(d[f"Cp_{k}"] for k in keys)
            orr, p = fisher_exact([[ci, d["Ci_n"] - ci], [cp, d["Cp_n"] - cp]])
            d[f"Ci_frac_{mode}"] = round(ci / d["Ci_n"], 4)
            d[f"Cp_frac_{mode}"] = round(cp / d["Cp_n"], 4)
            d[f"fisher_p_{mode}"] = p
        pres.append(d)
    pres = pd.DataFrame(pres)
    for mode in ("any", "annot"):
        pres[f"bh_q_{mode}"] = multipletests(pres[f"fisher_p_{mode}"], method="fdr_bh")[1]
    census = census.merge(pres, on="family")
    inr = lambda x: (x >= 0.05) & (x <= 0.95)
    census["poly_within_annot"] = inr(census.Ci_frac_annot) | inr(census.Cp_frac_annot)
    census["poly_between_annot"] = census.bh_q_annot < 0.05
    census["poly_within_any"] = inr(census.Ci_frac_any) | inr(census.Cp_frac_any)
    census["poly_between_any"] = census.bh_q_any < 0.05
    for lab in ("Ci", "Cp"):
        cn = []
        for fam in census.family:
            ps = Counter(m.split("|")[0] for m in fam_members[fam])
            v = [ps[s] for s in strains if species[s] == lab and ps[s] > 0]
            cn.append(round(sum(v) / len(v), 2) if v else 0)
        census[f"{lab}_mean_copies_carriers"] = cn
    census = census.sort_values(["family_class", "n_members"], ascending=[True, False]).reset_index(drop=True)
    census.insert(0, "label", [f"H{i+1:03d}" for i in range(len(census))])
    census.to_csv(A / "het_family_census.tsv", sep="\t", index=False)
    pm.loc[census.family].to_csv(A / "het_presence_matrix.tsv.gz", sep="\t", compression="gzip")
    lab = dict(zip(census.family, census.label))

    # --- positions (reuse hrmA contig lengths and telomere-end table)
    gp = zread(P / "gene_positions.tsv.zst")
    gp["member"] = gp.Short + "|" + gp.protein_id
    gp["family"] = gp.member.map(m2f)
    clen = pd.read_csv(H / "contig_lengths.tsv.gz", sep="\t", header=None, names=["Short", "contig", "contig_len"])
    gp = gp.merge(clen, on=["Short", "contig"], how="left")
    gp["dist_end"] = pd.concat([gp.start - 1, gp.contig_len - gp.end], axis=1).min(axis=1)
    gp = gp.sort_values(["Short", "contig", "start"]).reset_index(drop=True)
    gp["idx"] = gp.groupby(["Short", "contig"]).cumcount()
    gp["n_on_contig"] = gp.groupby(["Short", "contig"]).idx.transform("size")
    gp["bin"] = gp.family.map(freq.bin)
    tel = pd.read_csv(H / "contig_telomere_ends.tsv.gz", sep="\t")
    gp = gp.merge(tel[["Short", "contig", "left_telo", "right_telo"]], on=["Short", "contig"], how="left")
    gp["nearest_end_telo"] = gp.left_telo.where((gp.start - 1) <= (gp.contig_len - gp.end), gp.right_telo)
    gp["genes_to_end"] = pd.concat([gp.idx, gp.n_on_contig - 1 - gp.idx], axis=1).min(axis=1)

    hset = set(hp.target)
    hg = gp[gp.member.isin(hset)].copy()
    hg = hg.merge(hp[["target", "het_class", "architecture", "tlen"]].rename(columns={"target": "member"}), on="member")
    hg["species"] = hg.Short.map(species)
    hg["label"] = hg.family.map(lab)
    hg["family_class"] = hg.family.map(dict(zip(census.family, census.family_class)))
    hg.to_csv(A / "het_gene_locations.tsv.gz", sep="\t", index=False, compression="gzip")
    print("HET/NLR proteins:", len(hset), "with positions:", len(hg), file=sys.stderr)

    # all other classified proteins (repeat/effector-only) also saved for context
    pclass = dict(zip(prot.target, prot.het_class))
    parch = dict(zip(prot.target, prot.architecture))

    # --- neighbourhoods (vectorised over gp row numbers)
    gp["row"] = range(len(gp))
    rowof = dict(zip(gp.member, gp.row))
    parts = []
    for off in range(-FLANK, FLANK + 1):
        ok = (hg.idx + off >= 0) & (hg.idx + off < hg.n_on_contig)
        a = hg[ok]
        rr = a.member.map(rowof).values + off
        g = gp.iloc[rr]
        parts.append(pd.DataFrame(dict(anchor=a.member.values, anchor_family=a.family.values, Short=a.Short.values,
                                       species=a.species.values, contig=a.contig.values, offset=off,
                                       member=g.member.values, family=g.family.values, bin=g.bin.values,
                                       start=g.start.values, end=g.end.values)))
    neigh = pd.concat(parts).sort_values(["anchor", "offset"])
    neigh["het_class"] = neigh.member.map(pclass).fillna("")
    neigh["het_architecture"] = neigh.member.map(parch).fillna("")
    neigh.to_csv(A / "het_neighbourhoods.tsv.gz", sep="\t", index=False, compression="gzip")
    print("done", file=sys.stderr)


if __name__ == "__main__":
    main()
