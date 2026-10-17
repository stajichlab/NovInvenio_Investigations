#!/usr/bin/env python3
"""Shared-domain context, stage 1: gene-level tables.

Input: hmmsearch/all_proteomes.shared_domains.domtbl.idx.zst (hmmsearch --cut_ga
of 13 Pfam 38.2 models over every proteome in hmmsearch/manifest.tsv).

Writes (tables/):
  genes.afum.tsv.gz      one row per (protein, Pfam model): A. fumigatus study (295 proteomes)
  genes.cocci.tsv.gz     same, Coccidioides study (529 proteomes)
  genes.survey.tsv.gz    same, survey proteomes (no family/position context)
  genome_counts.tsv.gz   one row per proteome x label: number of proteins with the label
  reference_gene_order.tsv.gz   all genes of Asfu_Af293, RS, UCSF_Cp_Silveira in contig order
                                 (Af293 from full_v072 positions + unplaced Swiss-Prot coords;
                                 RS/Silveira from v5 positions)

Family context per run: cluster/tier1_cluster.tsv (member -> family),
frequency_table.tsv (bin, strain_count), significant_islands.tsv (member_families).
Positions: gene_positions.tsv.zst (Af: full_v072; Coccidioides: v5, v2 when absent in v5).
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
D = CP / "analysis/shared_domains_context"
HD = D / "hmmsearch"
T = D / "tables"
RUNS = {
    "full_v072": AF / "results/full_v072/config/pangenome",
    "immitis_in_posadasii_out_v5": CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome",
    "posadasii_in_immitis_out_v2": CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome",
}
RUN_TAG = {"full_v072": "af", "immitis_in_posadasii_out_v5": "v5", "posadasii_in_immitis_out_v2": "v2"}
SET_RUNS = {"afum": ["full_v072"], "cocci": ["immitis_in_posadasii_out_v5", "posadasii_in_immitis_out_v2"]}
EXTRA_AF293 = AF / "analysis/starship_island_overlap/inputs/af293_unplaced_swissprot_coords.tsv"
REFS = {"afum": ["Asfu_Af293"], "cocci": ["RS", "UCSF_Cp_Silveira"]}
# analysis label per Pfam model; the three HrmA models share one label
LABEL = {"AFUB_07903_YDR124W_hel": "PF11001", "HrmA": "HrmA", "HrmA_middle": "HrmA", "HrmA_N": "HrmA"}


def zread(path):
    path = Path(path)
    if path.suffix == ".zst":
        return subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True).stdout
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


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


def main():
    man = {}
    with open(HD / "manifest.tsv", newline="") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t"), 1):
            man[str(i)] = r
    # hits[(idx, pid)][model] = list of domain dicts
    hits = defaultdict(lambda: defaultdict(list))
    acc = {}
    for line in zread(HD / "all_proteomes.shared_domains.domtbl.idx.zst"):
        p = line.rstrip("\n").split(None, 23)
        idx, pid, tlen, model, macc = p[0], p[1], int(p[3]), p[4], p[5]
        acc[model] = macc
        hits[(idx, pid)][model].append(dict(tlen=tlen, full_E=p[7], full_score=float(p[8]),
                                            i_E=p[13], dom_score=float(p[14]),
                                            hmm_from=int(p[16]), hmm_to=int(p[17]),
                                            env_from=int(p[20]), env_to=int(p[21])))
    print(f"hit proteins: {len(hits)}", file=sys.stderr)

    # family context
    ctx = {}
    for run, pg in RUNS.items():
        m2r = {}
        with open(pg / "cluster/tier1_cluster.tsv") as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")
                m2r[mem] = rep
        freq = {}
        with open(pg / "frequency_table.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                freq[r["family"]] = (r["bin"], r["strain_count"])
        fam_isl = defaultdict(list)
        with open(pg / "significant_islands.tsv", newline="") as fh:
            for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
                for f in r["member_families"].split(","):
                    fam_isl[f].append((i, int(r["n_strains"])))
        ctx[run] = dict(m2r=m2r, freq=freq, fam_isl=fam_isl)

    needed_short = {man[i]["short"] for (i, _) in hits if man[i]["set"] in ("afum", "cocci")}
    needed_short |= {"Asfu_Af293", "RS", "UCSF_Cp_Silveira"}
    extra = {}
    with open(EXTRA_AF293, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            extra[r["protein_id"]] = (r["contig"], int(r["start"]), int(r["end"]))

    def positions(run, shorts):
        bys = defaultdict(list)
        for line in zread(RUNS[run] / "gene_positions.tsv.zst"):
            p = line.rstrip("\n").split("\t")
            if p[0] == "Short" or p[0] not in shorts:
                continue
            bys[(p[0], p[2])].append((int(p[3]), int(p[4]), p[1]))
        return bys

    rank = {}            # (short, pid) -> (contig, start, end, rank, source_run)
    ref_rows = []
    bys_af = positions("full_v072", needed_short)
    placed = {pid for (s, c), lst in bys_af.items() if s == "Asfu_Af293" for (_, _, pid) in lst}
    n_extra = 0
    for pid, (c, s, e) in extra.items():
        if pid not in placed:
            bys_af[("Asfu_Af293", c)].append((s, e, pid))
            n_extra += 1
    print(f"Af293 unplaced coords added: {n_extra}", file=sys.stderr)
    for run, bys in (("full_v072", bys_af),
                     ("immitis_in_posadasii_out_v5", positions("immitis_in_posadasii_out_v5", needed_short)),
                     ("posadasii_in_immitis_out_v2", positions("posadasii_in_immitis_out_v2", needed_short))):
        for (short, c), lst in bys.items():
            lst.sort()
            for i, (s, e, pid) in enumerate(lst):
                if (short, pid) not in rank:
                    rank[(short, pid)] = (c, s, e, i, RUN_TAG[run])
                    if (short == "Asfu_Af293" and run == "full_v072") or (
                            short in ("RS", "UCSF_Cp_Silveira") and run == "immitis_in_posadasii_out_v5"):
                        ref_rows.append(dict(short=short, contig=c, start=s, end=e, rank=i,
                                             n_genes_contig=len(lst), protein_id=pid, positions_run=RUN_TAG[run]))
        del bys
    write_tsv(T / "reference_gene_order.tsv.gz", ref_rows)

    # gene rows
    rows = defaultdict(list)
    counts = defaultdict(lambda: defaultdict(int))
    for (idx, pid), models in sorted(hits.items(), key=lambda x: (int(x[0][0]), x[0][1])):
        r = man[idx]
        s = r["set"]
        short = r["short"]
        arch = sorted(((d["env_from"], d["env_to"], m) for m, ds in models.items() for d in ds))
        arch_s = "|".join(f"{m}:{a}-{b}" for a, b, m in arch)
        labels = sorted({LABEL.get(m, m) for m in models})
        for lab in labels:
            counts[idx][lab] += 1
        loc = rank.get((short, pid))
        base = dict(set=s, group=r["group"], species=r["species"], short=short, protein_id=pid,
                    protein_length=next(iter(models.values()))[0]["tlen"],
                    protein_labels=",".join(labels), protein_architecture=arch_s)
        if s != "survey":
            if loc:
                base.update(contig=loc[0], start=loc[1], end=loc[2], gene_rank_on_contig=loc[3],
                            positions_run=loc[4])
            else:
                base.update(contig="NA", start="", end="", gene_rank_on_contig="", positions_run="")
            key = f"{short}|{pid}"
            for run in SET_RUNS[s]:
                t = RUN_TAG[run]
                cx = ctx[run]
                fam = cx["m2r"].get(key)
                if fam is None:
                    base[f"{t}_family"] = "not_in_run"
                    continue
                b, sc = cx["freq"].get(fam, ("NA", "NA"))
                isl = cx["fam_isl"].get(fam, [])
                base[f"{t}_family"] = fam
                base[f"{t}_family_bin"] = b
                base[f"{t}_family_strain_count"] = sc
                base[f"{t}_family_n_islands"] = len(isl)
                base[f"{t}_family_n_multi_strain_islands"] = sum(n >= 2 for _, n in isl)
        else:
            base["study"] = r["study"]
            base["taxon_group"] = r["taxon_group"]
        for m, ds in models.items():
            best = max(ds, key=lambda d: d["dom_score"])
            row = dict(base)
            row.update(domain=m, domain_acc=acc[m], label=LABEL.get(m, m), n_domain_copies=len(ds),
                       env=";".join(f"{d['env_from']}-{d['env_to']}" for d in ds),
                       hmm=";".join(f"{d['hmm_from']}-{d['hmm_to']}" for d in ds),
                       best_dom_score=best["dom_score"], best_i_Evalue=best["i_E"],
                       full_seq_score=best["full_score"], full_seq_Evalue=best["full_E"])
            rows[s].append(row)
    for s in ("afum", "cocci", "survey"):
        write_tsv(T / f"genes.{s}.tsv.gz", rows[s])
        print(f"{s}: rows={len(rows[s])} proteins={len({(r['short'], r['protein_id']) for r in rows[s]})}",
              file=sys.stderr)
    labels = sorted({LABEL.get(m, m) for m in acc})
    grows = []
    for idx, r in man.items():
        for lab in labels:
            grows.append(dict(set=r["set"], study=r["study"], group=r["group"], species=r["species"],
                              short=r["short"], taxon_group=r["taxon_group"], label=lab,
                              n_proteins=counts[idx][lab]))
    write_tsv(T / "genome_counts.tsv.gz", grows)


if __name__ == "__main__":
    main()
