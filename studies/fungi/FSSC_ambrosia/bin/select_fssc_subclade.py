#!/usr/bin/env python3
"""Select ingroup + outgroup strains for one FSSC sub-clade pangenome study.

Usage:
  /usr/bin/python3.12 select_fssc_subclade.py --study solani   --study-dir <.../FSSC_solani>
  /usr/bin/python3.12 select_fssc_subclade.py --study ambrosia --study-dir <.../FSSC_ambrosia>

The same file is copied into studies/fungi/FSSC_solani/bin and
studies/fungi/FSSC_ambrosia/bin. It reuses, read only, the BFD Fusarium
inventory, the all-vs-all mash matrix and the complex-anchor table built in
studies/fungi/fusarium_FOXY_vs_FSSC/selection/ (see that study's README).

Pool: every inventory row with status=ok whose nearest complex anchor is
F. vanettenii 77-13-4 (FSSC) at mash d <= 0.20 (nearest_complex_anchor.tsv).

Ingroup candidates
  solani:   the single-linkage cluster at mash d < SOLANI_LINK (0.03) that
            contains F. solani CBS 140079 (GCA_045529305.1). Membership is by
            mash only, whatever the BFD species label.
  ambrosia: Ambrosia Fusarium Clade (AFC). Candidates = AFC_NAMES, the
            unnamed "Neocosmospora sp. B6999-*" strains, the "Fusarium cf.
            solani" strains that are not EtdFoc/EthFoc isolates, and the
            conditional names COND_NAMES (papillatum, duplospermum).
            Placement test, applied to every candidate label (one label =
            one species name; each cf. solani strain is its own label):
              keep iff  min d(label, AFC strains of other labels)
                      <  min d(label, strains of named non-AFC FSSC species)
            "Named non-AFC FSSC species" = pool strains with a binomial
            species name that is not an AFC/conditional name, not "sp.",
            not "cf.", not "Fusarium oxysporum" (one mislabelled pool strain,
            14003), and not in MASH_NESTED_IN_AFC (F. drepaniforme: mash
            0.024 to F. tuaranense, 0.089 to any other named species; it is
            not on the AFC name list, so it is neither ingroup nor comparator
            nor outgroup).

QC (both studies, logged):
  1. BUSCO genome (fungi_odb12) complete >= MIN_BUSCO (90%).
  2. Protein count >= MIN_PROTEINS (11,000) and >= median - 3 * MAD_SCALED,
     median and MAD over the BUSCO-passing ingroup candidates;
     MAD_SCALED = 1.4826 * raw MAD (a robust SD).
  3. Same strain deposited twice (strain ID equal after removing
     non-alphanumerics): keep the highest N50.
  4. Near-identical clones: single-linkage at mash d < CLONE_D (0.0005);
     keep the highest N50 per group (ties: BUSCO, then ASMID).

Outgroup: from named FSSC species that are not ingroup (for solani: not
labelled "Fusarium solani"; those strains outside the cluster are other
lineages). For ambrosia, species identity of F. solani is by mash: the solani
cluster above is "F. solani (s.s. cluster)", and solani-labelled strains
outside it are not used. Strains must pass BUSCO >= 90, >= 11,000
proteins and N50 >= OUT_MIN_N50 (100 kb). Species are ranked by min mash d to the final ingroup; the N_OUT (4)
nearest species each contribute their highest-N50 passing strain.

Outputs (selection/): selected_strains.tsv, selection_log.tsv,
candidates.tsv (every candidate with QC fields and decision).
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics as st
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

FOXY_SEL = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOXY_vs_FSSC/selection")
BUSCO_DIR = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/results/genome_stats/BUSCO_genome")

POOL_MAX_D = 0.20
SOLANI_LINK = 0.03
SOLANI_SEED = "GCA_045529305.1_TL_NS_CBS140079"
MIN_BUSCO = 90.0
MIN_PROTEINS = 11000
MAD_K = 3.0
CLONE_D = 0.0005
N_OUT = 4
OUT_MIN_N50 = 100_000

AFC_NAMES = {"euwallaceae", "ambrosium", "floridanum", "obliquiseptatum", "tuaranense", "kuroshium",
             "oligoseptatum", "pseudensiforme", "protoensiforme", "rekanum", "awan", "akasia", "mekan",
             "variasi", "warna"}
COND_NAMES = {"papillatum", "duplospermum"}
MASH_NESTED_IN_AFC = {"drepaniforme"}


def epithet(sp: str) -> str:
    return sp.split()[-1]


def label(r: dict) -> str:
    sp = r["species"]
    if sp == "Fusarium cf. solani":
        return f"cf. solani {r['strain']}"
    return sp


def is_named(sp: str) -> bool:
    w = sp.split()
    return len(w) == 2 and w[1] not in ("sp.",) and "cf." not in w


def load_busco() -> dict[str, float]:
    out = {}
    for shard in BUSCO_DIR.iterdir():
        for f in shard.iterdir():
            if ".BUSCO_summary." in f.name:
                m = re.search(r"C:([\d.]+)%", f.read_text())
                if m:
                    out[f.name.split(".BUSCO_summary.")[0]] = float(m.group(1))
    return out


def count_proteins(path: str) -> int:
    p = subprocess.run(["grep", "-c", ">", path], capture_output=True, text=True, check=True)
    return int(p.stdout.strip())


def single_linkage(ids, d, thr):
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    ids = sorted(ids)
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if d[(a, b)] < thr:
                parent[find(a)] = find(b)
    g = defaultdict(list)
    for i in ids:
        g[find(i)].append(i)
    return sorted((sorted(v) for v in g.values()), key=lambda v: v[0])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True, choices=["solani", "ambrosia"])
    ap.add_argument("--study-dir", required=True)
    args = ap.parse_args()
    sel = Path(args.study_dir) / "selection"
    sel.mkdir(parents=True, exist_ok=True)

    inv = {r["asmid"]: r for r in csv.DictReader(open(FOXY_SEL / "bfd_fusarium_inventory.tsv"), delimiter="\t")
           if r["status"] == "ok"}
    near = {r["asmid"]: r for r in csv.DictReader(open(FOXY_SEL / "nearest_complex_anchor.tsv"), delimiter="\t")}
    pool = {a for a, r in near.items() if a in inv and r["nearest"] == "FSSC" and float(r["d_FSSC"]) <= POOL_MAX_D}
    d = {}
    p = subprocess.Popen(["zstd", "-dc", str(FOXY_SEL / "bfd_fusarium_all_vs_all.mash_dist.tsv.zst")],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        a, b, x, _, _ = line.split("\t")
        a, b = a.removesuffix(".fa.gz"), b.removesuffix(".fa.gz")
        if a in pool and b in pool:
            d[(a, b)] = float(x)
    p.wait()
    busco = load_busco()
    R = {}
    for a in pool:
        r = dict(inv[a])
        r["busco_genome_C"] = busco.get(a)
        r["n_proteins"] = count_proteins(r["proteins"])
        r["n50_i"] = int(r["n50"] or 0)
        r["d_FSSC_anchor"] = near[a]["d_FSSC"]
        R[a] = r
    print(f"pool: {len(pool)} FSSC-placed strains", file=sys.stderr)

    log = []
    cand_rows = []

    # ---- ingroup candidates
    solani_cl = next(g for g in single_linkage(pool, d, SOLANI_LINK) if SOLANI_SEED in g)
    if args.study == "solani":
        cl = solani_cl
        cand = set(cl)
        print(f"solani cluster (d<{SOLANI_LINK}): {len(cl)} strains, max within "
              f"{max(d[(a, b)] for a in cl for b in cl):.4f}", file=sys.stderr)
        for a in pool - cand:
            if R[a]["species"] in ("Fusarium solani", "Fusarium cf. solani"):
                dm = min(d[(a, b)] for b in cl)
                log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "not_in_solani_cluster",
                            f"labelled {R[a]['species']}; min mash to cluster {dm:.4f}"))
        comparator = None
    else:
        labs = defaultdict(list)
        for a in pool:
            labs[label(R[a])].append(a)

        def is_cand(lab, a):
            sp = R[a]["species"]
            if sp.startswith("Neocosmospora sp. B6999"):
                return True
            if sp == "Fusarium cf. solani":
                return not R[a]["strain"].startswith(("EtdFoc", "EthFoc"))
            return sp.startswith("Fusarium ") and epithet(sp) in (AFC_NAMES | COND_NAMES)
        cand_labs = {lab for lab, v in labs.items() if is_cand(lab, v[0])}
        core_labs = {lab for lab in cand_labs if epithet(R[labs[lab][0]]["species"]) in AFC_NAMES
                     or R[labs[lab][0]]["species"].startswith("Neocosmospora sp. B6999")}
        comparator = [a for a in pool if is_named(R[a]["species"]) and label(R[a]) not in cand_labs
                      and epithet(R[a]["species"]) not in MASH_NESTED_IN_AFC
                      and R[a]["species"] != "Fusarium oxysporum"]
        cand = set()
        for lab in sorted(cand_labs):
            A = labs[lab]
            other_core = [b for k in core_labs - {lab} for b in labs[k]]
            da, ba = min((d[(a, b)], b) for a in A for b in other_core)
            dn, bn = min((d[(a, b)], b) for a in A for b in comparator)
            keep = da < dn
            msg = (f"min d to other AFC {da:.4f} ({R[ba]['species']} {R[ba]['strain']}); "
                   f"min d to named non-AFC {dn:.4f} ({R[bn]['species']} {R[bn]['strain']})")
            print(f"placement {lab:40s} {'KEEP' if keep else 'DROP'}  {msg}", file=sys.stderr)
            for a in A:
                R[a]["placement"] = msg
                if keep:
                    cand.add(a)
                else:
                    log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "afc_placement", msg))
        for a in pool:
            if epithet(R[a]["species"]) in MASH_NESTED_IN_AFC or (
                    not is_named(R[a]["species"]) and a not in cand and
                    min(d[(a, b)] for b in cand) < min(d[(a, b)] for b in comparator)):
                log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "afc_by_mash_not_by_name",
                            f"mash nearer AFC ({min(d[(a, b)] for b in cand):.4f}) than named non-AFC "
                            f"({min(d[(a, b)] for b in comparator):.4f}); not on AFC name list, not used"))

    # ---- QC
    busco_ok = []
    for a in sorted(cand):
        bc = R[a]["busco_genome_C"]
        if bc is None or bc < MIN_BUSCO:
            log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "busco", f"BUSCO genome C={bc}"))
            R[a]["decision"] = "drop_busco"
        else:
            busco_ok.append(a)
    np_ = [R[a]["n_proteins"] for a in busco_ok]
    med = st.median(np_)
    mad = 1.4826 * st.median(abs(x - med) for x in np_)
    floor = max(MIN_PROTEINS, med - MAD_K * mad)
    print(f"protein counts (BUSCO-passing candidates n={len(np_)}): median {med:.0f}, scaled MAD {mad:.0f}, "
          f"median-3MAD {med - MAD_K * mad:.0f}, floor used {floor:.0f}", file=sys.stderr)
    qc_ok = []
    for a in busco_ok:
        n = R[a]["n_proteins"]
        if n < floor:
            why = f"{n} proteins < {'11,000' if n < MIN_PROTEINS else f'median-3MAD {med - MAD_K * mad:.0f}'}"
            log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "protein_count", why))
            R[a]["decision"] = "drop_proteins"
        else:
            qc_ok.append(a)
    # same strain deposited twice (e.g. "NRRL 20438" and "NRRL20438"): keep highest N50
    by_strain = defaultdict(list)
    for a in qc_ok:
        by_strain[(R[a]["species"], re.sub(r"[^A-Za-z0-9]", "", R[a]["strain"]).upper())].append(a)
    qc_ok = []
    for g in by_strain.values():
        g.sort(key=lambda a: (-R[a]["n50_i"], -(R[a]["busco_genome_C"] or 0), a))
        qc_ok.append(g[0])
        for a in g[1:]:
            log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "same_strain",
                        f"second assembly of strain {R[g[0]]['strain']} (mash {d[(a, g[0])]:.5f}); kept {g[0]} "
                        f"(N50 {R[g[0]]['n50_i']} vs {R[a]['n50_i']})"))
            R[a]["decision"] = "drop_same_strain"
    ingroup = []
    clone_groups = []
    for g in single_linkage(qc_ok, d, CLONE_D):
        g.sort(key=lambda a: (-R[a]["n50_i"], -(R[a]["busco_genome_C"] or 0), a))
        ingroup.append(g[0])
        R[g[0]]["decision"] = "ingroup"
        if len(g) > 1:
            clone_groups.append(g)
        for a in g[1:]:
            log.append((R[a]["dirname"], a, R[a]["species"], R[a]["strain"], "clone",
                        f"mash {d[(a, g[0])]:.5f} to kept {R[g[0]]['species']} {R[g[0]]['strain']} "
                        f"(N50 {R[a]['n50_i']} vs {R[g[0]]['n50_i']})"))
            R[a]["decision"] = "drop_clone"
    ingroup.sort()
    print(f"ingroup: {len(cand)} candidates -> {len(busco_ok)} BUSCO ok -> {len(qc_ok)} protein ok -> "
          f"{len(ingroup)} after clone removal", file=sys.stderr)

    # ---- outgroup
    if args.study == "solani":
        out_pool = [a for a in pool if a not in cand and is_named(R[a]["species"])
                    and R[a]["species"] not in ("Fusarium solani",)]
    else:
        # Outgroup species by mash, not by label: members of the F. solani s.s.
        # cluster (d < 0.03 to CBS 140079; includes "F. haematococcum" LQ1 and
        # "Fusarium sp." R1) count as F. solani; "F. solani"-labelled strains
        # outside that cluster (other, unnamed FSSC lineages) are not used.
        out_pool = []
        for a in comparator:
            if a in solani_cl:
                R[a]["species_out"] = "Fusarium solani (s.s. cluster)"
            elif R[a]["species"] == "Fusarium solani":
                continue
            out_pool.append(a)
        for a in solani_cl:
            if a not in out_pool:
                R[a]["species_out"] = "Fusarium solani (s.s. cluster)"
                out_pool.append(a)
    out_pool = [a for a in out_pool if (R[a]["busco_genome_C"] or 0) >= MIN_BUSCO
                and R[a]["n_proteins"] >= MIN_PROTEINS and R[a]["n50_i"] >= OUT_MIN_N50]
    by_sp = defaultdict(list)
    for a in out_pool:
        by_sp[R[a].get("species_out", R[a]["species"])].append(a)
    sp_d = {sp: min(d[(a, b)] for a in v for b in ingroup) for sp, v in by_sp.items()}
    ranked = sorted(sp_d, key=sp_d.get)
    print("nearest named non-ingroup species (QC-passing strains):", file=sys.stderr)
    for sp in ranked[:10]:
        print(f"  {sp:35s} {sp_d[sp]:.4f}  n={len(by_sp[sp])}", file=sys.stderr)
    outgroup = []
    for sp in ranked[:N_OUT]:
        best = sorted(by_sp[sp], key=lambda a: (-R[a]["n50_i"], -(R[a]["busco_genome_C"] or 0), a))[0]
        outgroup.append(best)
        R[best]["decision"] = "outgroup"
        R[best]["out_min_d"] = min(d[(best, b)] for b in ingroup)

    cols = ["GROUP", "dirname", "asmid", "locustag", "species", "species_in", "strain", "ncbi_taxid",
            "bioproject", "n50", "total_len", "n_contigs", "busco_genome_C", "n_proteins", "d_FSSC_anchor",
            "min_d_to_ingroup", "proteins", "gff3", "scaffolds"]
    with open(sel / "selected_strains.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for grp, ids in (("IN", ingroup), ("OUT", outgroup)):
            for a in ids:
                r = R[a]
                r["GROUP"] = grp
                r["min_d_to_ingroup"] = (f"{min(d[(a, b)] for b in ingroup if b != a):.4f}")
                w.writerow([r.get(c, "") for c in cols])
    with open(sel / "selection_log.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dirname", "asmid", "species", "strain", "step", "reason"])
        w.writerows(sorted(log, key=lambda x: (x[4], x[2], x[3])))
    with open(sel / "candidates.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["asmid", "dirname", "species", "strain", "busco_genome_C", "n_proteins", "n50", "decision",
                    "placement"])
        for a in sorted(cand | set(outgroup), key=lambda a: (R[a]["species"], R[a]["strain"])):
            r = R[a]
            w.writerow([a, r["dirname"], r["species"], r["strain"], r["busco_genome_C"], r["n_proteins"],
                        r["n50"], r.get("decision", "drop_placement"), r.get("placement", "")])
    with open(sel / "clone_groups.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["kept_asmid", "kept", "dropped"])
        for g in clone_groups:
            w.writerow([g[0], f"{R[g[0]]['species']} {R[g[0]]['strain']}",
                        "; ".join(f"{R[a]['species']} {R[a]['strain']}" for a in g[1:])])
    # full distance matrix of the selected strains
    ids = ingroup + outgroup
    with open(sel / "selected_mash_matrix.tsv", "w") as fh:
        fh.write("asmid\t" + "\t".join(ids) + "\n")
        for a in ids:
            fh.write(a + "\t" + "\t".join("0" if a == b else f"{d[(a, b)]:.5f}" for b in ids) + "\n")
    print(f"wrote {len(ingroup)} IN + {len(outgroup)} OUT to {sel / 'selected_strains.tsv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
