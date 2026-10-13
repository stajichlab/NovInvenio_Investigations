#!/usr/bin/env python3
"""Select ingroup (F. oxysporum f. sp. lycopersici, FOL) and outgroup strains.

Inputs, reused from the FOXY vs FSSC study (same BFD inventory and mash run):
  ../fusarium_FOXY_vs_FSSC/selection/bfd_fusarium_inventory.tsv
  ../fusarium_FOXY_vs_FSSC/selection/bfd_fusarium_all_vs_all.mash_dist.tsv.zst

Ingroup rule
  1. inventory status=ok and SPECIES_IN contains exactly "f. sp. lycopersici"
     (radicis-lycopersici excluded).
  2. BUSCO genome (fungi_odb12) C >= 90%.
  3. Protein-count floor: n_proteins >= PROT_FLOOR_ABS (11,000). The
     alternative "FOXY median - 3 MAD" (FOXY = the 209 FOXY strains of the
     FOXY vs FSSC study) is printed but not applied: MAD is 328, so that
     floor (13,972) would drop Fol59 (13,800 proteins, N50 3.0 Mb, BUSCO
     99.6%), a good assembly. The floor's job is to catch failed
     annotations (e.g. 416 proteins); none of the 28 FOL is below 13,800.
  4. Clone dedup: single-linkage groups at mash < DEDUP_D (0.0005); keep the
     member with the highest N50 (ties: BUSCO C, then ASMID).

Outgroup rule (4 strains)
  Candidates: FOSC-named strains (same name list as the FOXY study), not
  lycopersici / radicis-lycopersici, status ok, BUSCO >= 90, protein floor,
  N50 >= 1 Mb, N50 not within 0.1% of Fol4287's 4,589,937 (several BFD
  assemblies share that N50, a sign of scaffolding onto Fol4287), and mash
  distance to EVERY FOL strain >= OUT_MIN_D (0.010).
  Pick the 3 candidates nearest to FOL (min distance to any FOL strain) with a
  named forma specialis, one per f. sp.; plus the F. odoratissimum candidate
  with the highest N50 as a more distant FOSC outgroup.

Outputs: selection/fol_selected.tsv, selection/fol_selection_log.tsv,
selection/fol_mash_structure.tsv (pairwise mash among the FOL QC-passed set).
"""
from __future__ import annotations

import csv
import itertools
import re
import statistics
import subprocess
import sys
from pathlib import Path

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOL")
X = STUDY.parent / "fusarium_FOXY_vs_FSSC" / "selection"
BUSCO_DIR = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/results/genome_stats/BUSCO_genome")

MIN_BUSCO_C = 90.0
PROT_FLOOR_ABS = 11000
DEDUP_D = 0.0005
OUT_MIN_D = 0.010
N_OUT_FSP = 3
FOL4287_N50 = 4589937
FOSC = {"Fusarium oxysporum", "Fusarium odoratissimum", "Fusarium tardichlamydosporum",
        "Fusarium purpurascens", "Fusarium phialophorum", "Fusarium nirenbergiae",
        "Fusarium veterinarium"}


def key(p: str) -> str:
    return p.rsplit("/", 1)[-1].replace(".fa.gz", "")


def busco_c(asmid: str) -> float | None:
    hits = list(BUSCO_DIR.glob(f"*/{asmid}.BUSCO_summary.*.txt"))
    if not hits:
        return None
    m = re.search(r"C:([\d.]+)%", hits[0].read_text())
    return float(m.group(1)) if m else None


def nprot(path: str) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def fsp(s: str) -> str:
    m = re.search(r"f\. ?sp\. ?([A-Za-z-]+)", s or "")
    return m.group(1) if m else ""


def main() -> int:
    inv = {r["asmid"]: r for r in csv.DictReader(open(X / "bfd_fusarium_inventory.tsv"), delimiter="\t")
           if r["status"] == "ok"}
    fosc = {a for a, r in inv.items() if r["species"] in FOSC}
    d = {}
    p = subprocess.Popen(["zstd", "-dc", str(X / "bfd_fusarium_all_vs_all.mash_dist.tsv.zst")],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        a, b, x, _, _ = line.split("\t")
        a, b = key(a), key(b)
        if a in fosc and b in fosc:
            d[(a, b)] = float(x)
    p.wait()

    foxy_counts = [int(r["n_proteins"]) for r in csv.DictReader(open(X / "protein_counts.tsv"), delimiter="\t")
                   if r["group"] == "FOXY"]
    med = statistics.median(foxy_counts)
    mad = statistics.median(abs(x - med) for x in foxy_counts)
    floor = PROT_FLOOR_ABS
    print(f"FOXY protein count median={med} MAD={mad} (median-3MAD={med - 3 * mad}, not applied); floor={floor}", file=sys.stderr)

    log = []
    fol_all = sorted(a for a, r in inv.items()
                     if "f. sp. lycopersici" in r["species_in"] and "radicis" not in r["species_in"])
    info = {}
    for a in fol_all:
        r = inv[a]
        info[a] = {"busco": busco_c(a), "nprot": nprot(r["proteins"]), "n50": int(r["n50"] or 0)}
    qc = []
    for a in fol_all:
        i, r = info[a], inv[a]
        if i["busco"] is None or i["busco"] < MIN_BUSCO_C:
            log.append((r["dirname"], a, "IN", "busco", f"C={i['busco']}"))
        elif i["nprot"] < floor:
            log.append((r["dirname"], a, "IN", "protein_floor", f"n_proteins={i['nprot']} < {floor}"))
        else:
            qc.append(a)

    with open(STUDY / "selection" / "fol_mash_structure.tsv", "w") as fh:
        fh.write("strain_a\tstrain_b\tmash_dist\n")
        for a, b in itertools.combinations(qc, 2):
            fh.write(f"{inv[a]['dirname']}\t{inv[b]['dirname']}\t{d[(a, b)]:.6f}\n")

    parent = {a: a for a in qc}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in itertools.combinations(qc, 2):
        if d[(a, b)] < DEDUP_D:
            parent[find(a)] = find(b)
    groups = {}
    for a in qc:
        groups.setdefault(find(a), []).append(a)
    fol = []
    for g in groups.values():
        g.sort(key=lambda a: (-info[a]["n50"], -(info[a]["busco"] or 0), a))
        fol.append(g[0])
        for a in g[1:]:
            log.append((inv[a]["dirname"], a, "IN", "clone_dedup",
                        f"mash={d[(a, g[0])]:.6f} to kept {inv[g[0]]['dirname']}"))
    fol.sort()

    # outgroup
    fol_set_all = set(fol_all)
    cands = []
    for a in fosc:
        r = inv[a]
        if a in fol_set_all or "lycopersici" in r["species_in"]:
            continue
        n50 = int(r["n50"] or 0)
        if n50 < 1_000_000 or abs(n50 - FOL4287_N50) <= FOL4287_N50 * 0.001:
            continue
        mind = min(d[(a, b)] for b in fol_all)
        if mind < OUT_MIN_D:
            continue
        cands.append((mind, a))
    cands.sort()
    out, used_fsp = [], set()
    for mind, a in cands:
        f = fsp(inv[a]["species_in"])
        if not f or f in used_fsp or inv[a]["species"] != "Fusarium oxysporum":
            continue
        bc, npr = busco_c(a), nprot(inv[a]["proteins"])
        if bc is None or bc < MIN_BUSCO_C or npr < floor:
            log.append((inv[a]["dirname"], a, "OUT", "candidate_failed_qc", f"BUSCO={bc} n_proteins={npr}"))
            continue
        info[a] = {"busco": bc, "nprot": npr, "n50": int(inv[a]["n50"]), "min_d_fol": mind}
        out.append(a)
        used_fsp.add(f)
        if len(out) == N_OUT_FSP:
            break
    odor = sorted((a for m, a in cands if inv[a]["species"] == "Fusarium odoratissimum"),
                  key=lambda a: -int(inv[a]["n50"]))
    for a in odor:
        bc, npr = busco_c(a), nprot(inv[a]["proteins"])
        if bc is not None and bc >= MIN_BUSCO_C and npr >= floor:
            info[a] = {"busco": bc, "nprot": npr, "n50": int(inv[a]["n50"]),
                       "min_d_fol": min(d[(a, b)] for b in fol_all)}
            out.append(a)
            break
        log.append((inv[a]["dirname"], a, "OUT", "candidate_failed_qc", f"BUSCO={bc} n_proteins={npr}"))

    cols = ["group", "dirname", "asmid", "locustag", "species", "species_in", "strain", "ncbi_taxid",
            "n50", "total_len", "n_contigs", "busco_genome_C", "n_proteins", "min_mash_to_FOL",
            "proteins", "gff3", "scaffolds"]
    with open(STUDY / "selection" / "fol_selected.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for grp, ids in (("IN", fol), ("OUT", out)):
            for a in ids:
                r, i = inv[a], info[a]
                mind = i.get("min_d_fol", min((d[(a, b)] for b in fol_all if b != a), default=0))
                w.writerow([grp, r["dirname"], a, r["locustag"], r["species"], r["species_in"], r["strain"],
                            r["ncbi_taxid"], r["n50"], r["total_len"], r["n_contigs"], i["busco"],
                            i["nprot"], f"{mind:.4f}", r["proteins"], r["gff3"], r["scaffolds"]])
    with open(STUDY / "selection" / "fol_selection_log.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dirname", "asmid", "group", "step", "reason"])
        w.writerows(log)
    with open(STUDY / "selection" / "fol_candidate_qc.tsv", "w") as fh:
        fh.write("dirname\tasmid\tbusco_genome_C\tn_proteins\tn50\n")
        for a in fol_all:
            fh.write(f"{inv[a]['dirname']}\t{a}\t{info[a]['busco']}\t{info[a]['nprot']}\t{info[a]['n50']}\n")
    print(f"FOL: {len(fol_all)} labelled, {len(qc)} pass QC, {len(fol)} after clone dedup; "
          f"outgroup {len(out)}: {[inv[a]['dirname'] for a in out]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
