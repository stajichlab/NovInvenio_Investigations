#!/usr/bin/env python3
"""Select the FOXY and FSSC strain sets for the reciprocal pangenome runs.

Inputs (all under selection/):
  bfd_fusarium_inventory.tsv            (bin/inventory_bfd_fusarium.py)
  bfd_fusarium_all_vs_all.mash_dist.tsv.zst
      mash dist (k=21, s=10000; BFD's own results/mash_sketch/<ASMID>.fa.msh)
      over every inventory row with status=ok.

Rule, in order (every drop is written to selection/selection_log.tsv):
  1. status=ok in the inventory (proteins + gff3 + genome present, LOCUSTAG
     joins to exactly one samples.csv row).
  2. BUSCO (genome mode, fungi_odb12, BFD results/genome_stats/BUSCO_genome)
     complete >= MIN_BUSCO_C.
  3. Complex membership = taxonomic name AND mash placement:
       FOXY: SPECIES in FOXY_NAMES, nearest anchor = Fol4287
             (GCF_000149955.1) and mash distance to it <= FOXY_MAX_D.
       FSSC: SPECIES in FSSC_NAMES or genus Neocosmospora, nearest anchor =
             F. vanettenii 77-13-4 (GCF_000151355.1), distance <= FSSC_MAX_D.
     Anchors for the other complexes are in ANCHORS; "nearest" is over all.
  4. Near-clone dedup inside each complex: single-linkage groups at mash
     distance < DEDUP_D (the pipeline's own pangenome_mash_threshold). Keep
     the member with the highest N50 (ties: BUSCO C, then ASMID).
  5. FOXY only: stratified subsample to FOXY_TARGET. Strata = UPGMA clusters
     of the deduplicated FOXY strains cut at STRATUM_H. Each stratum gets
     max(1, floor(target * n_k / N)); the remainder goes to the largest
     fractional parts. Inside a stratum, strains are drawn with
     random.Random(SEED): first one per distinct label (forma specialis from
     SPECIES_IN for F. oxysporum; the species name for the FOSC segregate
     species), then at random.
     Then any FOXY species name with no strain drawn gets its highest-N50
     deduplicated strain added.
  FSSC is kept whole after steps 1-4.

Outputs: selection/selected_strains.tsv, selection/selection_log.tsv,
selection/foxy_strata.tsv. Run with the NovInvenio pixi python (needs scipy).
"""
from __future__ import annotations

import csv
import math
import random
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOXY_vs_FSSC")
SEL = STUDY / "selection"
BUSCO_DIR = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/results/genome_stats/BUSCO_genome")

MIN_BUSCO_C = 90.0
FOXY_MAX_D = 0.035
FSSC_MAX_D = 0.17
DEDUP_D = 0.001
FOXY_TARGET = 200
STRATUM_H = 0.01
SEED = 20260925

ANCHORS = {
    "FOSC": "GCF_000149955.1_ASM14995v2",       # F. oxysporum f. sp. lycopersici 4287
    "FSSC": "GCF_000151355.1_v2.0",             # F. vanettenii 77-13-4
    "FSAMSC": "GCA_037255895.1_ASM3725589v1",   # F. graminearum KSU23473
    "FFSC": "GCA_003316975.2_ASM331697v2",      # F. verticillioides
    "FIESC": "GCA_002894235.1_CS5819v1",        # F. equiseti
    "FLSC": "GCA_014898835.1_ASM1489883v1",     # F. lateritium
    "FTSC": "GCA_003045085.1_ASM304508v1",      # F. tricinctum
    "FRSC": "GCA_014899085.1_ASM1489908v1",     # F. redolens
    "FNSC": "GCA_013623555.1_ASM1362355v1",     # F. nisikadoi
}

# F. oxysporum species complex: F. oxysporum plus the segregate species of
# Lombard et al. 2019 (Persoonia 43) and Maryani et al. 2019 (Stud. Mycol. 92)
# that occur in BFD. F. foetens, F. commune, F. redolens and F. inflexum are
# NOT included (outside FOSC s. str. in those treatments).
FOXY_NAMES = {
    "Fusarium oxysporum", "Fusarium odoratissimum", "Fusarium tardichlamydosporum",
    "Fusarium purpurascens", "Fusarium phialophorum", "Fusarium nirenbergiae",
    "Fusarium veterinarium",
}

# F. solani species complex = Neocosmospora sensu Sandoval-Denis et al. 2019
# (Persoonia 43) / FSSC of O'Donnell et al. 2020 and Geiser et al. 2021, as
# named in BFD under Fusarium. Names marked (m) are ones whose FSSC placement
# I could not confirm from the literature from memory; they are kept only
# because step 3 also requires mash placement next to F. vanettenii.
FSSC_NAMES = {
    "Fusarium solani", "Fusarium cf. solani", "Fusarium vanettenii",
    "Fusarium falciforme", "Fusarium cf. falciforme", "Fusarium keratoplasticum",
    "Fusarium haematococcum", "Fusarium solani-melongenae", "Fusarium metavorans",
    "Fusarium liriodendri", "Fusarium paranaense", "Fusarium neocosmosporiellum",
    "Fusarium euwallaceae", "Fusarium ambrosium", "Fusarium floridanum",
    "Fusarium kuroshium", "Fusarium oligoseptatum", "Fusarium obliquiseptatum",
    "Fusarium papillatum", "Fusarium pseudensiforme", "Fusarium protoensiforme",
    "Fusarium tuaranense", "Fusarium piperis", "Fusarium akasia", "Fusarium awan",
    "Fusarium mekan", "Fusarium rekanum", "Fusarium variasi", "Fusarium warna",
    "Fusarium duplospermum", "Fusarium illudens",
    "Fusarium virguliforme", "Fusarium phaseoli", "Fusarium tucumaniae",
    "Fusarium cuneirostrum", "Fusarium brasiliense", "Fusarium azukiicola",
    "Fusarium phalaenopsidis",  # (m)
    "Fusarium drepaniforme",    # (m)
}


def asm_key(path: str) -> str:
    name = path.rsplit("/", 1)[-1]
    for suf in (".fa.gz", ".fna.gz", ".fa"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def busco_index() -> dict[str, float]:
    out = {}
    for shard in BUSCO_DIR.iterdir():
        for f in shard.iterdir():
            if ".BUSCO_summary." in f.name:
                out[f.name.split(".BUSCO_summary.")[0]] = f
    return out


def busco_c(f: Path) -> float | None:
    m = re.search(r"C:([\d.]+)%", f.read_text())
    return float(m.group(1)) if m else None


def fsp_label(species_in: str) -> str:
    m = re.search(r"f\. ?sp\. ?([A-Za-z-]+)", species_in or "")
    return m.group(1) if m else ""


def load_dist(keep: set[str]) -> dict[tuple[str, str], float]:
    d = {}
    p = subprocess.Popen(["zstd", "-dc", str(SEL / "bfd_fusarium_all_vs_all.mash_dist.tsv.zst")],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        a, b, dist, _, _ = line.split("\t")
        a, b = asm_key(a), asm_key(b)
        if a in keep and b in keep:
            d[(a, b)] = float(dist)
    p.wait()
    return d


def single_linkage(ids: list[str], d, thr: float) -> list[list[str]]:
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if d.get((a, b), 1.0) < thr:
                parent[find(a)] = find(b)
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)
    return list(groups.values())


def main() -> int:
    inv = [r for r in csv.DictReader(open(SEL / "bfd_fusarium_inventory.tsv"), delimiter="\t")]
    log = []
    ok = {}
    for r in inv:
        if r["status"] != "ok":
            log.append((r["dirname"], r["asmid"], r["species"], "inventory", r["status"]))
            continue
        ok[r["asmid"]] = r
    bidx = busco_index()
    for a, r in ok.items():
        r["busco_genome_C"] = busco_c(bidx[a]) if a in bidx else None
        r["fsp"] = fsp_label(r["species_in"])
        r["n50_i"] = int(r["n50"] or 0)

    anchor_ids = set(ANCHORS.values())
    d = load_dist(set(ok))
    for a, r in ok.items():
        ds = {k: d.get((a, v), 1.0) for k, v in ANCHORS.items()}
        r["nearest"] = min(ds, key=ds.get)
        r["d_FOSC"], r["d_FSSC"] = ds["FOSC"], ds["FSSC"]

    members = {"FOXY": [], "FSSC": []}
    for a, r in ok.items():
        sp = r["species"]
        name_foxy = sp in FOXY_NAMES
        name_fssc = sp in FSSC_NAMES or r["dirname"].startswith("Neocosmospora_")
        if not (name_foxy or name_fssc):
            if r["nearest"] in ("FOSC", "FSSC") and ((r["nearest"] == "FOSC" and r["d_FOSC"] <= FOXY_MAX_D) or
                                                     (r["nearest"] == "FSSC" and r["d_FSSC"] <= FSSC_MAX_D)):
                log.append((r["dirname"], a, sp, "name_not_in_complex_list",
                            f"mash places it in {r['nearest']} (d_FOSC={r['d_FOSC']:.4f} d_FSSC={r['d_FSSC']:.4f}); excluded by name rule"))
            continue
        grp = "FOXY" if name_foxy else "FSSC"
        bc = r["busco_genome_C"]
        if bc is None or bc < MIN_BUSCO_C:
            log.append((r["dirname"], a, sp, "busco", f"{grp}: BUSCO genome C={bc}"))
            continue
        if grp == "FOXY" and not (r["nearest"] == "FOSC" and r["d_FOSC"] <= FOXY_MAX_D):
            log.append((r["dirname"], a, sp, "mash_placement",
                        f"FOXY name but nearest={r['nearest']} d_FOSC={r['d_FOSC']:.4f}"))
            continue
        if grp == "FSSC" and not (r["nearest"] == "FSSC" and r["d_FSSC"] <= FSSC_MAX_D):
            log.append((r["dirname"], a, sp, "mash_placement",
                        f"FSSC name but nearest={r['nearest']} d_FSSC={r['d_FSSC']:.4f}"))
            continue
        members[grp].append(a)

    kept = {}
    for grp, ids in members.items():
        ids.sort()
        reps = []
        for g in single_linkage(ids, d, DEDUP_D):
            g.sort(key=lambda a: (-ok[a]["n50_i"], -(ok[a]["busco_genome_C"] or 0), a))
            reps.append(g[0])
            for a in g[1:]:
                log.append((ok[a]["dirname"], a, ok[a]["species"], "dedup",
                            f"{grp}: mash<{DEDUP_D} to kept {ok[g[0]]['dirname']}"))
        kept[grp] = sorted(reps)

    # FOXY stratified subsample
    fx = kept["FOXY"]
    n = len(fx)
    m = np.zeros((n, n))
    for i, a in enumerate(fx):
        for j in range(i + 1, n):
            v = d.get((a, fx[j]), 1.0)
            m[i, j] = m[j, i] = v
    Z = linkage(squareform(m, checks=False), method="average")
    lab = fcluster(Z, t=STRATUM_H, criterion="distance")
    strata = defaultdict(list)
    for a, k in zip(fx, lab):
        strata[int(k)].append(a)
    target = min(FOXY_TARGET, n)
    alloc = {k: max(1, math.floor(target * len(v) / n)) for k, v in strata.items()}
    alloc = {k: min(alloc[k], len(strata[k])) for k in alloc}
    rem = target - sum(alloc.values())
    order = sorted(strata, key=lambda k: -((target * len(strata[k]) / n) % 1))
    while rem > 0:
        progressed = False
        for k in order:
            if rem == 0:
                break
            if alloc[k] < len(strata[k]):
                alloc[k] += 1
                rem -= 1
                progressed = True
        if not progressed:
            break
    rng = random.Random(SEED)
    chosen = []
    with open(SEL / "foxy_strata.tsv", "w") as fh:
        fh.write("stratum\tn_strains\tn_selected\tfsp_labels\n")
        for k in sorted(strata):
            pool = sorted(strata[k])
            rng.shuffle(pool)
            first, seen = [], set()
            for a in pool:
                f = (ok[a]["fsp"] or "(none)") if ok[a]["species"] == "Fusarium oxysporum" else ok[a]["species"]
                if f not in seen:
                    seen.add(f)
                    first.append(a)
            rest = [a for a in pool if a not in first]
            pick = (first + rest)[: alloc[k]]
            chosen.extend(pick)
            fh.write(f"{k}\t{len(pool)}\t{len(pick)}\t{','.join(sorted(seen))}\n")
    # every FOXY species name present after dedup gets at least one strain
    have = {ok[a]["species"] for a in chosen}
    for sp in sorted({ok[a]["species"] for a in fx} - have):
        best = sorted((a for a in fx if ok[a]["species"] == sp),
                      key=lambda a: (-ok[a]["n50_i"], a))[0]
        chosen.append(best)
        log.append((ok[best]["dirname"], best, sp, "species_floor",
                    "FOXY: added so every FOSC species name keeps >=1 strain"))
    chosen_set = set(chosen)
    for a in fx:
        if a not in chosen_set:
            log.append((ok[a]["dirname"], a, ok[a]["species"], "subsample", "FOXY: not drawn in stratified subsample"))
    kept["FOXY"] = sorted(chosen)

    with open(SEL / "selected_strains.tsv", "w", newline="") as fh:
        cols = ["complex", "dirname", "asmid", "locustag", "species", "species_in", "strain", "fsp",
                "ncbi_taxid", "bioproject", "n50", "total_len", "n_contigs", "busco_genome_C",
                "d_FOSC", "d_FSSC", "proteins", "gff3", "scaffolds"]
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for grp in ("FOXY", "FSSC"):
            for a in kept[grp]:
                r = ok[a]
                w.writerow([grp] + [r.get(c, "") if not isinstance(r.get(c), float) else f"{r[c]:.4f}"
                                    for c in cols[1:]])
    with open(SEL / "selection_log.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dirname", "asmid", "species", "step", "reason"])
        w.writerows(log)
    steps = Counter(x[3] for x in log)
    print(f"FOXY members after QC/placement={len(members['FOXY'])} after dedup={n} strata={len(strata)} "
          f"selected={len(kept['FOXY'])}", file=sys.stderr)
    print(f"FSSC members after QC/placement={len(members['FSSC'])} after dedup={len(kept['FSSC'])}", file=sys.stderr)
    print(f"log steps: {dict(steps)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
