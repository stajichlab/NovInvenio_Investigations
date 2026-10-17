#!/usr/bin/env python3
"""Summary numbers for the PF11001 note: per-strain hits split by family class,
species-level survey (deduplicated by species name), hits per family class."""
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AFT = NII / "studies/fungi/Afumigatus_pangenome/analysis/pf11001/tables"
CPT = NII / "studies/fungi/coccidioides_pangenome/analysis/pf11001/tables"
out = []


def per_strain(path, famcol, bincol, label):
    rows = list(csv.DictReader(gzip.open(path, "rt"), delimiter="\t"))
    by = defaultdict(Counter)
    for r in rows:
        b = r.get(bincol, "NA")
        by[(r["short"], r["species"], r["group"])]["core" if b in ("core", "soft_core") else "noncore"] += 1
    bins = Counter(r.get(bincol, "NA") for r in rows)
    out.append(f"{label}: hit proteins by family bin: {dict(bins)}")
    for sp in sorted({k[1] for k in by}):
        ks = [k for k in by if k[1] == sp]
        core = [by[k]["core"] for k in ks]
        nonc = [by[k]["noncore"] for k in ks]
        out.append(f"{label} {sp}: strains={len(ks)} core-family hits per strain median={median(core)} "
                   f"range={min(core)}-{max(core)}; non-core hits per strain median={median(nonc)} "
                   f"range={min(nonc)}-{max(nonc)}; strains with >=1 non-core hit={sum(x > 0 for x in nonc)}")


per_strain(AFT / "genes.afum_all_strains.tsv.gz", "family", "family_bin", "afum(full_v072)")
per_strain(CPT / "genes.cocci_all_strains.tsv.gz", "immitis_family", "immitis_family_bin", "cocci(v5 families)")

# survey dedup by species
sp = {}
with open(AFT / "survey_genomes.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sp.setdefault(r["species"], r)
by_tg = defaultdict(list)
for s, r in sp.items():
    by_tg[r["taxon_group"]].append((s, int(r["n_pf11001_proteins"])))
rows = []
for tg, lst in sorted(by_tg.items()):
    rows.append(dict(taxon_group_first_seen=tg, n_species=len(lst), n_with_hit=sum(n > 0 for _, n in lst),
                     counts=";".join(f"{s}={n}" for s, n in sorted(lst))))
with open(AFT / "survey_by_species_dedup.tsv", "w") as fo:
    fo.write("taxon_group_first_seen\tn_species\tn_with_hit\tcounts\n")
    for r in rows:
        fo.write("\t".join(str(v) for v in r.values()) + "\n")
out.append(f"survey: {len(sp)} unique species names, {sum(int(r['n_pf11001_proteins'])>0 for r in sp.values())} with >=1 hit")
for line in out:
    print(line)
