#!/usr/bin/env python3.12
"""Cap the cyanobacteria ingroup at N genomes per genus, keeping the highest-N50 assemblies.
Outgroups (GROUP=OUT) are all kept. Genus = first word of Species.
Usage: select_capped.py --cap 8  -> config_capped<N>.csv + analysis/capped<N>_selection.tsv
Genome DNA is read from data_dir/dna (plain FASTA)."""
import argparse, csv, collections
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("--cap", type=int, default=8); a = ap.parse_args()
H = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/bacteria/cyanobacteria")
rows = list(csv.DictReader(open(H / "config.csv"))); fields = list(rows[0])
def n50(p):
    L, cur = [], 0
    for l in open(p):
        if l[0] == ">":
            if cur: L.append(cur)
            cur = 0
        else: cur += len(l.strip())
    if cur: L.append(cur)
    L.sort(reverse=True); tot = sum(L); acc = 0
    for x in L:
        acc += x
        if acc >= tot / 2: return len(L), x, tot
for r in rows:
    if r["GROUP"] == "IN": r["_q"] = n50(H / "data_dir" / "dna" / r["DNA"])
by = collections.defaultdict(list)
for r in rows:
    if r["GROUP"] == "IN": by[r["Species"].split()[0]].append(r)
keep, sel = [], []
for g, v in sorted(by.items()):
    v.sort(key=lambda r: (-r["_q"][1], r["_q"][0], r["Short"]))
    for i, r in enumerate(v):
        sel.append((g, r["Short"], *r["_q"], i < a.cap)); 
        if i < a.cap: keep.append(r["Short"])
kept = set(keep)
with open(H / f"config_capped{a.cap}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
    for r in rows:
        if r["GROUP"] == "OUT" or r["Short"] in kept: w.writerow({k: r[k] for k in fields})
with open(H / "analysis" / f"capped{a.cap}_selection.tsv", "w") as f:
    f.write("genus\tShort\tn_contigs\tN50\ttotal_bp\tkept\n")
    for s in sel: f.write("\t".join(map(str, s)) + "\n")
print(len(kept), "IN kept of", sum(len(v) for v in by.values()), "in", len(by), "genera")
