#!/usr/bin/env python3
"""Summarise lichen_taxcheck.sh: best hit per protein, by lineage group, per proteome.

Input: analysis/lichen_check/hits.tsv.zst (diamond outfmt 6: qseqid sseqid pident length
evalue bitscore staxids skingdoms sphylums sscinames) and queries.txt (one query name per
line). Queries are named <Short>|<protein id>.
The best hit is the row with the highest bitscore for the query (first on ties).
Each best hit is assigned to a group by walking its NCBI taxon id up nodes.dmp (the taxonomy the
nr database was built with, NODES below): Fungi, Viridiplantae, Metazoa, other Eukaryota,
Bacteria, Archaea, Viruses, or other. If staxids lists several ids, the first one is used.
The diamond kingdom columns are not used: many bacterial taxa have no kingdom rank.
Output: analysis/lichen_check/summary.tsv, one row per (proteome, group), and per-protein
calls in best_hits.tsv.gz. Proteins with no hit are counted as "no_hit".
"""
import csv, gzip, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

D = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/yeast_convergence/analysis/lichen_check")
NODES = Path("/srv/projects/db/ncbi/preformatted/20260128/nodes.dmp")
GROUPS = {4751: "Fungi", 33090: "Viridiplantae", 33208: "Metazoa", 2759: "other_Eukaryota",
          2: "Bacteria", 2157: "Archaea", 10239: "Viruses"}   # order of search = nearest ancestor first


def load_parents(path=NODES):
    parent = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t|\t")
            parent[int(f[0])] = int(f[1])
    return parent


def group_of(taxid, parent, cache):
    """Nearest ancestor (or self) in GROUPS; 'other' if none (e.g. unclassified)."""
    if taxid in cache:
        return cache[taxid]
    t, seen = taxid, []
    result = "other"
    while t in parent:
        seen.append(t)
        if t in GROUPS:
            result = GROUPS[t]
            break
        if parent[t] == t:
            break
        t = parent[t]
    for s in seen:
        cache[s] = result
    cache[taxid] = result
    return result


def main():
    parent = load_parents()
    cache = {}
    best = {}
    p = subprocess.Popen(["zstd", "-dc", str(D / "hits.tsv.zst")], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        q, bits = f[0], float(f[5])
        if q not in best or bits > best[q][0]:
            best[q] = (bits, f[6], f[8], f[9])
    totals = Counter(l.split("|", 1)[0] for l in (D / "queries.txt").read_text().split())
    calls = defaultdict(Counter)
    with gzip.open(D / "best_hits.tsv.gz", "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["query", "bitscore", "staxid", "group", "phylum", "sciname"])
        for q, (b, tx, ph, sn) in best.items():
            first = tx.split(";")[0]
            g = group_of(int(first), parent, cache) if first.isdigit() else "other"
            sp = q.split("|", 1)[0]
            calls[sp][g] += 1
            w.writerow([q, b, first, g, ph, sn])
    with open(D / "summary.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["proteome", "group_of_best_hit", "n_proteins", "fraction_of_all_proteins", "n_proteins_total"])
        for sp, tot in totals.items():
            hit = sum(calls[sp].values())
            for k, n in sorted(calls[sp].items(), key=lambda x: -x[1]):
                w.writerow([sp, k, n, f"{n / tot:.4f}", tot])
            w.writerow([sp, "no_hit", tot - hit, f"{(tot - hit) / tot:.4f}", tot])


if __name__ == "__main__":
    main()
