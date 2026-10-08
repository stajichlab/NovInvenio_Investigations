#!/usr/bin/env python3.12
"""Cluster the Bd genomes from the Mash triangle (average linkage) and write a Newick tree."""
import sys, os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, to_tree
from scipy.spatial.distance import squareform
H = os.path.dirname(os.path.abspath(__file__))
lines = open(f"{H}/bd_mash_triangle.tsv").read().splitlines()
n = int(lines[0]); names = []; D = np.zeros((n, n))
for i, l in enumerate(lines[1:]):
    p = l.split("\t"); b = os.path.basename(p[0]).removeprefix("Batrachochytrium_dendrobatidis_").removesuffix(".scaffolds.fa.gz")
    names.append(b)
    for j, v in enumerate(p[1:]):
        D[i, j] = D[j, i] = float(v)
Z = linkage(squareform(D), "average")
def nw(node, ids):
    if node.is_leaf(): return ids[node.id]
    return "(" + nw(node.left, ids) + "," + nw(node.right, ids) + ")"
open(f"{H}/bd_mash_upgma.nwk", "w").write(nw(to_tree(Z), names) + ";\n")
print("distance min/median/max", D[np.triu_indices(n, 1)].min(), np.median(D[np.triu_indices(n, 1)]), D.max())
for t in (0.0005, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008):
    c = fcluster(Z, t, "distance"); sizes = sorted(np.bincount(c)[1:], reverse=True)
    print(t, len(sizes), [int(x) for x in sizes[:10]])
t = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
c = fcluster(Z, t, "distance")
with open(f"{H}/bd_mash_clusters_t{t}.tsv", "w") as f:
    f.write("strain\tcluster\n")
    for nm, cl in zip(names, c): f.write(f"{nm}\t{cl}\n")
anch = ["JEL423", "JEL197", "CLFT024-02", "CLFT001", "CLFT071", "CLFT065", "CLFT144", "TST75", "RTP6", "KBO_317", "KRBOOR_317", "RC5.1", "SRS812", "MexMkt", "Hung_2014", "NBRC106979"]
for a in anch:
    if a in names: print(a, c[names.index(a)])
