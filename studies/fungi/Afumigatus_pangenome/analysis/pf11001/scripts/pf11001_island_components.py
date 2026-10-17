#!/usr/bin/env python3
"""Group PF11001 families that share an island (union-find on PF11001 families
only; two families are linked when one island contains both). Counts how many
separate groups ("independent" PF11001 island units) each run has, using only
the families counted by the pipeline test (pipeline_counted == True)."""
import csv
from collections import defaultdict
from pathlib import Path

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
T = {"full_v072": NII / "studies/fungi/Afumigatus_pangenome/analysis/pf11001/tables",
     "immitis_in_posadasii_out_v5": NII / "studies/fungi/coccidioides_pangenome/analysis/pf11001/tables",
     "posadasii_in_immitis_out_v2": NII / "studies/fungi/coccidioides_pangenome/analysis/pf11001/tables"}
out = open(T["full_v072"] / "pf11001_family_island_groups.tsv", "w")
out.write("run\tfamily_set\tn_families_counted\tn_counted_in_islands\tn_groups\tgroup_sizes\tgroups\n")
for run, d in T.items():
    for fset in ("pipeline_counted", "any_source"):
        fams = {}
        with open(d / f"pf11001_families.{run}.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if fset == "any_source" or r["pipeline_counted"] == "True":
                    fams[r["family"]] = r
        par = {f: f for f in fams}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        with open(d / f"pf11001_islands.{run}.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                fs = [f for f in r["pf11001_families"].split(",") if f in fams]
                for f in fs[1:]:
                    a, b = find(fs[0]), find(f)
                    if a != b:
                        par[a] = b
        in_isl = [f for f, r in fams.items() if int(r["n_islands"]) > 0]
        groups = defaultdict(list)
        for f in in_isl:
            groups[find(f)].append(f)
        sizes = sorted((len(g) for g in groups.values()), reverse=True)
        out.write(f"{run}\t{fset}\t{len(fams)}\t{len(in_isl)}\t{len(groups)}\t{','.join(map(str, sizes))}\t"
                  + " | ".join(",".join(sorted(g)) for g in groups.values() if len(g) > 1) + "\n")
out.close()
print(open(T["full_v072"] / "pf11001_family_island_groups.tsv").read())
