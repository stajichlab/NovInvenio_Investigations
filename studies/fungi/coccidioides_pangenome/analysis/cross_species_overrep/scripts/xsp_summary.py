#!/usr/bin/env python3
"""Long-format summary.tsv (section, item, metric, value, source) from summary_domains.tsv,
family_tests.tsv, family_partner_coverage.tsv and cross_genus_island_families.tsv."""
from __future__ import annotations

import csv
from pathlib import Path

OUT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/"
           "coccidioides_pangenome/analysis/cross_species_overrep")


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    out = []
    for r in rows(OUT / "summary_domains.tsv"):
        for k, v in r.items():
            if k != "definition":
                out.append(("domains", r["definition"], k, v, "summary_domains.tsv"))
    for r in rows(OUT / "family_tests.tsv"):
        for k, v in r.items():
            if k not in ("cocci_run", "test"):
                out.append(("families", f"{r['cocci_run']}|{r['test']}", k, v, "family_tests.tsv"))
    for r in rows(OUT / "family_partner_coverage.tsv"):
        for k in ("n_af_families", "n_with_best_hit", "frac_with_best_hit", "n_with_rbh", "frac_with_rbh"):
            out.append(("partner_coverage", f"{r['cocci_run']}|af_{r['af_set']}", k, r[k],
                        "family_partner_coverage.tsv"))
    xg = rows(OUT / "cross_genus_island_families.tsv")
    out.append(("cross_genus_island_families", "all", "n_families", len(xg),
                "cross_genus_island_families.tsv"))
    out.append(("cross_genus_island_families", "all", "n_all_three_shell_or_cloud",
                sum(r["af_bin"] in ("shell", "cloud") and r["v5_bin"] in ("shell", "cloud")
                    and r["v2_bin"] in ("shell", "cloud") for r in xg),
                "cross_genus_island_families.tsv"))
    out.append(("cross_genus_island_families", "all", "n_af_in_af293_starship_element",
                sum(r["af_starship_elements"] not in ("", "not_listed") for r in xg),
                "cross_genus_island_families.tsv"))
    out.append(("cross_genus_island_families", "all", "n_af_captain_yes",
                sum(r["af_captain"] == "yes" for r in xg), "cross_genus_island_families.tsv"))
    out.append(("cross_genus_island_families", "all", "n_af_without_pfam",
                sum(r["af_pfam"] == "" for r in xg), "cross_genus_island_families.tsv"))
    with open(OUT / "summary.tsv", "w") as fo:
        fo.write("section\titem\tmetric\tvalue\tsource\n")
        for r in out:
            fo.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
