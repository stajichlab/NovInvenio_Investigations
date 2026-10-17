#!/usr/bin/env python3
"""Build a manifest of proteomes to search with PF11001.

Columns: set, study, group, species, short, taxon_group, pep_path.
set = 'afum' (A. fumigatus pangenome config.csv), 'cocci' (coccidioides
config.csv, 530 rows), or 'survey' (other studies, deduplicated by pep
file name).
"""
import csv
import sys
from pathlib import Path

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
SURVEY = [
    "eukaryotes/zoosporic_holozoa_vs_dikarya",
    "fungi/agaricomycetes_pairwise",
    "fungi/chytrid_vs_mucoromycota",
    "fungi/malassezia_vs_ustilaginomycotina",
    "fungi/metschnikowiaceae_cah",
    "fungi/serinales_vs_pichiales_v2",
    "fungi/neocallimastigomycota_vs_early_fungi",
    "fungi/pezizo_set1",
    "fungi/sordariales_shallow",
    "fungi/yeast_filamentous",
    "fungi/mushrooms_tremella",
    "fungi/zoosporic_dikarya",
]


def rows(study, config="config.csv"):
    d = NII / "studies" / study
    with open(d / config, newline="") as fh:
        for r in csv.DictReader(fh):
            pep = d / "data_dir" / "pep" / r["Protein"]
            yield r, pep


def main(out):
    seen = set()
    with open(out, "w") as fo:
        fo.write("set\tstudy\tgroup\tspecies\tshort\ttaxon_group\tpep_path\n")
        for setname, study in (("afum", "fungi/Afumigatus_pangenome"),
                               ("cocci", "fungi/coccidioides_pangenome")):
            for r, pep in rows(study):
                if not pep.exists():
                    print(f"missing {pep}", file=sys.stderr)
                    continue
                fo.write(f"{setname}\t{study}\t{r['GROUP']}\t{r['Species']}\t{r['Short']}\t"
                         f"{r['TaxonGroup']}\t{pep}\n")
        for study in SURVEY:
            for r, pep in rows(study):
                if not pep.exists():
                    print(f"missing {pep}", file=sys.stderr)
                    continue
                if pep.name in seen:
                    continue
                seen.add(pep.name)
                fo.write(f"survey\t{study}\t{r['GROUP']}\t{r['Species']}\t{r['Short']}\t"
                         f"{r['TaxonGroup']}\t{pep}\n")


if __name__ == "__main__":
    main(sys.argv[1])
