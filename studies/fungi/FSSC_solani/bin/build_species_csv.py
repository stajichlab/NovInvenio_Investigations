#!/usr/bin/env python3
"""Write <study>/species.csv from <study>/selection/selected_strains.tsv.

Adapted from studies/fungi/fusarium_FOXY_vs_FSSC/bin/build_species_csv.py.
Differences: the study dir is an argument, and Group is IN/OUT directly (one
direction per study, so config.csv is the samplesheet; no reciprocal step).
Every row is local BFD funannotate output (proteins, scaffolds, gff3 in
predict_results/). TaxonGroup is left empty: pangenome.nf fills ingroup rows
with its own mash clade labels (FILL_TAXON_GROUP).
Short = <first 4 letters of epithet>_<strain>, non [A-Za-z0-9._-] -> '_'.

Usage: /usr/bin/python3.12 build_species_csv.py --study-dir <study>
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

FIELDS = ["Short", "Species", "Strain", "Group", "TaxonGroup", "Protein_Source", "Protein_Accession",
          "Taxon_ID", "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession"]


def short_for(r: dict) -> str:
    words = r["species"].split()
    epi = next((w for w in words[1:] if w not in ("cf.", "sp.", "aff.")), words[-1])
    strain = r["strain"].split(";")[0].strip() or "_".join(r["asmid"].split("_")[:2])
    return re.sub(r"[^A-Za-z0-9._-]+", "_", f"{epi[:4]}_{strain}").strip("_")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-dir", required=True)
    study = Path(ap.parse_args().study_dir)
    rows = list(csv.DictReader(open(study / "selection" / "selected_strains.tsv"), delimiter="\t"))
    out, seen = [], {}
    for r in rows:
        s = short_for(r)
        if s in seen:
            s = f"{s}_{r['locustag']}"
        if s in seen:
            sys.exit(f"ERROR: Short collision {s}")
        seen[s] = r["dirname"]
        out.append({
            "Short": s, "Species": r["species"], "Strain": r["strain"], "Group": r["GROUP"],
            "TaxonGroup": "", "Protein_Source": "local_faa", "Protein_Accession": r["proteins"],
            "Taxon_ID": r["ncbi_taxid"], "Genome_Source": "local_genome", "Genome_Accession": r["scaffolds"],
            "GFF3_Source": "local_gff3", "GFF3_Accession": r["gff3"],
        })
    with open(study / "species.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)
    print(f"wrote {study / 'species.csv'}: {len(out)} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
