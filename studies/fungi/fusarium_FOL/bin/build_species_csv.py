#!/usr/bin/env python3
"""Write species.csv for this study from selection/fol_selected.tsv
(copied from ../fusarium_FOXY_vs_FSSC/bin/build_species_csv.py).

Every row is local: Protein_Source=local_faa, Genome_Source=local_genome,
GFF3_Source=local_gff3, all pointing into the BFD funannotate
predict_results/ directory of that strain (one annotation method for all
strains). The genome is predict_results/<tag>.scaffolds.fa, the same FASTA
the GFF3 coordinates refer to.

Group = IN (FOL) or OUT (the 4 FOSC outgroup strains).
TaxonGroup is left EMPTY on purpose: pangenome.nf then fills ingroup rows
with its own Mash+PCoA+k-means clade label (FILL_TAXON_GROUP), so pair
classification has more than one clade to test. (Coccidioides used the
species name here, which gives one ingroup clade.)

Short = <first 4 letters of the epithet>_<strain>, non [A-Za-z0-9._-]
characters replaced by '_'; ASMID accession used when the strain is empty.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOL")
FIELDS = ["Short", "Species", "Strain", "Group", "TaxonGroup", "Protein_Source", "Protein_Accession",
          "Taxon_ID", "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession"]


def short_for(r: dict) -> str:
    words = r["species"].split()
    epi = next((w for w in words[1:] if w not in ("cf.", "sp.", "aff.")), words[-1])
    strain = r["strain"].split(";")[0].strip() or r["asmid"].split("_")[0] + "_" + r["asmid"].split("_")[1]
    return re.sub(r"[^A-Za-z0-9._-]+", "_", f"{epi[:4]}_{strain}").strip("_")


def main() -> int:
    rows = list(csv.DictReader(open(STUDY / "selection" / "fol_selected.tsv"), delimiter="\t"))
    out, seen = [], {}
    for r in rows:
        s = short_for(r)
        if s in seen:
            s = f"{s}_{r['locustag']}"
        if s in seen:
            sys.exit(f"ERROR: Short collision {s}")
        seen[s] = r["dirname"]
        out.append({
            "Short": s, "Species": r["species"], "Strain": r["strain"], "Group": r["group"],
            "TaxonGroup": "", "Protein_Source": "local_faa", "Protein_Accession": r["proteins"],
            "Taxon_ID": r["ncbi_taxid"], "Genome_Source": "local_genome", "Genome_Accession": r["scaffolds"],
            "GFF3_Source": "local_gff3", "GFF3_Accession": r["gff3"],
        })
    with open(STUDY / "species.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)
    print(f"wrote species.csv: {len(out)} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
