#!/usr/bin/env python3.12
"""Write species.csv for Bd_pangenome from ~/projects/Bd/Pangenome/{aa,dna,gff3}.

Study-specific (hardcoded paths, strain naming rule). All 344 genomes are
Batrachochytrium dendrobatidis, funannotate output; no metadata table exists, so
Group is IN for every row (no outgroup) and TaxonGroup is blank until lineage
labels are supplied.
"""
import csv
import re
import sys
from pathlib import Path

BD = Path("/bigdata/stajichlab/jstajich/projects/Bd/Pangenome")
STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Bd_pangenome")
PREFIX = "Batrachochytrium_dendrobatidis_"

rows, seen = [], {}
for pep in sorted((BD / "aa").glob(f"{PREFIX}*.proteins.fa.gz")):
    strain = pep.name[len(PREFIX):-len(".proteins.fa.gz")]
    dna = BD / "dna" / f"{PREFIX}{strain}.scaffolds.fa.gz"
    gff = BD / "gff3" / f"{PREFIX}{strain}.gff3"
    for p in (dna, gff):
        if not p.exists():
            sys.exit(f"ERROR: missing {p}")
    clean = re.sub(r"[^A-Za-z0-9._-]", "_", strain)
    short = clean if clean.startswith("Bd") else f"Bd_{clean}"
    if short in seen:
        sys.exit(f"ERROR: Short {short!r} collides: {strain!r} and {seen[short]!r}")
    seen[short] = strain
    rows.append({
        "Short": short, "Species": "Batrachochytrium dendrobatidis", "Strain": strain,
        "Group": "IN", "TaxonGroup": "",
        "Protein_Source": "local_faa", "Protein_Accession": str(pep.resolve()), "Taxon_ID": "",
        "Genome_Source": "local_genome", "Genome_Accession": str(dna.resolve()),
        "GFF3_Source": "local_gff3", "GFF3_Accession": str(gff.resolve()),
    })

with open(STUDY / "species.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(f"wrote {len(rows)} rows")
