#!/usr/bin/env python3
"""Generate species.csv for the Coccidioides pangenome study from the local
funannotate annotation_freeze directory -- no fetch, every row is
Protein_Source=local_faa / Genome_Source=local_genome / GFF3_Source=local_gff3.

Every row gets Group=IN: this study has no ingroup/outgroup design (see
notes/superpowers/specs/2026-09-15-coccidioides-pangenome-local-input.md,
"New study-specific script" section) -- pangenome.nf's GROUP column is an
ingroup/outgroup split for its own internal Mash/clade-sketch step, not a
species-subset selector. The three planned runs (whole-set/immitis/posadasii)
are produced downstream by filtering config.csv (see
filter_config_by_taxon.py, Task 5), not by this file or by GROUP.

Reusability: FREEZE_ROOT and the filename regex (in
audit_coccidioides_inputs.py) are the two things a future study with the
same "pre-existing funannotate annotation_freeze directory" shape would need
to change -- copy this file and audit_coccidioides_inputs.py, edit those two
things, done. Not a generic CLI tool (only one study needs this shape so
far); promote to shared bin/ if a third study needs it (root CLAUDE.md's
own rule for when to do that).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from audit_coccidioides_inputs import (  # noqa: E402
    FREEZE_ROOT, NII_ROOT, audit, list_freeze_strains,
)

STUDY_DIR = NII_ROOT / "studies" / "fungi" / "coccidioides_pangenome"
TAXON_ID_BY_SPECIES = {
    "Coccidioides immitis": "5501",
    "Coccidioides posadasii": "199306",
}
FIELDNAMES = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession",
    "GFF3_Source", "GFF3_Accession",
]


def row_for_strain(strain: str, species: str, pep_path: Path, dna_path: Path, gff_path: Path) -> dict:
    return {
        "Short": strain,
        "Species": species,
        "Strain": strain,
        "Group": "IN",
        "TaxonGroup": species,
        "Protein_Source": "local_faa",
        "Protein_Accession": str(pep_path),
        "Taxon_ID": TAXON_ID_BY_SPECIES[species],
        "Genome_Source": "local_genome",
        "Genome_Accession": str(dna_path),
        "GFF3_Source": "local_gff3",
        "GFF3_Accession": str(gff_path),
    }


def build_rows(freeze_root: Path = FREEZE_ROOT) -> tuple[list[dict], set[str]]:
    """Returns (rows, excluded_strains). Excludes audit()'s low-BUSCO set;
    does NOT exclude strains missing a QC row entirely (included with a
    warning printed to stderr, per spec Open Question 2)."""
    recon = audit(freeze_root=freeze_root)
    species_by_strain = list_freeze_strains(freeze_root / "pep")
    rows = []
    for strain, species in sorted(species_by_strain.items()):
        if strain in recon.strains_excluded_low_busco:
            continue
        if strain in recon.strains_missing_qc:
            print(f"WARNING: {strain} has no asm_stats.tsv QC row -- including anyway", file=sys.stderr)
        species_word = species.split()[1]
        pep = freeze_root / "pep" / f"Coccidioides_{species_word}_{strain}.proteins.fa"
        dna = freeze_root / "DNA" / f"Coccidioides_{species_word}_{strain}.scaffolds.fa"
        gff = freeze_root / "GFF" / f"Coccidioides_{species_word}_{strain}.gff3"
        rows.append(row_for_strain(strain, species, pep, dna, gff))
    return rows, recon.strains_excluded_low_busco


def main() -> int:
    rows, excluded = build_rows()
    out_path = STUDY_DIR / "species.csv"
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out_path}: {len(rows)} strains ({len(excluded)} excluded for BUSCO_Complete < 90%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
