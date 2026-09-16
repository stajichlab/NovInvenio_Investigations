#!/usr/bin/env python3
"""Filter a pangenome config.csv down to one species, by TaxonGroup, for the
per-species sub-runs -- reuses the SAME data_dir as the whole-set run
(pangenome.nf resolves FASTA/GFF3 by basename against --pangenome_data_dir,
so a smaller samplesheet needs no extra copying). See
notes/superpowers/specs/2026-09-15-coccidioides-pangenome-local-input.md,
"Correction" note under "New study-specific script" for why this exists
instead of a GROUP-column filter (GROUP is pangenome.nf's ingroup/outgroup
split, not a species selector).
"""
from __future__ import annotations

import csv
from pathlib import Path

STUDY_DIR = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations") / "studies" / "fungi" / "coccidioides_pangenome"


def filter_rows(rows: list[dict], taxon_group: str) -> list[dict]:
    return [r for r in rows if r["TaxonGroup"] == taxon_group]


def main() -> int:
    config_csv = STUDY_DIR / "config.csv"
    with config_csv.open() as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)

    for taxon_group, out_name in [
        ("Coccidioides immitis", "config_immitis.csv"),
        ("Coccidioides posadasii", "config_posadasii.csv"),
    ]:
        filtered = filter_rows(rows, taxon_group)
        out_path = STUDY_DIR / out_name
        with out_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(filtered)
        print(f"Wrote {out_path}: {len(filtered)} strains")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
