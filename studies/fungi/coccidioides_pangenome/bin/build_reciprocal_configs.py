#!/usr/bin/env python3
"""Build the two reciprocal pangenome samplesheets for this study from config.csv:

  config_immitis_in_posadasii_out.csv   GROUP = IN for C. immitis, OUT for C. posadasii
  config_posadasii_in_immitis_out.csv   GROUP = IN for C. posadasii, OUT for C. immitis

Every other column and the row order are copied from config.csv unchanged; only
GROUP is set, by the Species column. Rows of any other species are dropped (there
are none in config.csv as of 2026-09-24). The reciprocal comparison is this study's
main pangenome comparison (notes/pangenome-method-investigations/
2026-09-24-rescued-rerun-and-outgroup-polarity.md).

Writes a provenance record for each output into this study's DATA_MANIFEST.yaml
(lib/provenance.py), with config.csv's sha256 as the source checksum.

Hardcoded paths (not derived from __file__) per CLAUDE.md "Study-specific vs.
shared scripts".
"""
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
STUDY_DIR = NII_ROOT / "studies" / "fungi" / "coccidioides_pangenome"
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import append_manifest, build_record  # noqa: E402

IMMITIS = "Coccidioides immitis"
POSADASII = "Coccidioides posadasii"
OUTPUTS = {
    "config_immitis_in_posadasii_out.csv": {IMMITIS: "IN", POSADASII: "OUT"},
    "config_posadasii_in_immitis_out.csv": {POSADASII: "IN", IMMITIS: "OUT"},
}


def reciprocal_rows(rows: list[dict], group_by_species: dict[str, str]) -> list[dict]:
    out = []
    for r in rows:
        g = group_by_species.get(r["Species"])
        if g is None:
            continue
        out.append({**r, "GROUP": g})
    return out


def main() -> int:
    source = STUDY_DIR / "config.csv"
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    with open(source, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)

    records = []
    for name, rule in OUTPUTS.items():
        out_rows = reciprocal_rows(rows, rule)
        path = STUDY_DIR / name
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)  # CRLF, same as config.csv
            w.writeheader()
            w.writerows(out_rows)
        n_in = sum(r["GROUP"] == "IN" for r in out_rows)
        records.append(build_record(
            source_url=f"(internal -- {source.relative_to(NII_ROOT)}, sha256 {source_sha})",
            source_release="derived from this study's config.csv",
            license="Internal / unpublished (not yet released outside this project)",
            local_path=path.relative_to(NII_ROOT),
            checksum=hashlib.sha256(path.read_bytes()).hexdigest(),
            derived_by=(f"studies/fungi/coccidioides_pangenome/bin/build_reciprocal_configs.py: "
                        f"GROUP set by Species ({', '.join(f'{k} -> {v}' for k, v in rule.items())})"),
        ))
        print(f"wrote {name}: {n_in} IN, {len(out_rows) - n_in} OUT", file=sys.stderr)

    append_manifest(records, STUDY_DIR / "DATA_MANIFEST.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
