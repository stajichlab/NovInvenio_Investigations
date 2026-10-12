#!/usr/bin/env python3
"""Build the two reciprocal pangenome samplesheets from config.csv:

  config_FOXY_in_FSSC_out.csv   GROUP=IN for Group FOXY, OUT for FSSC
  config_FSSC_in_FOXY_out.csv   GROUP=IN for Group FSSC, OUT for FOXY

config.csv carries GROUP=FOXY/FSSC (species.csv Group, copied through by
bin/build_study_config.py). Every other column and the row order are kept.
Modeled on studies/fungi/coccidioides_pangenome/bin/build_reciprocal_configs.py.
Appends a provenance record per output to this study's DATA_MANIFEST.yaml.
"""
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
STUDY_DIR = NII_ROOT / "studies" / "fungi" / "fusarium_FOXY_vs_FSSC"
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import append_manifest, build_record  # noqa: E402

OUTPUTS = {
    "config_FOXY_in_FSSC_out.csv": {"FOXY": "IN", "FSSC": "OUT"},
    "config_FSSC_in_FOXY_out.csv": {"FSSC": "IN", "FOXY": "OUT"},
}


def main() -> int:
    src = STUDY_DIR / "config.csv"
    sha = hashlib.sha256(src.read_bytes()).hexdigest()
    with open(src, newline="") as fh:
        rd = csv.DictReader(fh)
        fields, rows = rd.fieldnames, list(rd)
    recs = []
    for name, rule in OUTPUTS.items():
        out = [{**r, "GROUP": rule[r["GROUP"]]} for r in rows if r["GROUP"] in rule]
        path = STUDY_DIR / name
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(out)
        n_in = sum(1 for r in out if r["GROUP"] == "IN")
        print(f"wrote {path.name}: {len(out)} rows ({n_in} IN, {len(out) - n_in} OUT)", file=sys.stderr)
        recs.append(build_record(
            source_url="(derived, no external source)", source_release="n/a", license="n/a",
            local_path=path, checksum=None,
            derived_by=f"bin/build_reciprocal_configs.py from config.csv (sha256 {sha})",
        ))
    append_manifest(recs, STUDY_DIR / "DATA_MANIFEST.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
