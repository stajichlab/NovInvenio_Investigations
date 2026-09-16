#!/usr/bin/env python3
"""Reconcile Coccidioides local-input sources before species.csv generation.

Cross-checks the annotation_freeze protein directory (the file-count ground
truth for "strains we can actually onboard") against samples.csv (sequencing
runs) and asm_stats.tsv (BUSCO/assembly QC), and applies the BUSCO_Complete
< 90% exclusion threshold decided in
notes/superpowers/specs/2026-09-15-coccidioides-pangenome-local-input.md
(Open Question 2) -- 8/497 strains fall below this threshold, all clearly
failed/fragmented assemblies (e.g. 0.2%-83.8% complete), with a clean gap
before the next-lowest strain at 90.3%.

This is a study-specific script (hardcoded freeze directory shape and
NII_ROOT), following studies/fungi/coccidioides_pangenome/bin/'s convention
per root CLAUDE.md's "study-specific vs shared scripts" rule -- see the spec
above for why this one dataset's naming convention is deterministic enough
to justify a script at all (vs. the new-study skill's normal per-species
judgment call).
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FREEZE_ROOT = Path(
    "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci"
    "/Assembly/annotation_freeze/20260112"
)
SAMPLES_CSV = FREEZE_ROOT.parent.parent / "samples.csv"
ASM_STATS_TSV = FREEZE_ROOT.parent.parent / "asm_stats.tsv"
BUSCO_MIN_COMPLETE = 90.0

_FREEZE_NAME_RE = re.compile(r"^Coccidioides_(immitis|posadasii)_(.+)\.proteins\.fa$")


def list_freeze_strains(pep_dir: Path) -> dict[str, str]:
    """{strain: "Coccidioides {species}"} parsed from *.proteins.fa filenames."""
    out: dict[str, str] = {}
    for p in sorted(pep_dir.glob("Coccidioides_*.proteins.fa")):
        m = _FREEZE_NAME_RE.match(p.name)
        if not m:
            continue
        species, strain = m.group(1), m.group(2)
        out[strain] = f"Coccidioides {species}"
    return out


def missing_annotation_strains(samples_csv: Path, freeze_strains: set[str]) -> set[str]:
    """Strains in samples.csv with no matching annotation_freeze output."""
    with samples_csv.open() as fh:
        sample_strains = {row["Strain"] for row in csv.DictReader(fh)}
    return sample_strains - freeze_strains


def load_busco_complete(asm_stats_tsv: Path) -> dict[str, float]:
    """{strain: BUSCO_Complete} -- strips the '.AAFTF' SampleID suffix (this
    is a genome-processing job-label suffix, not a distinct sample; see
    Open Question 7 in the spec)."""
    out: dict[str, float] = {}
    with asm_stats_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sample_id = row["SampleID"]
            strain = sample_id[: -len(".AAFTF")] if sample_id.endswith(".AAFTF") else sample_id
            val = row.get("BUSCO_Complete", "").strip()
            if val and val.upper() != "NA":
                out[strain] = float(val)
    return out


def low_busco_strains(busco_by_strain: dict[str, float], min_complete: float = BUSCO_MIN_COMPLETE) -> set[str]:
    return {s for s, v in busco_by_strain.items() if v < min_complete}


@dataclass
class Reconciliation:
    strains_in_freeze: set[str] = field(default_factory=set)
    strains_in_samples_csv: set[str] = field(default_factory=set)
    strains_missing_annotation: set[str] = field(default_factory=set)
    strains_missing_qc: set[str] = field(default_factory=set)
    strains_excluded_low_busco: set[str] = field(default_factory=set)
    busco_by_strain: dict[str, float] = field(default_factory=dict)


def audit(freeze_root: Path = FREEZE_ROOT, samples_csv: Path = SAMPLES_CSV,
          asm_stats_tsv: Path = ASM_STATS_TSV) -> Reconciliation:
    freeze_strains_map = list_freeze_strains(freeze_root / "pep")
    freeze_strains = set(freeze_strains_map)
    with samples_csv.open() as fh:
        sample_strains = {row["Strain"] for row in csv.DictReader(fh)}
    busco = load_busco_complete(asm_stats_tsv)
    return Reconciliation(
        strains_in_freeze=freeze_strains,
        strains_in_samples_csv=sample_strains,
        strains_missing_annotation=sample_strains - freeze_strains,
        strains_missing_qc=freeze_strains - set(busco),
        strains_excluded_low_busco=low_busco_strains(busco),
        busco_by_strain=busco,
    )


def main() -> int:
    r = audit()
    print(f"Strains in annotation_freeze: {len(r.strains_in_freeze)}")
    print(f"Strains in samples.csv: {len(r.strains_in_samples_csv)}")
    print(f"Missing annotation ({len(r.strains_missing_annotation)}): "
          f"{sorted(r.strains_missing_annotation)}")
    print(f"Missing QC/asm_stats row ({len(r.strains_missing_qc)}): "
          f"{sorted(r.strains_missing_qc)}")
    print(f"Excluded, BUSCO_Complete < {BUSCO_MIN_COMPLETE}% "
          f"({len(r.strains_excluded_low_busco)}): "
          f"{sorted(r.strains_excluded_low_busco)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
