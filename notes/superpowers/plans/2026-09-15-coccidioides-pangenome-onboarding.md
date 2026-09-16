# Coccidioides Pangenome Onboarding Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Onboard the 535-strain local Coccidioides annotation freeze into a new NII
study, run `nf_NovInvenio --pipeline pangenome` on it (whole-set + per-species splits),
and leave behind a documented, reusable local-pangenome-dataset onboarding pattern.

**Architecture:** One study-specific `species.csv` generator (deterministic filename
parsing, no fetch) feeds the existing `bin/build_study_config.py` `local_*` dispatch to
produce one shared `config.csv`/`data_dir` (527 of 535 strains after the
BUSCO_Complete < 90% exclusion, `Group=IN`, no outgroup).
Two additional samplesheets (immitis-only, posadasii-only) are filtered copies of that
one `config.csv`, reusing the same `data_dir` — no extra data copying. A study-specific
run script invokes `pangenome.nf` (not the shared `run_study.sh`, which is wired to
`main.nf`'s different param names). The pipeline-side GFF3 parser fix is a **separate,
external, already-in-progress task** on `nf_NovInvenio`'s `pangenome-profiling-module`
branch, owned by another agent — this plan's Task 1 tracks it as a dependency-only
checkpoint, not something this plan's worker implements.

**Tech Stack:** Python 3 (stdlib `csv`/`argparse`/`pathlib`, no new deps), Bash,
Nextflow (`nf_NovInvenio`, external pipeline, invoked by local-checkout path),
`pixi run pytest` for the NII repo's own test conventions (check `pixi.toml` — if no
test runner is configured yet for `bin/` scripts, plain `python3 -m pytest` in-place is
fine; this repo has no established Python test suite pattern for study-specific
scripts, so keep tests lightweight and colocated, matching this repo's existing
scripts which are not currently under pytest).

**Spec:** `notes/superpowers/specs/2026-09-15-coccidioides-pangenome-local-input.md`
(read this in full before starting — it has the verified fact base, thresholds, and
rejected alternatives this plan is built from).

**Note on commit templates below:** each task's commit message includes attribution
lines reflecting *this* plan-writing session's active attribution instructions. If
this plan is executed later by a different session/agent, use whatever attribution
convention is active for that execution instead of copying these lines verbatim if
they've gone stale.

## Global Constraints

- No study-specific script derives paths from `__file__`/cwd; hardcode
  `NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")`
  (per root `CLAUDE.md` and `~/.claude/CLAUDE.md`'s BASH_SOURCE/SLURM warning, same
  reasoning applied to directory distance instead of a SLURM work dir).
- All new scripts for this study go in `studies/fungi/coccidioides_pangenome/bin/`,
  never shared `bin/` (this is the first study needing this directory shape; CLAUDE.md
  promotes a pattern to shared `bin/` only on a *third* study needing it).
- No tracked data file is committed without a provenance record — `local_faa`/
  `local_genome`/`local_gff3` sources go through `build_study_config.py`'s existing
  provenance path unchanged; do not hand-write manifest entries.
- Every large or generated intermediate stays under `studies/fungi/
  coccidioides_pangenome/{data_dir,results}/`, both gitignored (class 1/3 per root
  CLAUDE.md) — never committed regardless of size.
- Test small before running full: no task that invokes SLURM under the `stajichlab`
  partition happens before the small-subset validation tasks (Task 6) pass.
- `stajichlab` is a real SLURM **partition** (`sinfo -p stajichlab`: 30-day time
  limit, 4 nodes — confirmed 2026-09-15), not an `--account` flag; nothing in
  `nf_NovInvenio`'s `conf/ucr_hpcc_slurm.config` currently routes any process to it
  (only `short`/`epyc`/`preempt` are wired) — see Task 8.
- Do not edit anything under `/bigdata/stajichlab/jstajich/projects/NovInvenio`
  (the `nf_NovInvenio` checkout) as part of this plan — that repo's fix is being
  applied by another agent on branch `pangenome-profiling-module`; this plan only
  *depends* on it landing (Task 1 is a checkpoint, not an implementation task).

---

### Task 1: Confirm upstream GFF3 parser fix has landed (checkpoint, not implementation)

**Files:** none in this repo — this task reads `/bigdata/stajichlab/jstajich/projects/NovInvenio` (a separate checkout, owned by another agent's session) to confirm status only.

**Context:** `nf_NovInvenio`'s `bin/pangenome_build_gene_positions.py` currently parses
only `protein_id=` (NCBI-style GFF3). Coccidioides' funannotate GFF3s have no
`protein_id=` at all (0/535 confirmed) — every CDS row uses `Parent=` instead. Another
agent, working on branch `pangenome-profiling-module` in that checkout, has confirmed
this diagnosis independently and said they will apply the fix themselves (prefer
`protein_id=`, fall back to `Parent=` split on comma, cross-check resolved IDs against
the protein FASTA, hard-error below 50%/warn above 2% unresolved, add two new pytest
fixtures). They also raised whether `pangenome_build_family_positions.py`'s downstream
join needs adjustment — verified in this session: it does not, because that function
does a plain dict lookup keyed on whatever string is in `gene_positions.tsv`'s
`protein_id` column (no regex, no `protein_id`-shape assumption anywhere in
`pangenome_build_family_positions.py`), so as long as the upstream fix emits the
FASTA-header token (not the raw GFF3 attribute) regardless of dialect, downstream is
already dialect-agnostic.

- [ ] **Step 1: Check branch status before starting any task below that needs the fix**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git log --oneline -5 -- bin/pangenome_build_gene_positions.py
git diff --stat HEAD -- bin/pangenome_build_gene_positions.py bin/pangenome_build_family_positions.py tests/
```

Expected once landed: a new commit touching `bin/pangenome_build_gene_positions.py`
and a new/modified test file under `tests/` covering the funannotate-dialect and
mismatched-dialect fixtures described above.

- [ ] **Step 2: If not yet landed, do NOT proceed past Task 6 (small-subset
  validation) or Task 8 (full run)** — Tasks 2, 3, 4, 5, and 7 (audit script,
  species.csv generator, `config.csv`/`data_dir` materialization, per-species
  filter, run script) have no dependency on the pipeline fix and can proceed
  regardless. Re-check this task's Step 1 before starting Task 6. (As of
  2026-09-15, checked directly: the fix has **not** landed yet — the only commit
  touching this area, `3de5252`, is the original file creation, still only
  `_PROTEIN_ID_RE`, no `Parent=` fallback, no test file. Confirm status again
  before relying on this.)

- [ ] **Step 3: Once landed, run the new upstream pytest fixtures directly to confirm
  they pass** (do not just trust the commit message):

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run pytest tests/test_pangenome_build_gene_positions.py -v
```

Expected: PASS, including the new funannotate-dialect and mismatched-dialect cases.

- [ ] **Step 4: No commit needed in this repo for this task** — it's a status
  checkpoint only. Move to Task 2 (or Task 5 if Tasks 2-4 are already done and this
  was the last blocker).

---

### Task 2: Data-quality audit script

**Files:**
- Create: `studies/fungi/coccidioides_pangenome/bin/audit_coccidioides_inputs.py`
- Test: `studies/fungi/coccidioides_pangenome/bin/test_audit_coccidioides_inputs.py`

**Interfaces:**
- Produces: a `Reconciliation` dataclass / plain dict with keys
  `strains_in_freeze: set[str]`, `strains_in_samples_csv: set[str]`,
  `strains_missing_annotation: set[str]` (in samples.csv, not in freeze),
  `strains_missing_qc: set[str]` (in freeze, not in `asm_stats.tsv`),
  `strains_excluded_low_busco: set[str]` (BUSCO_Complete < 90.0),
  `busco_by_strain: dict[str, float]` — Task 3's generator script imports and calls
  this module's `audit(freeze_root, samples_csv, asm_stats_tsv) -> Reconciliation`
  function directly (no subprocess/JSON round-trip).

This is a **new, freestanding script**, not an extension of anything existing —
nothing in this repo currently reconciles `samples.csv` against `annotation_freeze`
or `asm_stats.tsv`.

- [ ] **Step 1: Write the failing test for strain-set parsing from the freeze
  directory**

```python
# studies/fungi/coccidioides_pangenome/bin/test_audit_coccidioides_inputs.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from audit_coccidioides_inputs import list_freeze_strains

def test_list_freeze_strains_parses_species_and_strain(tmp_path):
    pep_dir = tmp_path / "pep"
    pep_dir.mkdir()
    (pep_dir / "Coccidioides_immitis_1M0.proteins.fa").touch()
    (pep_dir / "Coccidioides_posadasii_B3224.proteins.fa").touch()
    result = list_freeze_strains(pep_dir)
    assert result == {
        "1M0": "Coccidioides immitis",
        "B3224": "Coccidioides posadasii",
    }
```

- [ ] **Step 2: Run it to verify it fails**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 -m pytest studies/fungi/coccidioides_pangenome/bin/test_audit_coccidioides_inputs.py -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'audit_coccidioides_inputs'`.

- [ ] **Step 3: Write the minimal implementation for strain-set parsing**

```python
# studies/fungi/coccidioides_pangenome/bin/audit_coccidioides_inputs.py
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


@dataclass
class Reconciliation:
    strains_in_freeze: set[str] = field(default_factory=set)
    strains_in_samples_csv: set[str] = field(default_factory=set)
    strains_missing_annotation: set[str] = field(default_factory=set)
    strains_missing_qc: set[str] = field(default_factory=set)
    strains_excluded_low_busco: set[str] = field(default_factory=set)
    busco_by_strain: dict[str, float] = field(default_factory=dict)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 -m pytest studies/fungi/coccidioides_pangenome/bin/test_audit_coccidioides_inputs.py -v
```

Expected: PASS.

- [ ] **Step 5: Write the failing test for the samples.csv reconciliation**

```python
def test_missing_annotation_strains_reported(tmp_path):
    from audit_coccidioides_inputs import missing_annotation_strains
    freeze_strains = {"1M0", "B3224"}
    samples_csv = tmp_path / "samples.csv"
    samples_csv.write_text(
        "RunAcc,Strain,BioSample,Center,Experiment,Project,Organism,FileBase,Notes,LocusTag\n"
        "r1,1M0,,,,,,,,\n"
        "r2,B3224,,,,,,,,\n"
        "r3,NOTANNOT1,,,,,,,,\n"
    )
    result = missing_annotation_strains(samples_csv, freeze_strains)
    assert result == {"NOTANNOT1"}
```

- [ ] **Step 6: Run it to verify it fails**, then implement:

```python
def missing_annotation_strains(samples_csv: Path, freeze_strains: set[str]) -> set[str]:
    """Strains in samples.csv with no matching annotation_freeze output."""
    with samples_csv.open() as fh:
        sample_strains = {row["Strain"] for row in csv.DictReader(fh)}
    return sample_strains - freeze_strains
```

- [ ] **Step 7: Run test to verify it passes.**

- [ ] **Step 8: Write the failing test for the BUSCO join (the `.AAFTF` suffix bug
  Fable review caught) and the exclusion threshold**

```python
def test_busco_join_handles_aaftf_suffix(tmp_path):
    from audit_coccidioides_inputs import load_busco_complete, low_busco_strains
    asm_stats = tmp_path / "asm_stats.tsv"
    asm_stats.write_text(
        "SampleID\tBUSCO_Complete\n"
        "1M0.AAFTF\t97.2\n"
        "NM_9861.AAFTF\t0.2\n"
    )
    busco = load_busco_complete(asm_stats)
    assert busco == {"1M0": 97.2, "NM_9861": 0.2}
    assert low_busco_strains(busco, min_complete=90.0) == {"NM_9861"}
```

- [ ] **Step 9: Run it to verify it fails**, then implement:

```python
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
            if val:
                out[strain] = float(val)
    return out


def low_busco_strains(busco_by_strain: dict[str, float], min_complete: float = BUSCO_MIN_COMPLETE) -> set[str]:
    return {s for s, v in busco_by_strain.items() if v < min_complete}
```

- [ ] **Step 10: Run test to verify it passes.**

- [ ] **Step 11: Write the top-level `audit()` function and a CLI entry point that
  prints a human-readable reconciliation report**

```python
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
```

- [ ] **Step 12: Run the full test file and the CLI against real data**

```bash
python3 -m pytest studies/fungi/coccidioides_pangenome/bin/test_audit_coccidioides_inputs.py -v
python3 studies/fungi/coccidioides_pangenome/bin/audit_coccidioides_inputs.py
```

Expected test: all PASS. Expected CLI output (values from the spec's verified fact
base — confirm they still match, data may have changed since 2026-09-15): 535 strains
in freeze, 559 in samples.csv, 24 missing annotation, ~30 missing QC row, 8 excluded
for low BUSCO.

- [ ] **Step 13: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add studies/fungi/coccidioides_pangenome/bin/audit_coccidioides_inputs.py \
        studies/fungi/coccidioides_pangenome/bin/test_audit_coccidioides_inputs.py
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: add local-input data-quality audit script

Reconciles annotation_freeze against samples.csv and asm_stats.tsv,
applying the BUSCO_Complete < 90% exclusion threshold decided in the
2026-09-15 onboarding spec.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 3: `species.csv` generator script

**Files:**
- Create: `studies/fungi/coccidioides_pangenome/bin/build_coccidioides_species_csv.py`
- Test: `studies/fungi/coccidioides_pangenome/bin/test_build_coccidioides_species_csv.py`

**Interfaces:**
- Consumes: `audit_coccidioides_inputs.audit()` (Task 2) for the exclusion set.
- Produces: `studies/fungi/coccidioides_pangenome/species.csv` with header
  `Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`
  — this is the exact input `bin/build_study_config.py` (Task 4) reads.

- [ ] **Step 1: Write the failing test for one row's field mapping**

```python
# studies/fungi/coccidioides_pangenome/bin/test_build_coccidioides_species_csv.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_coccidioides_species_csv import row_for_strain, TAXON_ID_BY_SPECIES

def test_row_for_strain_immitis():
    row = row_for_strain(
        strain="1M0",
        species="Coccidioides immitis",
        pep_path=Path("/x/pep/Coccidioides_immitis_1M0.proteins.fa"),
        dna_path=Path("/x/DNA/Coccidioides_immitis_1M0.scaffolds.fa"),
        gff_path=Path("/x/GFF/Coccidioides_immitis_1M0.gff3"),
    )
    assert row["Short"] == "1M0"
    assert row["Species"] == "Coccidioides immitis"
    assert row["Strain"] == "1M0"
    assert row["Group"] == "IN"
    assert row["TaxonGroup"] == "Coccidioides immitis"
    assert row["Taxon_ID"] == "5501"
    assert row["Protein_Source"] == "local_faa"
    assert row["Protein_Accession"] == "/x/pep/Coccidioides_immitis_1M0.proteins.fa"
    assert row["Genome_Source"] == "local_genome"
    assert row["GFF3_Source"] == "local_gff3"

def test_taxon_id_table_matches_verified_values():
    # 5501/199306 verified via `taxonkit lineage` against the installed NCBI
    # taxdump on 2026-09-15 -- see spec Open Question 8.
    assert TAXON_ID_BY_SPECIES == {
        "Coccidioides immitis": "5501",
        "Coccidioides posadasii": "199306",
    }
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Implement**

```python
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
filter_config_by_taxon.py, Task 6), not by this file or by GROUP.

Reusability: FREEZE_ROOT and the filename regex below are the two things a
future study with the same "pre-existing funannotate annotation_freeze
directory" shape would need to change -- copy this file, edit those two
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
    recon = audit()
    species_by_strain = list_freeze_strains(freeze_root / "pep")
    rows = []
    for strain, species in sorted(species_by_strain.items()):
        if strain in recon.strains_excluded_low_busco:
            continue
        if strain in recon.strains_missing_qc:
            print(f"WARNING: {strain} has no asm_stats.tsv QC row -- including anyway", file=sys.stderr)
        pep = freeze_root / "pep" / f"Coccidioides_{species.split()[1]}_{strain}.proteins.fa"
        dna = freeze_root / "DNA" / f"Coccidioides_{species.split()[1]}_{strain}.scaffolds.fa"
        gff = freeze_root / "GFF" / f"Coccidioides_{species.split()[1]}_{strain}.gff3"
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
```

- [ ] **Step 4: Run test to verify it passes.**

- [ ] **Step 5: Run the CLI against real data and sanity-check the output**

```bash
mkdir -p studies/fungi/coccidioides_pangenome
python3 studies/fungi/coccidioides_pangenome/bin/build_coccidioides_species_csv.py
wc -l studies/fungi/coccidioides_pangenome/species.csv     # expect 528 (527 rows + header; 535 - 8 excluded, assuming no missing-QC exclusions and no filename-parse misses)
cut -d, -f4 studies/fungi/coccidioides_pangenome/species.csv | sort -u   # expect only "Group" and "IN"
cut -d, -f5 studies/fungi/coccidioides_pangenome/species.csv | sort | uniq -c  # expect ~171 immitis / ~364 posadasii minus however many of the 8 BUSCO-excluded strains fall in each species (not yet determined -- read the actual counts here, don't assume an exact split)
```

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/coccidioides_pangenome/bin/build_coccidioides_species_csv.py \
        studies/fungi/coccidioides_pangenome/bin/test_build_coccidioides_species_csv.py \
        studies/fungi/coccidioides_pangenome/species.csv
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: generate species.csv from local annotation_freeze

527 of 535 freeze strains (8 excluded for BUSCO_Complete < 90%) dispatch to
build_study_config.py's local_faa/local_genome/local_gff3 sources, Group=IN
throughout (no ingroup/outgroup design for this study).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 4: Materialize `config.csv` + `data_dir` (real ~20.7 GB copy — confirm before running)

**Files:**
- No new files to write by hand — this task runs the existing
  `bin/build_study_config.py` against Task 3's `species.csv`.
- Modify (generated): `studies/fungi/coccidioides_pangenome/config.csv`,
  `studies/fungi/coccidioides_pangenome/data_dir/{pep,dna,gff3}/`,
  `studies/fungi/coccidioides_pangenome/DATA_MANIFEST.yaml`.

**This step performs a real, ~20.7 GB disk write to shared `/bigdata` storage** (user
confirmed acceptable at this size in the 2026-09-15 spec review — copy, not symlink).
Confirm free space before running, since this is not easily reversible mid-run.

- [ ] **Step 1: Check available space on the target filesystem**

```bash
df -h /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
```

Confirm at least ~25 GB free (20.7 GB data + headroom) before proceeding.

- [ ] **Step 2: Run the existing shared script (no changes needed to it)**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 bin/build_study_config.py --study-dir studies/fungi/coccidioides_pangenome
```

Expected: no `ERROR:` lines; ends with
`Wrote studies/fungi/coccidioides_pangenome/data_dir (pep/dna/gff3 -- pass as --data_dir to nf_NovInvenio)`.

- [ ] **Step 3: Verify row/file counts match `species.csv`**

```bash
wc -l studies/fungi/coccidioides_pangenome/config.csv
ls studies/fungi/coccidioides_pangenome/data_dir/pep | wc -l
ls studies/fungi/coccidioides_pangenome/data_dir/dna | wc -l
ls studies/fungi/coccidioides_pangenome/data_dir/gff3 | wc -l
du -sh studies/fungi/coccidioides_pangenome/data_dir
```

Expected: `config.csv` row count == `species.csv` row count; each `data_dir` subdir
file count matches; total size ≈ 20.7 GB (allow some variance — this is real copied
data, not a fixed number).

- [ ] **Step 4: Spot-check one GFF3 landed correctly (not silently truncated/empty)**

```bash
grep -c $'\tCDS\t' studies/fungi/coccidioides_pangenome/data_dir/gff3/1M0.gff3
```

Expected: a large positive count (thousands), not 0.

- [ ] **Step 5: Commit `config.csv` and `DATA_MANIFEST.yaml`** — correction from an
  earlier draft of this plan, which wrongly assumed these were gitignored. Verified
  directly against `.gitignore` (only `studies/*/*/data_dir` is excluded, no
  trailing slash specifically because some studies symlink a shared `data_dir`) and
  against an existing study (`pezizo_set1` has both `config.csv` and
  `DATA_MANIFEST.yaml` tracked in git). `data_dir/` itself stays gitignored — do not
  add it.

```bash
git status studies/fungi/coccidioides_pangenome/    # sanity-check: data_dir/ should show as ignored, config.csv/DATA_MANIFEST.yaml as untracked-but-trackable
git add studies/fungi/coccidioides_pangenome/config.csv studies/fungi/coccidioides_pangenome/DATA_MANIFEST.yaml
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: materialize config.csv + provenance manifest

535-strain local annotation_freeze, dispatched through build_study_config.py's
local_faa/local_genome/local_gff3 sources (no fetch). data_dir/ itself stays
gitignored (~20.7 GB copied data).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 5: `config.csv` filter script for the immitis/posadasii sub-runs

**Files:**
- Create: `studies/fungi/coccidioides_pangenome/bin/filter_config_by_taxon.py`
- Test: `studies/fungi/coccidioides_pangenome/bin/test_filter_config_by_taxon.py`

**Interfaces:**
- Consumes: `studies/fungi/coccidioides_pangenome/config.csv` (Task 4's output,
  columns `GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup` per
  `build_study_config.py`'s docstring).
- Produces: `config_immitis.csv` / `config_posadasii.csv` in the same directory, same
  header, filtered rows, pointing at the same `data_dir` (no new files copied).

- [ ] **Step 1: Write the failing test**

```python
# studies/fungi/coccidioides_pangenome/bin/test_filter_config_by_taxon.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from filter_config_by_taxon import filter_rows

def test_filter_rows_by_taxon_group():
    rows = [
        {"TaxonGroup": "Coccidioides immitis", "Short": "1M0"},
        {"TaxonGroup": "Coccidioides posadasii", "Short": "B3224"},
    ]
    result = filter_rows(rows, taxon_group="Coccidioides immitis")
    assert [r["Short"] for r in result] == ["1M0"]
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Implement**

```python
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
import sys
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
```

- [ ] **Step 4: Run test to verify it passes.**

- [ ] **Step 5: Run against real `config.csv` (requires Task 4 done first) and verify
  counts**

```bash
python3 studies/fungi/coccidioides_pangenome/bin/filter_config_by_taxon.py
wc -l studies/fungi/coccidioides_pangenome/config_immitis.csv studies/fungi/coccidioides_pangenome/config_posadasii.csv
```

Expected: row counts sum (minus 2 header lines) to the whole-set `config.csv` row
count minus 1 header line.

- [ ] **Step 6: Commit the script, and separately the generated filtered configs**
  (correction: `config_immitis.csv`/`config_posadasii.csv` are the same class as
  `config.csv` — small, derived, tracked, per Task 4 Step 5's correction, not
  gitignored. Committing only the script and leaving these two untracked would
  violate this repo's provenance rule — a tracked-workflow-adjacent file with no
  record of what produced it.)

```bash
git add studies/fungi/coccidioides_pangenome/bin/filter_config_by_taxon.py \
        studies/fungi/coccidioides_pangenome/bin/test_filter_config_by_taxon.py
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: add per-species config.csv filter for sub-runs

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
git add studies/fungi/coccidioides_pangenome/config_immitis.csv \
        studies/fungi/coccidioides_pangenome/config_posadasii.csv
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: generate per-species config.csv variants

Derived from config.csv via bin/filter_config_by_taxon.py (TaxonGroup filter),
same data_dir reused for all three run variants.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 6: Small-subset validation (test small before full run)

**Files:** no new files — this task exercises Tasks 1-5's outputs.

**Depends on:** Task 1 (pipeline fix landed) — do not start this task until Task 1's
Step 3 (pytest fixtures passing on `nf_NovInvenio`) is confirmed.

- [ ] **Step 1: Direct parser validation against ALL 535 real GFF3+FASTA pairs (no
  Nextflow) — the primary correctness check, per spec §4**

**Correction from an earlier draft of this plan**: `pangenome_build_gene_positions.py`
is not a per-strain CLI — verified directly, its real interface is one call against
the whole samplesheet: `--config <samplesheet> --gff3_dir <dir> --groups
<comma-list, default IN,OUT> --output <one TSV, all strains, with a Short column>`.
Also: `build_study_config.py` names `data_dir` outputs `<Short>.pep.fa`/
`<Short>.dna.fa`/`<Short>.gff3` (verified directly against its source, lines ~140/
161/175/210-234) — **not** `.proteins.fa`/`.scaffolds.fa` (those are the *source*
freeze filenames, already handled correctly in Tasks 2/3, which read from the freeze
directory, not `data_dir`).

**Note: `$SCRATCH` is only set inside a SLURM job**, not on the login node (per
`~/.claude/CLAUDE.md`'s HPCC guidance) — this check is pure Python over already-local
files, fast enough to run on the login node, so use the session's own scratchpad
directory (or any writable path under the study dir) for its output, not `$SCRATCH`:

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
CHECK_OUT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/.check/gene_positions_check.tsv"
mkdir -p "$(dirname "$CHECK_OUT")"
python3 bin/pangenome_build_gene_positions.py \
    --config /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/config.csv \
    --gff3_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir/gff3 \
    --groups IN \
    --output "$CHECK_OUT"
```

(Re-check `bin/pangenome_build_gene_positions.py --help` before running — this is the
interface as of 2026-09-15's WIP commit `3de5252`; the landed fix (Task 1) may add
flags, e.g. a `--protein_fasta_dir` for the FASTA cross-check, without changing this
shape.)

- [ ] **Step 2: Verify row counts match protein counts for every strain (group the
  one output TSV by its `Short` column, not per-strain files)**

```python
# quick check script, run inline or as a one-off
import csv
from collections import Counter
from pathlib import Path

pep_dir = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir/pep")
positions_tsv = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/.check/gene_positions_check.tsv")

n_positions_by_strain = Counter()
with positions_tsv.open() as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        n_positions_by_strain[row["Short"]] += 1

mismatches = []
for pep_file in sorted(pep_dir.glob("*.pep.fa")):
    strain = pep_file.stem  # "<Short>.pep" -> strip once more if needed; confirm actual stem shape against a real file first
    strain = strain[:-4] if strain.endswith(".pep") else strain
    n_proteins = sum(1 for line in open(pep_file) if line.startswith(">"))
    n_positions = n_positions_by_strain.get(strain, 0)
    if n_positions != n_proteins:
        mismatches.append((strain, n_positions, n_proteins))
print(f"{len(mismatches)} strains with position/protein count mismatch")
for m in mismatches[:20]:
    print(m)
```

Expected: 0 mismatches (or a small number explainable by real multi-isoform genes, not
a systematic dialect failure) — this is the pass criterion, per spec §4, not any
downstream Nextflow stage output.

- [ ] **Step 3: 5-10 strain Nextflow smoke test (pipeline-wiring check, not the
  primary correctness check)** — pick 5 immitis + 5 posadasii strains spanning a range
  of BUSCO completeness. Build the tiny samplesheet with a one-off Python filter (11
  lines of `config.csv`, same header, by `Short` — do not hand-type CSV rows, copy
  them programmatically to avoid a transcription error). Run under an **interactive
  SLURM allocation**, not directly on the login node — mmseqs2 and Nextflow's own
  polling are real compute, not appropriate for a shared login node even at this
  small scale — and write outputs to the session scratchpad directory, not `/tmp`
  (per this session's own environment rules):

```bash
srun --partition=short --time=1:00:00 --cpus-per-task=4 --mem=8G --pty bash
# inside the allocation:
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
nextflow run pangenome.nf \
    --pangenome_samplesheet /path/to/config_smoketest.csv \
    --pangenome_data_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir \
    --pangenome_project coccidioides_smoketest \
    --outdir "$SCRATCH/coccidioides_smoketest_results" \
    -profile local
```

Expected: completes without error; `gene_positions.tsv` in the output has non-empty
rows per strain (spot-check, don't rely on `PAIR_CLASSIFICATION` output as pass/fail —
per spec §4, it will legitimately show `insufficient_data` at this sample size
regardless of correctness).

- [ ] **Step 4: Confirm the all-`IN`/no-`OUT` samplesheet doesn't hit any hidden
  assumption** (spec flagged this as unverified, not assumed safe) — check the smoke
  test's Nextflow trace/log for any warning/error related to an empty outgroup
  channel. If none appears and the run completes, this is resolved; if something
  breaks, that's a new, real finding to bring back before the full run.

- [ ] **Step 5: No commit for this task** (temporary smoke-test outputs under
  `$SCRATCH`, not committed).

---

### Task 7: Study-specific run script + result-folder naming

**Files:**
- Create: `studies/fungi/coccidioides_pangenome/run_pangenome.sh`

**Context:** `bin/run_study.sh` cannot be reused as-is — it invokes `nextflow run
"$PIPELINE"` with `--config`/`--data_dir` (main.nf's param names) and no way to select
`pangenome.nf` as the entry script. This study needs its own run script, same pattern
as `Afumigatus_pangenome`'s existing `run_*.sh` scripts.

- [ ] **Step 1: Write the script**

```bash
#!/usr/bin/bash
# Run nf_NovInvenio's pangenome.nf against this study, for one of three
# variants (whole-set / immitis / posadasii). bin/run_study.sh isn't usable
# here -- it's wired to main.nf's --config/--data_dir params, not
# pangenome.nf's --pangenome_samplesheet/--pangenome_data_dir.
#
# Usage: studies/fungi/coccidioides_pangenome/run_pangenome.sh <wholeset|immitis|posadasii> [extra nextflow args...]
#
# NII_PIPELINE_DIR must point at a local nf_NovInvenio checkout (pangenome.nf
# isn't runnable via a bare `stajichlab/nf_NovInvenio` git-fetch the way
# main.nf is, per this study's spec -- confirm this is still true before
# relying on it, or set NII_PIPELINE_DIR explicitly either way):
#   NII_PIPELINE_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio \
#       studies/fungi/coccidioides_pangenome/run_pangenome.sh wholeset

set -euo pipefail

VARIANT="${1:?Usage: run_pangenome.sh <wholeset|immitis|posadasii> [extra nextflow args...]}"
shift || true

STUDY_DIR="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome"
PIPELINE_DIR="${NII_PIPELINE_DIR:-/bigdata/stajichlab/jstajich/projects/NovInvenio}"

case "$VARIANT" in
    wholeset)  SAMPLESHEET="$STUDY_DIR/config.csv" ;;
    immitis)   SAMPLESHEET="$STUDY_DIR/config_immitis.csv" ;;
    posadasii) SAMPLESHEET="$STUDY_DIR/config_posadasii.csv" ;;
    *) echo "ERROR: variant must be wholeset, immitis, or posadasii" >&2; exit 1 ;;
esac

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: $SAMPLESHEET not found -- run build_coccidioides_species_csv.py, build_study_config.py, and filter_config_by_taxon.py first" >&2
    exit 1
fi

# Named per README's ask: keep mmseqs results distinguishable from a future
# diamond rerun (diamond backend is currently hard-disabled in pangenome.nf).
OUTDIR="$STUDY_DIR/results/mmseqs_${VARIANT}"
mkdir -p "$OUTDIR"

LAUNCH_DIR="$STUDY_DIR/.nf_launch/${VARIANT}"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

# PIPELINE_DIR is a live branch another agent may still be committing to --
# record exactly which commit this run used.
echo "== pipeline commit: $(git -C "$PIPELINE_DIR" rev-parse HEAD) ==" >&2

nextflow run "$PIPELINE_DIR/pangenome.nf" \
    --pangenome_samplesheet "$SAMPLESHEET" \
    --pangenome_data_dir "$STUDY_DIR/data_dir" \
    --pangenome_project "coccidioides_${VARIANT}" \
    --outdir "$OUTDIR" \
    "$@"
```

- [ ] **Step 2: `chmod +x` and check only the argument-parsing/error paths — do NOT
  let this actually reach `nextflow run`** (if Task 5 already ran, `config_immitis.csv`
  exists, and an unguarded `run_pangenome.sh immitis` call would launch a real
  Nextflow run on the login node):

```bash
chmod +x studies/fungi/coccidioides_pangenome/run_pangenome.sh
studies/fungi/coccidioides_pangenome/run_pangenome.sh badvariant   # expect the usage error, exits before reaching nextflow
```

Do not exercise the `wholeset`/`immitis`/`posadasii` branches here by actually
running them — if Task 5 already produced `config_immitis.csv`/
`config_posadasii.csv` (or Task 4 produced `config.csv`), calling this script with a
real variant name will reach the `nextflow run` line and launch a real run on the
login node. Confirm the script builds the right `nextflow run` command by reading it,
not by executing it; Task 6 Step 3's actual smoke test (under an `srun` allocation)
is the real end-to-end check of this invocation shape.

- [ ] **Step 3: Commit**

```bash
git add studies/fungi/coccidioides_pangenome/run_pangenome.sh
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: add pangenome.nf run script (wholeset/immitis/posadasii)

bin/run_study.sh can't drive this -- it's wired to main.nf's --config/
--data_dir param names, not pangenome.nf's --pangenome_samplesheet/
--pangenome_data_dir. Result dirs are named mmseqs_<variant> per the
README's ask to keep a future diamond rerun distinguishable.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 8: Full run under `stajichlab` account (only after Task 6 passes)

**Files:** none new — invokes Task 7's script with SLURM/account overrides.

- [ ] **Step 1: Confirm Task 6's validation fully passed** (0 position/protein
  mismatches across all 535 strains, smoke test completed cleanly) before running
  anything under the `stajichlab` account. Do not skip this gate.

- [ ] **Step 2: `stajichlab` is a SLURM *partition*, not an account — verified
  directly (`sinfo -p stajichlab`: `up 30-00:00:0`, 4 nodes).** Nothing in
  `nf_NovInvenio`'s `conf/ucr_hpcc_slurm.config` currently routes any pangenome
  process there (only `short`/`epyc`/`preempt` are wired — see that file's
  `.*CLUSTER_TIER1`/`.*TBLASTN_PER_STRAIN`/`.*COOCCURRENCE`/`.*HMMSEARCH_CHUNK`
  blocks). Write a small override config for this study rather than passing a
  nonexistent `--account` flag:

```groovy
// studies/fungi/coccidioides_pangenome/stajichlab_queue.config
// Route the three heaviest pangenome processes to the stajichlab partition
// for this full run, per the study's "expedite with stajichlab" instruction --
// NOT the whole subworkflow (GENE_POSITIONS/PRESENCE_MATRIX/MASH_SKETCH/
// PREFIX_*/RESCUE_PASS/PAIR_CLASSIFICATION etc. stay on short/epyc, which is
// fine for their much smaller resource needs). Does not touch nf_NovInvenio
// itself -- passed via -c at invocation time.
//
// CLUSTER_TIER1 explicitly overrides `time` here: it inherits label 'high_cpu'
// (cpus=max_cpus, memory=32GB) from conf/ucr_hpcc_slurm.config, but that
// label's default time (2h first attempt, 24h on retry) would burn a wasted
// 2h attempt on the stajichlab partition's 30-day limit before escalating --
// there's no reason to accept that on a partition sized for long jobs.
// (No HMMSEARCH_CHUNK here -- that process doesn't exist in
// workflows/pangenome_profile.nf; it was a leftover from the novelty/loss
// workflow's own SLURM config, verified absent from this subworkflow.)
process {
    withName: '.*CLUSTER_TIER1' {
        queue = 'stajichlab'
        time  = '24.h'
    }
    withName: '.*TBLASTN_PER_STRAIN|.*COOCCURRENCE' {
        queue = 'stajichlab'
    }
}
```

**Correction (caught in Fable re-review): `CLUSTER_TIER1` is not unresourced.** It
carries `label 'high_cpu'` (`modules/pangenome/prefix_and_cluster.nf`), which
`conf/ucr_hpcc_slurm.config`'s `withLabel: 'high_cpu'` block sets to
`cpus = params.max_cpus`, `memory = '32.GB'` — plus the process-level default
`time = { task.attempt > 1 ? '24.h' : '2.h' }`, on top of the `-C
ryzen|broadwell|cascade` node constraint already set for `.*CLUSTER_TIER1`
specifically. **The real gap**: routing to `stajichlab` without an explicit `time`
override keeps the 2h-first-attempt cap, which defeats the point of using a
30-day-limit partition for a long clustering job at 4.6M proteins — the job would
just fail at 2h and burn a wasted retry before escalating. Add an explicit `time`
to the override config (below) instead of relying on the inherited retry ladder.

Also: `sinfo -p stajichlab` shows 5 nodes total, but 3 (`c[01-03]`) are `down*` and
are the `abu_dhabi` node class mmseqs2 already can't run on (same reason `.*CLUSTER_TIER1`
already carries `-C ryzen|broadwell|cascade` upstream) — only `r11`/`r12` (`ryzen`)
are both usable and architecture-compatible, so this reroute has 2 real nodes to
land on, not "4" as stated in an earlier draft of this section.

- [ ] **Step 3: Run whole-set first**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
NII_PIPELINE_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio \
    studies/fungi/coccidioides_pangenome/run_pangenome.sh wholeset \
    -profile slurm \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio/conf/ucr_hpcc_slurm.config \
    -c studies/fungi/coccidioides_pangenome/stajichlab_queue.config
```

- [ ] **Step 4: Once whole-set completes successfully, run immitis and posadasii**

```bash
NII_PIPELINE_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio \
    studies/fungi/coccidioides_pangenome/run_pangenome.sh immitis \
    -profile slurm \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio/conf/ucr_hpcc_slurm.config \
    -c studies/fungi/coccidioides_pangenome/stajichlab_queue.config

NII_PIPELINE_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio \
    studies/fungi/coccidioides_pangenome/run_pangenome.sh posadasii \
    -profile slurm \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio/conf/ucr_hpcc_slurm.config \
    -c studies/fungi/coccidioides_pangenome/stajichlab_queue.config
```

- [ ] **Step 4b: Commit `stajichlab_queue.config`** (small, tracked, study-specific —
  same class as the run script):

```bash
git add studies/fungi/coccidioides_pangenome/stajichlab_queue.config
git commit -m "$(cat <<'EOF'
coccidioides_pangenome: route full-run SLURM jobs to the stajichlab partition

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

- [ ] **Step 5: No commit for run outputs** (`results/mmseqs_*/` is gitignored,
  class-3 published-asset-or-nothing data per root CLAUDE.md — publish via the
  repo's existing `bin/publish_report_release.sh`/`bin/publish_alignment_release.sh`
  once figures/reports are built, which is out of scope for this plan per the spec's
  "Explicitly not in scope" section).

---

## Explicitly out of scope for this plan (per spec)

- Figures/reports (open/closed pangenome plots, co-occurrence stats).
- Starship/starfish cluster analysis.
- Diamond-backend validation.
- Any edits inside `/bigdata/stajichlab/jstajich/projects/NovInvenio` — owned by the
  other agent's session on `pangenome-profiling-module`.
