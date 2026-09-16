# Running log: Coccidioides pangenome onboarding

This is a step-by-step log of what was actually run to onboard this study,
kept up to date as work proceeds — not a design document (see `README.md` for
the study's goals, and
`notes/superpowers/specs/2026-09-15-coccidioides-pangenome-local-input.md` +
`notes/superpowers/plans/2026-09-15-coccidioides-pangenome-onboarding.md` in
the repo root for the full spec/plan this executes). This file exists so the
same local-annotation-freeze onboarding pattern can be reproduced for a future,
different local pangenome dataset without re-deriving every decision from
scratch.

## Why this study needed its own onboarding path

Every prior NII study's `species.csv` pointed at a UniProt proteome ID or an
NCBI Datasets accession — `bin/build_study_config.py` fetches those. This
study's input is 535 strains of **already-assembled, already-annotated
funannotate output** sitting on shared storage
(`/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation_freeze/20260112/{pep,DNA,GFF,CDS}/`)
— nothing to fetch, only files to point at. `build_study_config.py` already
supports this (`Protein_Source=local_faa` / `Genome_Source=local_genome` /
`GFF3_Source=local_gff3` — "path to an existing file → copy + provenance
record, no fetch"), so the only missing piece was a `species.csv` *generator*
for this dataset's specific, fully-deterministic filename convention
(`Coccidioides_{species}_{strain}.{proteins.fa,scaffolds.fa,gff3}`).

**Reuse guide for a future local pangenome dataset**: if a new dataset has the
same shape — one directory tree of pre-annotated protein/genome/GFF3 files,
one deterministic filename convention encoding species+strain — copy
`bin/audit_coccidioides_inputs.py` and `bin/build_coccidioides_species_csv.py`
into the new study's own `bin/`, and change exactly:
1. `FREEZE_ROOT` (and `SAMPLES_CSV`/`ASM_STATS_TSV` if a quality-audit source
   exists for the new dataset — if not, the audit script can be trimmed down
   to just `list_freeze_strains()` and skip the BUSCO-exclusion logic).
2. `_FREEZE_NAME_RE` (the filename regex) and the corresponding path
   reconstruction in `build_rows()`.
3. `TAXON_ID_BY_SPECIES` (verify any new taxon IDs the same way this study
   did — see below — not from memory).

Nothing else needs to change; both scripts are otherwise dataset-agnostic
Python using only the stdlib. This is a **copy-and-adapt** pattern, not a
generic CLI tool — see `CLAUDE.md`'s "study-specific vs. shared scripts" rule
for why (promote to shared `bin/` only once a *third* study needs this exact
shape).

## What was actually run, in order

### 1. Data-quality audit (`bin/audit_coccidioides_inputs.py`)

Reconciles three independent sources before generating `species.csv`:
- `annotation_freeze/20260112/pep/*.proteins.fa` filenames → the ground-truth
  strain list (535 strains, parsed via `Coccidioides_(immitis|posadasii)_(.+)\.proteins\.fa`).
- `Assembly/samples.csv` (559 rows, one per sequencing run) — cross-checked to
  catch strains that were sequenced but never made it through
  assembly/annotation.
- `Assembly/asm_stats.tsv` (BUSCO/assembly QC, keyed by `SampleID`) — **the
  join key is `<strain>.AAFTF`, not the bare strain name** (`.AAFTF` is a
  genome-processing job-label suffix baked into `SampleID`, not a distinct
  sample — confirmed against the freeze naming, and the user separately added
  a README note directly in the `Assembly/BUSCO/` results directory
  documenting this same `<strainID>.AAFTF` convention).

```bash
python3 studies/fungi/coccidioides_pangenome/bin/audit_coccidioides_inputs.py
```

**Real output** (2026-09-15):
- 535 strains in `annotation_freeze`, 559 in `samples.csv`.
- 24 `samples.csv` strains have no matching annotation output at all (dropped
  somewhere between sequencing and annotation — not investigated further,
  out of scope for this onboarding).
- 50 freeze strains have no `asm_stats.tsv` QC row (no BUSCO data available)
  — **included anyway, with a warning**, not excluded: absence of a QC record
  is not evidence of bad quality (this repo's design decision, not a default
  assumption).
- 8 strains have `BUSCO_Complete < 90%` (threshold decided from the actual
  distribution — see the spec, "Open Question 2" — not a round-number guess):
  `GT-153, NM_459, NM_9861, SD7, SD8, SJV_3, SOIL_582-1_S_NEW_CPA0065,
  VFC054_7SA`.

**Data-munging surprise worth flagging for reuse**: 2 of those 8 low-BUSCO
strains (`NM_9861` at 0.2% complete, `SD8` at 44.8%) have **no annotation_freeze
output at all** — their assemblies failed so badly they were never carried
through to annotation. So the BUSCO-exclusion set and the annotation-freeze
set aren't nested the way you'd first assume; a naive `535 - 8 = 527` estimate
is wrong. The actual generator output is `535 - 6 = 529` (the other 6 excluded
strains *are* in the freeze and get correctly dropped). **When reusing this
pattern on a new dataset: don't assume the QC-exclusion count and the
annotated-strain count are independent — always compute the actual
intersection, don't do exclusion-count arithmetic by hand.**

### 2. `species.csv` generation (`bin/build_coccidioides_species_csv.py`)

```bash
python3 studies/fungi/coccidioides_pangenome/bin/build_coccidioides_species_csv.py
```

Every row: `Protein_Source=local_faa`, `Genome_Source=local_genome`,
`GFF3_Source=local_gff3`, `Group=IN` (this study has **no ingroup/outgroup
design** — see below), `Taxon_ID` from a 2-entry table
(*C. immitis*=5501, *C. posadasii*=199306, genus *Coccidioides*=5500).

**Taxon ID provenance**: supplied directly by the user, then independently
verified in-session via:
```bash
module load taxonkit
echo -e "5501\n199306\n5500" | taxonkit lineage
```
confirming all three resolve to the expected
`...Onygenaceae;Coccidioides;Coccidioides {immitis,posadasii}` lineages against
the installed NCBI taxdump on this cluster. **When reusing this pattern: don't
hardcode a taxon ID from memory — either verify it the same way, or query it
programmatically (NCBI Datasets/E-utilities), per this repo's "no assumptions"
rule.**

**Result**: 529 rows, 169 *C. immitis* / 360 *C. posadasii*.

### 3. `config.csv` + `data_dir` materialization (existing shared script, unmodified)

```bash
python3 bin/build_study_config.py --study-dir studies/fungi/coccidioides_pangenome
```

This is the existing, already-shared `bin/build_study_config.py` — no changes
needed for this study; its `local_faa`/`local_genome`/`local_gff3` dispatch
already did exactly what was needed. Real disk cost: **~20 GB copied** into
`studies/fungi/coccidioides_pangenome/data_dir/{pep,dna,gff3}/` (529 files
each — confirmed via `du -sh`/`ls | wc -l` after running, matching a
pre-computed estimate of 20.7 GB from summing the source freeze directory's
real file sizes with `stat -L` before running). Copy, not symlink — a
deliberate choice (see the spec) given `/bigdata` had 11 TB free at the time.
`data_dir/` itself is gitignored; `config.csv` and `DATA_MANIFEST.yaml` are
normal tracked files (same as every other NII study — **an earlier draft of
the plan wrongly assumed these were gitignored too; they are not, only
`data_dir/` is**).

### 4. Per-species `config.csv` filtering (`bin/filter_config_by_taxon.py`)

The three planned runs (whole-set / immitis-only / posadasii-only) are
**not** implemented via `pangenome.nf`'s `GROUP` column — that column is the
pipeline's own ingroup/outgroup split for an internal Mash/clade-sketch step,
unrelated to species selection, and every row in this study is `Group=IN`
(no outgroup design here). Instead:

```bash
python3 studies/fungi/coccidioides_pangenome/bin/filter_config_by_taxon.py
```

filters the one whole-set `config.csv` by its `TaxonGroup` column into
`config_immitis.csv` (169 rows) and `config_posadasii.csv` (360 rows), both
pointing at the **same** `data_dir` — `pangenome.nf` resolves FASTA/GFF3 by
basename against `--pangenome_data_dir`, so no additional copying is needed
for the two sub-runs.

### 5. Run script (`run_pangenome.sh`)

`bin/run_study.sh` (the repo's normal study-runner) can't drive this — it's
wired to `main.nf`'s `--config`/`--data_dir` param names, not `pangenome.nf`'s
`--pangenome_samplesheet`/`--pangenome_data_dir`, and has no way to select a
different entry script. This study has its own runner instead:

```bash
NII_PIPELINE_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio \
    studies/fungi/coccidioides_pangenome/run_pangenome.sh <wholeset|immitis|posadasii> [extra nextflow args]
```

Result directories are named `results/mmseqs_<variant>/` (not bare
`results/<variant>/`) so a future `diamond`-backend rerun (currently
hard-disabled in `pangenome.nf` pending validation) stays distinguishable, per
the study README's own request. The script logs the pipeline's exact git
commit hash on every invocation, since `nf_NovInvenio` is a live branch
(`pangenome-profiling-module`) another agent may still be committing to.

## Upstream pipeline fix — landed 2026-09-15

**The pipeline's GFF3 gene-position parser only handled `protein_id=`
(NCBI-style GFF3).** This study's funannotate GFF3s have **no `protein_id=`
attribute at all** — every CDS row uses `Parent=` instead (verified: 0/535
GFF3s carry `protein_id=`; all carry `Parent=`; for one spot-checked strain,
the `Parent=` ID set matches the protein FASTA header set exactly, 8655/8655).
Unfixed, this doesn't error — it silently produces empty gene-position tables
for every strain, which cascades to empty family-position data and collapses
every co-gain/co-loss candidate to `insufficient_data` downstream.

Fixed directly in `nf_NovInvenio` (branch `pangenome-profiling-module`,
commit `9b7e39e`): `bin/pangenome_build_gene_positions.py` now prefers
`protein_id=` when present, falls back to CDS `Parent=` (split on comma for
the rare multi-transcript-shared-CDS case) otherwise, and cross-checks every
resolved ID against that strain's actual protein FASTA headers (new required
`--protein_dir` arg, wired through `modules/pangenome/positions.nf` and
`workflows/pangenome_profile.nf` as `"<data_dir_abs>/pep"`) rather than
trusting attribute presence alone. Hard-errors below 50%
resolved-and-FASTA-matching (genuine dialect mismatch); warns above 2%
unresolved. 7 new tests added
(`nf_NovInvenio/tests/test_pangenome_build_gene_positions.py`); full pipeline
suite still passes (471 passed, 3 skipped, `pixi run pytest`). Also updated
`Afumigatus_pangenome/NEXTFLOW_MIGRATION_NOTES.md` to mark this resolved,
since that's where the gap was originally flagged.

## Small-subset validation — passed 2026-09-15

**Direct parser check (no Nextflow) against all 529 real GFF3+FASTA pairs**:
```bash
python3 bin/pangenome_build_gene_positions.py \
    --config studies/fungi/coccidioides_pangenome/config.csv \
    --gff3_dir studies/fungi/coccidioides_pangenome/data_dir/gff3 \
    --protein_dir studies/fungi/coccidioides_pangenome/data_dir/pep \
    --groups IN --output <output.tsv>
```
`build_gene_positions: 529 strains parsed, 0 skipped (no GFF3)` — **zero
warnings, zero hard errors**: every strain's `Parent=` fallback resolved and
matched its own protein FASTA cleanly. Row-count-vs-protein-count check
(grouping the output by `Short`): **0/529 strains with a mismatch** — the
pass criterion from the spec, met exactly.

**10-strain Nextflow smoke test** (5 *C. immitis* + 5 *C. posadasii*, run
directly on a compute node already inside a SLURM allocation — no separate
`srun` needed since the session itself was already running under the
`stajichlab` partition, not the login node): **57/57 processes completed, 0
failed.** `gene_positions.tsv` and `family_positions.tsv` both fully
populated (85,850 rows each, exactly matching — confirms
`pangenome_build_family_positions.py`'s downstream join really is
dialect-agnostic, as predicted from code inspection before this fix). Row
count for strain `1M0` (8655) matches the direct-check count exactly.
`pair_classification.tsv` correctly has 0 candidate pairs at this sample
size — expected (statistical minimums like `min_co_carrying`/`min_clades`
can't be met with only 10 strains), not a failure signal, per the spec.
The all-`IN`/no-`OUT` samplesheet (this study has no ingroup/outgroup
design) ran through `MASH_SKETCH_INGROUP`/`ASSIGN_CLADES` without incident,
resolving the one previously-unverified assumption from the spec.

Next: resource sizing for the full-scale mmseqs clustering step (~4.6M
total proteins across 529 strains, not the 10-strain smoke-test scale), a
study-specific `stajichlab_queue.config`, then the three full runs
(whole-set, immitis, posadasii).

## Commits so far

```
4b121cb coccidioides_pangenome: add local-input data-quality audit script
326f5e6 coccidioides_pangenome: generate species.csv from local annotation_freeze
fb8856d coccidioides_pangenome: materialize config.csv + provenance manifest
a7eeb45 coccidioides_pangenome: add per-species config.csv filter for sub-runs
c13301d coccidioides_pangenome: generate per-species config.csv variants
88888ad coccidioides_pangenome: add pangenome.nf run script (wholeset/immitis/posadasii)
```

## Next steps (once the upstream fix lands)

1. Direct parser validation against all 529 real GFF3+FASTA pairs (no
   Nextflow) — pass criterion is `gene_positions.tsv` row count per strain
   equal to that strain's protein count, not any downstream pipeline stage
   output.
2. A 5-10 strain Nextflow smoke test under an interactive SLURM allocation
   (pipeline-wiring check, not the correctness check).
3. Resource sizing for the full-scale mmseqs clustering step (~4.6M total
   proteins across 529 strains — real scale, not the smoke-test scale), via a
   study-specific `stajichlab_queue.config` routing the heaviest processes to
   the `stajichlab` SLURM partition.
4. Full runs: whole-set, then immitis, then posadasii.

See `notes/superpowers/plans/2026-09-15-coccidioides-pangenome-onboarding.md`
for the full task-by-task detail.
