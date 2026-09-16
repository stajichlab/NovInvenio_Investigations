# Pangenome pipeline -> Nextflow migration notes

Status: **planning only, no migration work started.** This captures an independent
pipeline review (Fable model, 2026-09-15) done specifically to prep for moving this
study's `bin/` script chain into `nf_NovInvenio` as a reusable "pangenome profiling"
Nextflow DSL2 subworkflow — usable against any pangenome dataset NII points it at, not
just this one 295-strain *A. fumigatus* study. Per this repo's `CLAUDE.md`, pipeline
code belongs in `nf_NovInvenio`, not here — this file is the design record for that
future work, kept in NII next to the study it was derived from.

Read alongside `PANGENOME_CLUSTER_PROFILE_NOTES.md` (the study's own findings/results)
and `notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md` (the
general method design this was all built from).

## A. Process decomposition (DSL2)

Mostly a linear chain with two scatter points and a side branch. Proposed processes:

1. `PREFIX_PROTEOMES` (per-strain, scatter over config rows): Short-prefix headers
   (`><Short>|<id>`); optional `collapse_isoforms`. Currently an unscripted README
   one-liner (awk with hardcoded column indices) — needs a real script.
2. `CONCAT_PROTEOMES` (collect) -> `CLUSTER_TIER1` (whole-study; `cluster_backend.py`)
   -> emits `tier1_cluster.tsv`, `tier1_rep_seq.fasta`.
3. `PRESENCE_MATRIX` (`build_presence_matrix.py`).
4. Rescue branch (**redesigned 2026-09-15**, see the correctness-bug note in
   `PANGENOME_CLUSTER_PROFILE_NOTES.md` — the process shape below reflects the NEW
   per-strain design, not the original combined-DB one Fable reviewed):
   `EXTRACT_ABSENT_QUERIES` (whole-study; `extract_absent_family_queries.py`) ->
   per-strain scatter -> `PREFIX_GENOME` + `MAKEBLASTDB` (per strain, small) ->
   `TBLASTN_PER_STRAIN` (per strain, own DB + own absent-family query set) -> collect
   -> `RESCUE_PASS` (`rescue_pass.py` already accepts repeated `--tblastn_tsv`, so no
   separate merge process is needed).
5. `MASH_SKETCH` (whole-study, once) -> `DEREPLICATE` and `ASSIGN_CLADES` both
   consuming the same `.msh`/dist matrix (today each script re-sketches independently —
   real, fixable duplication).
6. `FREQUENCY_BINS` -> `COOCCURRENCE` (whole-study; see D for sharding notes).
7. Positions branch, independent of 4-6: `GENE_POSITIONS` (per-strain GFF3 parse, then
   collect) -> `FAMILY_POSITIONS` (needs cluster TSV).
8. `CAPTAIN_HMMSEARCH` (whole-study; hmmfetch + hmmsearch; currently no script, only
   `results/captain_gene/` outputs from ad hoc commands — needs a real script).
9. `PAIR_CLASSIFICATION` joins 6, 7, 8.
10. `synteny_windows.py` is a library, not a process (its CLI is a deliberate stub) —
    move to `lib/`, not a process on its own.

## B. What's migration-friendly vs. not

**Ready largely as-is** (clean argparse, no cwd writes): `build_presence_matrix.py`,
`rescue_pass.py`, `frequency_bins.py`, `build_gene_positions.py`,
`build_family_positions.py`, `cooccurrence.py`, `pair_classification.py`,
`split_fasta_chunks.py`, `extract_absent_family_queries.py` (new 2026-09-15).

**Not ready:**
- `cluster_backend.py`: identity/coverage thresholds hardcoded (0.9/0.8 tier-1,
  0.4/0.8 tier-2 mentioned in the design spec); writes `tmp_mmseqs/` to cwd; no
  `--threads`/`-p` plumbing to wire up to `task.cpus`.
- `dereplicate_strains.py` / `assign_clades.py`: assume `data_dir/dna/<file>` layout;
  map mash results by full path string (fragile if Nextflow staging changes path
  spelling); no mash `-p` threading flag exposed.
- Every script does `sys.path.insert(... parent.parent/"lib")` plus an
  `NOVINVENIO_ROOT` env-var lookup for `config_parser` — in `nf_NovInvenio` this
  collapses to the pipeline's own `lib/`, a real but mechanical rename/move.
- Manual steps not encoded in any script today: proteome/genome Short-prefix concat
  (currently a README awk one-liner with hardcoded `$7,$4` column indices); the
  DUF3435 captain-gene hmmsearch (ad hoc commands, no script); the `TaxonGroup` fill
  into `config.csv` (priority DAPC > Table S21 > Mash-clade fallback) — note
  `cooccurrence.py` reads clades from `config.csv` directly, not from
  `clade_assignments.tsv`, so this fill has to happen and be committed before
  co-occurrence can run.
- Launcher scripts: hardcoded `/rhome/.../NII` paths, `cd $NII_ROOT` before `pixi run`,
  hardcoded `--out` log paths — all fine for one study's SLURM scripts, all need to
  become Nextflow process directives / `params.*` instead.
- **Not documented anywhere reproducible**: the production rescue-pass run (before the
  2026-09-15 fix) combined 11 tblastn sources — a "safe-completed" partial output from
  a cancelled single job plus 10 chunked array outputs (see commit `bf15230`). That
  particular provenance is not reconstructable from the scripts alone; a Nextflow
  version needs one canonical execution path, not a hand-assembled patchwork.

## C. The tblastn scatter-gather pattern as a template

*(Fable's original review here covered the pre-2026-09-15 combined-genome-DB design,
which the rescue-pass correctness fix has now superseded with a per-strain design —
see `PANGENOME_CLUSTER_PROFILE_NOTES.md`'s "Rescue pass CORRECTNESS BUG" entry. The
findings below are kept for their general scatter-gather lessons, but the specific
`max_target_seqs` finding is now already fixed in the current per-strain script, not
still-open.)*

**Keep**: round-robin query/strain splitting for load balance (length- or
query-count-balanced, not just equal-count chunks); one shared resource built once,
not rebuilt per chunk; per-chunk zstd-compressed output; gather via repeated
`--tblastn_tsv` rather than a separate merge step.

**Do not copy into Nextflow:**
- `SLURM_ARRAY_TASK_ID` / `%02d` filename coupling — use `.splitFasta()` (or a Groovy
  channel operator for per-strain scatter) + `.collect()` instead of hand-rolled array
  indexing.
- A fixed `--array=N` count — derive chunk/task count from actual query or strain
  count, sized toward the HPCC skill's ~1-1.5h-per-task guidance rather than a
  hardcoded number (the study's launcher scripts do this by hand today: 10 -> 16 -> 24
  tasks, each a manual re-tuning).
- **[RESOLVED 2026-09-15, keep as a migration lesson]** `-max_target_seqs 5` against
  one combined multi-strain DB caps hits to 5 contigs across ALL strains per query, not
  per strain — a real correctness bug in the original design, not just a style
  concern. The per-strain-DB redesign (`run_rescue_pass_per_strain.sh`) fixes this
  structurally: a Nextflow version should scatter by STRAIN (one small DB per strain)
  from the start, not by query chunk against one shared DB, so this class of bug can't
  recur.
- `rescue_pass.py`: `fh.readlines()` loads each tblastn file fully into memory:; fine
  at current per-strain-chunk sizes, worth flagging if a future dataset's per-chunk
  output grows much larger. The `.zst` decompression path (`compressed_io.py`) never
  checks the `zstd` subprocess's exit status — a truncated/corrupted compressed file
  would currently fail silently with partial hits rather than erroring loudly.
- Rescue pass's `--min_pident 90` must track whatever identity threshold tier-1
  clustering itself used, or the two stages disagree about what counts as "the same
  gene" — should become one shared `params.*` value, not two independently-set flags.

## D. Resource/scaling concerns for a general dataset

Measured anchors from this study (293+2 strains, 2.79M proteins, 47,983 tier-1
families): mmseqs tier-1 clustering ~14 min; the original combined-DB tblastn (before
the per-strain redesign) needed days at 32 cores and was abandoned for a chunked
approach (~4h20m/chunk at 8 cores x 10 chunks); co-occurrence OOM-killed at an 8GB
memory cgroup, comfortable under 1h at 32GB (pre-exact-test-rewrite numbers); the
rescued matrix's larger eligible-family count (13k->26k) pushed candidate pairs ~3.9x
higher, needing a 64GB/48h allocation for the (now superseded) Monte Carlo version —
since 2026-09-15's exact-hypergeometric rewrite, the same computation is projected at
under 1h total, which should also simplify this allocation considerably; pair
classification ~4.5 min once its inputs exist.

**What should be parameterized as a function of dataset size, not hardcoded:**
- `COOCCURRENCE` memory/time should scale with eligible-family-count squared (the
  candidate-pair count driver) — worth exposing as a computed Nextflow resource
  directive rather than a fixed number per study. `--n_perms` is now dead code after
  the exact-test rewrite (kept only for CLI back-compat) and should be dropped
  entirely in a fresh Nextflow port rather than carried forward.
- tblastn (per-strain design): task count ~ strain count / target task duration;
  per-task memory ~ that batch's total genome size (small now, since each DB is
  single-strain).
- mmseqs/diamond cpus and memory ~ total protein count across the dataset.
- Mash sketching should run ONCE and feed both `DEREPLICATE` and `ASSIGN_CLADES`
  (today's scripts each re-sketch independently — real duplicated work to fix in the
  Nextflow version, not just a nice-to-have).

## E. Config surface for a general "run pangenome profiling on any species set" module

- `params.pangenome_samplesheet` — the `config.csv` equivalent (GROUP, Short, Protein,
  DNA, GFF3, TaxonGroup columns).
- `params.ingroup_label` / `params.outgroup_label` — currently hardcoded literal
  `"IN"`/`"OUT"` strings throughout the Python scripts.
- `params.cluster_backend` (`mmseqs`|`diamond`), `tier1_min_id`, `tier1_cov`,
  `tier2_min_id`, `id_sep`.
- `params.rescue.{enable, evalue, max_target_seqs, min_pident, min_qcov}` — note
  `min_pident` here must match whatever tier1_min_id clustering used (see C above).
- `params.mash_threshold`; `params.clades.{source: samplesheet|mash, k_range,
  k_fixed, n_pcoa}`.
- `params.freq_bins.{core, softcore, shell}` frequency cutoffs;
  `params.cooccurrence.{min_strain_count, fdr_alpha, screen_alpha}`.
- `params.captain_hmm` (Pfam accession/HMM path — DUF3435 is Starship/fungi-specific,
  must be optional/swappable per dataset, not hardcoded) and `params.pfam_hmm` path.
- `params.pair_class.{k, physical_threshold, trans_threshold, min_co_carrying,
  perm_alpha, min_clades}`.
- `params.collapse_isoforms` (bool).
- ~~The GFF3 `protein_id=` CDS attribute is assumed (NCBI style) by the
  gene-positions step — should add a hard error/warning...~~ **RESOLVED
  2026-09-15** (`nf_NovInvenio` `bin/pangenome_build_gene_positions.py`,
  branch `pangenome-profiling-module`): the gene-positions parser now prefers
  `protein_id=` when present, falls back to CDS `Parent=` (split on comma for
  the rare multi-transcript-shared-CDS case) when it isn't, and cross-checks
  every resolved ID against that strain's actual protein FASTA headers
  (`--protein_dir`, newly required) rather than trusting attribute presence
  alone. Hard-errors below 50% resolved-and-FASTA-matching (genuine dialect
  mismatch), warns above 2% unresolved. Found and fixed against this exact
  study's sibling: the Coccidioides pangenome study's funannotate GFF3s have
  0/535 `protein_id=` attributes, all `Parent=` instead — this fix is what
  unblocks that study's small-subset validation (its own
  `RUNNING_README.md`/onboarding spec in
  `studies/fungi/coccidioides_pangenome/` has the full diagnosis). Covered by
  7 new tests in `nf_NovInvenio`'s `tests/test_pangenome_build_gene_positions.py`.

## Open questions for whoever picks this up

- Does `nf_NovInvenio` already have a per-strain scatter pattern (e.g. from
  funannotate/genotyping stages) this can reuse directly, or does the scatter
  machinery need to be built fresh for this module?
- Where does the Mash sketch / dereplication / clade-assignment trio belong relative
  to other `nf_NovInvenio` stages that might already compute strain-level distance
  metrics (e.g. for phylogenomics) — is there duplication to avoid at the pipeline
  level, not just within this module?
- The TaxonGroup fill (DAPC > published-table match > Mash fallback) is currently a
  manual, study-specific research decision (which published tables to match against,
  with what priority) — does a general module need this to be a required
  `params.taxon_group_source`, or should it ship with a sensible generic default
  (Mash-only) and treat the richer priority scheme as this study's own customization?
