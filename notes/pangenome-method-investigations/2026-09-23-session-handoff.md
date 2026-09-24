# Session handoff — 2026-09-23

Companion to `2026-09-21-state-of-the-analysis.md` (the prior state of the analysis).
This note covers what changed since: issue #133 implemented, validated, and (as of
this writing) still being validated by a live rerun; a chain of infrastructure bugs
found and fixed along the way; and what's still open.

## 1. Code landed (all merged to `NovInvenio` main)

| PR | What |
|---|---|
| #152 | Issue #133: structural rescue filter — rejects a rescue hit that overlaps an existing gene of a different family, is under 150 aa, or falls in a ≥3-family repeat-hotspot window. This is the primary fix; see #133's issue body / `2026-09-20-rescue-pass-characterisation.md` for the evidence base. |
| #153 | Issue #131: `docs/pangenome-assumptions.md` — evidence-graded register of every `pangenome_*` param. Issue closed manually 2026-09-23 (PR didn't auto-link). |
| #154 | Issue #130: `ASSEMBLY_QUALITY_QC` process, runs unconditionally every pangenome run — Spearman/partial-Spearman of accessory-family-count vs assembly N50/contig-count, contig-terminus enrichment, auto-warning at \|rho\|>0.3. |
| #155 | Issue #134: `Diagnostic` framework (`diagnostics.tsv`/banners), `--pangenome_strict`. Only `rescue_redundancy` is wired to real data so far; `assembly_quality_confound`/`clade_structure_validity`/`gain_loss_skew` are honestly `not_computed` stubs — wiring #130's real output into this is unstarted follow-up work. |
| #156 | Vectorized `plot_presence_absence_matrix`/`accumulation_curve` (40x speedup, 2.7s→0.07s at 5k×500 benchmark). |
| #157 | Ungated core report figures (frequency dist, heatmap, accumulation curve, classification breakdown) from the `--pangenome_island_pfam_hmm`-only branch — they render on every run now. |
| #163 | Two report-viz bugs found via a live synthetic-data demo render: heatmap had no present/absent legend (fixed); domain-enrichment chart's `fdr_q=0` sentinel was a fixed 300 that could visually erase every other bar (now caps at 1.3× the max finite value). |

Data/infra fixes, merged to `NovInvenio_Investigations` main:
- Corrected two mislabeled strains (`485B-1_L_OLD_CPA0023` immitis→posadasii,
  `B3476` posadasii→immitis) across `species.csv` + 5 derived config CSVs, evidence
  in `CORRECTIONS.md` (mash-distance based, from `2026-09-20-phylogrouping-and-relatedness.md`).
- Rerouted `CLUSTER_TIER1`/`TBLASTN_PER_STRAIN`/`COOCCURRENCE` from `stajichlab`
  (stopped mid-run) to `preempt` (`stajichlab_queue.config`).
- **The actual root cause of 4 consecutive failed validation reruns**: not the
  queue move, not node instability, not a Nextflow status-polling bug (all
  suspected and ruled out in turn) — `EXTRACT_RESCUE_POSITIONS` has a **fixed
  `memory = '8.GB'`** in `NovInvenio/conf/ucr_hpcc_slurm.config` with no
  `task.attempt` scaling, so every retry OOM-kills identically (exit 137), and
  once retries exhaust, `errorStrategy` falls to `'finish'` — which produces
  "Execution cancelled -- Finishing pending tasks before exit," previously
  misread as an external kill. Fixed in `stajichlab_queue.config` (`32.GB *
  task.attempt`, same tier as `RESCUE_PASS`/`FREQUENCY_BINS`). Also: the user
  added `executor { queueGlobalStatus = true }` to the pipeline's own
  `ucr_hpcc_slurm.config` (landed directly, commit `24cf075`) for cluster-wide
  SLURM visibility — real fix, but not what actually resolved the repeat
  failures; the memory fix was.

## 2. Live run — the actual validation of everything above

**SLURM job 29026388**, submitted with `-p exfab -A exfab` (per instruction),
pipeline pinned to commit `3277bdb` (the #133 merge, NOT full main — deliberately,
to protect Nextflow's `-resume` cache; see below), `-resume`, running at
`studies/fungi/coccidioides_pangenome/.nf_launch/genus_vs_ureesii_rescue_structural/`.

As of this note: **~9.5h elapsed, still RUNNING, no failures since the memory fix
landed** (attempt 5; attempts 1-4 all failed, see the "Failed attempts" section
below if debugging a recurrence). Output so far confirms `COOCCURRENCE` and
`RESCUE_PASS` succeeded (`cooccurring_pairs.tsv.zst`, `presence_matrix.rescued.tsv`,
`pfam.domtblout` all present); no `report/`, `trans_modules/`, or
`pair_classification.tsv` yet — still working through `PAIR_CLASSIFICATION`/
islands/captain/report stages.

**Check status**: `sacct -j 29026388 --format=JobID,Partition,State,Elapsed`
**Output dir**: `studies/fungi/coccidioides_pangenome/results/rescue_structural_genus_vs_ureesii/output/pangenome/`
**Resubmit if needed** (same launch dir, `-resume` reuses everything already
completed): `sbatch /rhome/jstajich/.claude/jobs/a404179f/tmp/submit_nextflow_head.sh`
— this script lived in a job-scratch tmp dir. **Update (later 2026-09-23):
a durable copy is now at `studies/fungi/coccidioides_pangenome/.nf_launch/
genus_vs_ureesii_rescue_structural/submit_nextflow_head.sh`; resubmit with that
path.**

**The number that matters once it finishes**: the gain:loss ratio. Prior
(unrescued, broken) run: 99.3% gain / 0.7% loss (143:1). A. fumigatus, where
rescue worked correctly: 77.5%/10.9%/11.6% (7.1:1, ambiguous non-zero). If this
run's ratio moves from 143:1 toward ~7:1 with non-zero ambiguous, the whole
diagnosis is confirmed. If it doesn't, something else is still wrong — go back
to `2026-09-21-state-of-the-analysis.md`'s section 5 before assuming the fix
failed; check the rescue funnel stats first (now in `diagnostics.tsv` per #134).
**Correction (later 2026-09-23):** this run is pinned at `3277bdb`, which
predates #154 (`ASSEMBLY_QUALITY_QC`) and #155 (`DIAGNOSTICS`). It will NOT
write `diagnostics.tsv` or the `assembly_quality_*` tables. Get the rescue
funnel counts from the `RESCUE_PASS` task's `.command.err` in the launch dir's
`work/` instead, or run `pangenome_diagnostics.py` / the QC script by hand on
this run's outputs.

Funnel counts from that `.command.err` (task `40/e82441`, cached):
530 files, 131,025,237 hit rows parsed, 20,005,269 passed identity/coverage.
Structural filter rejected 15,604,028 (78.0%) as overlapping a different
family's gene, 3,180,585 (15.9%) as query rep < 150 aa, and 475,558 (2.4%) as
repeat-hotspot. 745,098 rows (3.7%) remained, which gave 572,880 unique
(family, strain) ABSENT→GENOME_ONLY cells. Under #134's rule (overlap/passed
> 0.5), `rescue_redundancy` would be TRIGGERED for this run. The filter did
remove those hits; the trigger says most rescuable hits were double-counts,
not that they reached the matrix. Whether 572,880 applied cells is a correct
number is not yet tested -- the gain:loss ratio is the test.

### Failed attempts 1-4, for context if this recurs
1. Job 28984843: failed when `stajichlab` was stopped mid-run, 439 tblastn tasks
   exhausted retries → prompted the `preempt` queue reroute.
2. Job 28994241 (`-resume`): reached `COOCCURRENCE` (58M/145M pairs) before
   failing — same signature, cause not yet found.
3. Job 29000889 (`-resume`, worktree recreated after original got deleted by
   another agent's cleanup): cache reuse regressed badly (`cached=4` vs
   attempt 2's `cached=1060`) — `git worktree add` resets mtimes, and
   Nextflow's default resume cache keys on path+size+mtime, not content
   checksum. Still failed, same signature, 21.7M/145M pairs this time.
4. Job 29019928 (`-resume`, with `queueGlobalStatus=true` picked up via `-c`
   from the main checkout — config files aren't part of the resume-cache hash,
   so this was safe): got to 72.5M/145M pairs, failed **again**, same
   signature. This is what led to actually reading `EXTRACT_RESCUE_POSITIONS`'s
   exit code (137, not just "cancelled") and finding the real bug.

**Lesson for next time**: don't delete/recreate a pipeline worktree mid-validation
if you care about `-resume` cache reuse. `-c` config files are always safe to
swap between attempts; the pipeline checkout itself is not, unless pinned to the
exact same commit AND never re-checked-out.

## 3. What's still open

- **#132 remaining sensitivity sweeps** (core/shell cutoffs, `pair_class_k`,
  Leiden resolution) — only priority 1 (tier-1 clustering) has been swept, and
  that was a negative result (keep 0.9/0.8/cov-mode 0, nothing proven optimal).
  Deliberately not started this session — each sweep point is itself a
  multi-hour cluster run, and the queue's already caused 4 failures; needs an
  explicit go-ahead given the compute cost, ideally after the corrected numbers
  from job 29026388 are in hand.
- **Wiring #130's real assembly-quality-QC output into #134's diagnostics
  framework** — done in NovInvenio PR #167 (open, not merged, 2026-09-23):
  `assembly_quality_confound` is now computed from
  `assembly_quality_correlations.tsv`. `clade_structure_validity` and
  `gain_loss_skew` are still `not_computed`. Not yet exercised on real data.
- **HET/Starship/co-gain-loss report needs regenerating** — `detect_trans_
  modules.py` is a standalone post-hoc script (not a Nextflow step), last run
  2026-09-18 against the *unrescued* matrix
  (`results/mmseqs_genus_vs_ureesii/output/pangenome/trans_modules/
  TRANS_MODULE_REPORT.md` — real findings, HET/NACHT domains in 312 islands,
  Starship DUF3435 recurring across the 3 largest trans-modules, an intact
  HR-PKS cluster moving across 315 islands — but every number needs
  re-deriving once 29026388's `pair_classification.tsv` exists).
- **No publish path from `pangenome.nf` reports to the public site
  (`nii.fungalgenomes.org`)** — `bin/sync_reports.sh`/`bin/generate_docs.py`
  only understand the *other* pipeline pathway (`main.nf`'s novelty-discovery/
  loss-search: `presence_matrix.function.tsv`, `novelties.html`/`core.html`/
  `losses.html`). `docs/fungi/{Afumigatus_pangenome,coccidioides_pangenome}/`
  are empty placeholder cards ("data staged, pipeline not yet run") in
  `docs/fungi/index.html` — `generate_docs.py`'s only test for "complete" is
  whether `docs/<domain>/<set>/report.html` exists, and nothing currently
  produces that file for a `pangenome.nf` run. Building that publish step is
  new work, not a config flip.
- **No Parquet anywhere in this codebase** — all tabular output is TSV, several
  already `.zst`-compressed (`cooccurring_pairs.tsv.zst`, `gene_positions.tsv.zst`,
  `family_positions.tsv.zst`). If Parquet distribution is wanted, that's new
  scope.

## 4. Quick reference

- Flagship study: `studies/fungi/coccidioides_pangenome/`, config
  `config_genus_vs_ureesii.csv` (529 ingroup + *U. reesii* outgroup).
- Pipeline repo: `/bigdata/stajichlab/jstajich/projects/NovInvenio` (main branch
  has everything in section 1). Live run's pipeline copy:
  `/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/structural-rescue-criterion`
  (pinned at `3277bdb` — do not `git pull` or recreate this until the run finishes).
- Prior (pre-fix) completed run for browsing figures/report shape right now:
  `studies/fungi/coccidioides_pangenome/results/mmseqs_genus_vs_ureesii/output/pangenome/`
  — numbers are stale, structure/method is representative.
