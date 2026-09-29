# Session handoff — 2026-09-28 (#212, #214, Afumigatus test sets)

Companion to `2026-09-28-session-handoff.md`, written by a parallel session on the same day
(island locus view, SANS, #132 notes, cleanup). This session fixed report figures,
found and fixed a mislabelled frequency class (NovInvenio #212), republished the
Coccidioides runs, found why both v2 full runs failed, and set up small
*A. fumigatus* test sets to iterate on before any more full runs.

## 1. Code landed in `NovInvenio` (all merged to main)

| PR / commit | What |
|---|---|
| #205 | Island-size plot: numeric x axis, a tick for every size up to 10, then every 5. |
| #206 | Top islands table: each Pfam domain links to its InterPro/Pfam page (URLs from `island_pfam_enrichment.tsv`). |
| #207, #210 | Report breadcrumbs: broadest first (All studies > Group > Study > Run), paths fixed for the per-run layout `docs/<domain>/<study>/<run>/`. |
| #211 | Per-genome class figure: one stacked bar per genome (core ... singleton), ordered GROUP > Species > total; names up to 150 genomes. |
| `8c55750` (#212 PR 1, merged directly) | Per-group frequency bins. `frequency_table.tsv` keeps `bin` (ingroup) and adds `frequency_out` / `strain_count_out` / `bin_out`; new classes `nonrep_only`, `outgroup_only`, `ingroup_only`, `absent`; param `pangenome_outgroup_min_bin_strains` (3). `per_strain_summary.tsv` gains `group`, `is_representative`, `bins_from`, `nonrep_only`, `outgroup_only`; outlier z within group over representatives. Pie replaced by per-group composition bars. Spec/plan: `docs/superpowers/{specs,plans}/2026-09-28-pangenome-group-bins*`. |
| #214 (`f87fd1e`, #212 PR 2) | `report_tables/group_class_overlap.tsv` (49 cells) + "Ingroup vs outgroup content" section: class x class heatmap, shared / ingroup-only / outgroup-only bar, matplotlib UpSet. Spec/plan: `docs/superpowers/{specs,plans}/2026-09-28-pangenome-group-overlap-figures*`. |
| #213 (`cb5d027`, the parallel session) | Island steps (`ISLAND_DNA_TARGETS`, `ISLAND_LOCI`, `ISLAND_SYNTENY`) take the merged Pfam table, not the per-chunk channel. This is the v2 failure fix (section 3). |

#212 is closed. Why it mattered: the `singleton` class held families absent from
the ingroup and families present only in dereplicated near-duplicate genomes.
In `rescue_freqpol_immitis_in_posadasii_out`, 21,694 of 29,502 "singletons" had
ingroup count 0; outgroup genomes' median "singleton" count was 1,582 (ingroup
37). After the fix: `singleton` 7,808, `outgroup_only` 20,296, `nonrep_only`
1,398; outgroup per-genome singleton median 25. No family counted in the ingroup
changed class; co-occurrence, islands and enrichment inputs are unchanged.

## 2. `NovInvenio_Investigations` changes

- Site: "studyies" -> "studies" on the gallery cards (`lib/site_pages.py`).
- U. reesii excluded from pangenome work; recorded in the Coccidioides study
  README; `run_tier1_sweep.sh` marked retired. No data deleted.
- The 3 published Coccidioides runs were updated offline (report steps only, no
  pipeline rerun), in two steps:
  - #212 PR 1 (`889abdf`): frequency bins -> REPORT_TABLES -> assembly QC -> DIAGNOSTICS -> render.
  - #214 (`83c164a`): REPORT_TABLES + render.
  - Backups: `results/<run>/output/pangenome/pre212_20260928/` and `pre214_20260928/`.
  - `run.json` records each update under `report_rerender` (sync support: `c8788c9`).
- `bin/pangenome_outgroup_polarity.py` reads `frequency_table.tsv` by column name (`cd2adf4`).
- `bin/sync_pangenome_report.py` archives `group_class_overlap.tsv` and skips 0-byte tables (`9d688bd`).
- Published pages keep hand-patched `island_synteny.html` breadcrumbs (release
  asset). A re-sync copies the run's own file back; re-apply the text patch, or
  regenerate the page with the current pipeline.

## 3. Failed runs and the cause

Both v2 full runs at NovInvenio `a394ace` (v0.7.0) failed on 2026-09-28 at the same step:

| Run | Ran | Failure |
|---|---|---|
| Afumigatus `full_v070` (SLURM 29171115) | 18 h 17 min, 1,248 tasks done | `ISLAND_LOCI` input file-name collision on `dna_calls_batch_00N.tsv` |
| Coccidioides `immitis_in_posadasii_out_v2` (29171116) | 20 h 44 min, 2,166 tasks done | the same |

Cause (confirmed from the Cocci work dirs): at `a394ace` the island steps got
`FAMILY_PFAM_SCAN.out.domtblout`, one file per Pfam chunk. With 4 chunks,
`ISLAND_DNA_TARGETS` ran 4 times, each writing `batch_001..003.tsv`; the collect
into `ISLAND_LOCI` then collided. Fixed by #213 (`cb5d027`).

Relaunch constraint: `bin/ni pangenome run` refuses to reuse a run name at a
different commit. A relaunch on the fixed commit needs new run names, so it starts
from zero (no Nextflow cache). The other three Cocci v2 runs were never launched.

## 4. Live run — small test set

Plan agreed with the PI: test on a small, informative set first, then relaunch the full runs.

- **Set 1** `studies/fungi/Afumigatus_test10`, run `in10_v1`: 10 *A. fumigatus*, all `IN`,
  no outgroup. Af293, A1163, the fewest-contig genome (BUSCO >= 98.9%) of each Barber
  cluster 1-7, and a second genome of cluster 1.
  - Launched 2026-09-28 at NovInvenio `f87fd1e`.
  - Head job **29196843** (exfab). At handoff: 31 tasks completed, 0 failed.
  - `CLUSTER_TIER1` (**29196894**) pending in `preempt`.
  - First real run with no outgroup: watch for any step that assumes one.
- **Set 2** `studies/fungi/Afumigatus_test_c1c7`, run `c1_vs_c7_v1`: 5 Barber cluster 1 (IN)
  vs 5 cluster 7 (OUT). These are the most distant clusters: mean mash 0.0075 between,
  0.0009 / 0.0021 within. Ready; not launched.
- Both sets:
  - Genome choice by `Afumigatus_test10/bin/make_test_subsets.py`; provenance in each `DATA_MANIFEST.yaml`.
  - Data read from `Afumigatus_pangenome/data_dir`; nothing is copied.
  - No `publish.yaml`. Commit `b3444c7`.
- Nextflow GitHub auth: `~/.nextflow/scm` now holds the `gh` token (user hyphaltip, mode 600).
  The first launch had failed on the unauthenticated API rate limit in `nextflow pull`.

## 5. Next steps

1. Watch `in10_v1` to the end. Check the report, the per-genome figure, the island/locus view
   and clinker pages, and the "outgroup not binned" skip lines.
2. Run set 2 (`c1_vs_c7_v1`). Check per-group bins, `bin_out` and the overlap figures on a case
   where the expected pattern is clear.
3. Only then relaunch the full Afumigatus and Coccidioides runs under new run names at the
   tested commit. Those runs will carry #212 and #214 natively, so no offline update is needed.
4. After the new Cocci runs are published, retire `cocci_locusview_preview`. It is the only
   published run with the locus view; the two `rescue_freqpol_*` runs (pipeline `91e3157`)
   predate it.

## 6. Deferred minors (from the two branch reviews)

- The report's samplesheet reader reads raw `GROUP` (no alias mapping or strip) for axis labels
  and per-genome block order.
- No test runs `pangenome_report_tables.main()` with the new `--samplesheet` / `--strain_inventory`.
- Heatmap colour scale starts at the smallest non-zero value, so a count-1 cell is near-white.
- A header-only overlap file draws empty figures (REPORT_TABLES never writes one).
- UpSet: empty space above the set panel. Shared bar: nested parentheses in the label.
- `OVERLAP_CLASSES` is defined in both the tables and the render script.
- Cross-group dereplication: an ingroup genome deduplicated against an outgroup representative
  leaves its own group's denominator (predates #212; not in any spec).
