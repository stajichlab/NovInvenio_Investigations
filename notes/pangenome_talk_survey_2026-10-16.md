# Pangenome talk survey (2026-10-16)

Read-only survey. Scope: Afumigatus pangenome, C. immitis and C. posadasii pangenomes
(studies/fungi/coccidioides_pangenome), Afumigatus test sets, Bd_pangenome (listed only).
Every number has a source. "Not found" means I looked and it does not exist.
Abbreviation: `CP` = `studies/fungi/coccidioides_pangenome`, `AF` = `studies/fungi/Afumigatus_pangenome`.

## 1. Runs, status, current run

### Afumigatus_pangenome (AF)

- Registry: `AF/pangenome_runs.yaml`. Runs: `full_v070` (comment: "failed at ISLAND_LOCI"),
  `full_v071` ("baseline ... not published"), `full_v072` (`current: true`, `publish: true`,
  pipeline 05fcf3a, `pangenome_island_reference_strains: Asfu_Af293,Asfu_A1163`).
- Samplesheet `AF/config.csv`: 293 IN + 2 OUT (`cut -d, -f1 | sort | uniq -c`).
  OUT = `Aslen_ref` (A. lentulus), `Neofi_ref` (A. fischeri).
- Representatives (`full_v072/.../report_tables/per_strain_summary.tsv`, cols 2-4):
  IN 121 representative + 172 non-representative; OUT 2.
- Annotation: `AF/species.csv` Protein_Source: 292 `ncbi`, 3 `uniprot`.
- Status of full_v072: `AF/.nf_launch/full_v072/ni_run.json` `"exit_status": 0`, SLURM 29648677;
  `sacct`: COMPLETED 2026-10-08 15:53 -> 21:18 (5 h 25 min; resumed full_v071 cache);
  slurm-29648677.out: `[SUCCESS] completed=19 failed=0 cached=913`.
- Published: `docs/fungi/Afumigatus_pangenome/full_v072/run.json` (published 2026-10-11,
  n_families 47954, n_strains 295); `docs/fungi/Afumigatus_pangenome/study.json` current = full_v072.
- Older, superseded outputs (not from pangenome.nf): `AF/results/SUMMARY.md`, `REPORT.md`,
  `full_293run/`, `nf_latest/`, `rescue_pass/` (Sep 15-25, bin/ chain).
- Results dir of current run: `AF/results/full_v072/config/pangenome/`.

### Afumigatus test sets

| Set | Run | Genomes | Status | Published |
|---|---|---|---|---|
| Afumigatus_test45 | isl45_v1 | 45 IN + 2 OUT (`config.csv`) | sacct 29211980 COMPLETED 1:28:29 | yes, docs/fungi/Afumigatus_test45/isl45_v1 (n_families 29213, run.json) |
| Afumigatus_test45 | isl45_nopfam_v1 | same | results dir exists; not checked further | no |
| Afumigatus_test45 | genome-count-sweep | subsample analysis | published page | yes |
| Afumigatus_test10 | in10_v1 | 10 IN (memory) | completed (memory: 77/77 tasks), 0 islands | no |
| Afumigatus_test_c1c7 | c1_vs_c7_v1 | 5 IN vs 5 OUT (memory) | not launched (no results/ or .nf_launch/) | no |

isl45_v1 bins (`frequency_table.tsv` col 4, `cut -f4 | sort | uniq -c`): core 7071, soft_core 698,
shell 4348, cloud 4874, singleton 6386, outgroup_only 5836 (29213 families).

### coccidioides_pangenome (CP)

- Registry: `CP/pangenome_runs.yaml`.
  - immitis-in: `_v2` (failed, `ni_run.json` exit_status 1), `_v3` (baseline, not published),
    `_v4` (05fcf3a, published), `_v5` (`current: true`, a2836ad = NovInvenio PR #235 branch,
    `pangenome_stratum_source: mash_clade`, reference strains `RS,WA_211`).
  - `posadasii_in_immitis_out_v2`, `immitis_only_v2`, `posadasii_only_v2`: registered but
    **never launched** (no `CP/results/<run>` and no `CP/.nf_launch/<run>` directory).
- Samplesheets: `config_immitis_in_posadasii_out.csv` 169 IN + 360 OUT;
  `config_posadasii_in_immitis_out.csv` 360 IN + 169 OUT; `config_immitis.csv` 169 IN;
  `config_posadasii.csv` 360 IN.
- Representatives in v5 (`per_strain_summary.tsv`): IN (immitis) 133 rep + 36 non-rep;
  OUT (posadasii) 278 rep + 82 non-rep.
- Annotation: `CP/species.csv` Protein_Source 529 `local_faa`; README: funannotate annotation
  freeze `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation_freeze/20260112`.
- Status of v5: `ni_run.json` exit_status 0, SLURM 29724129; sacct COMPLETED
  2026-10-09 09:13 -> 2026-10-10 06:20 (21 h 07 min); `[SUCCESS] completed=28 failed=0 cached=1602`.
- Published: `docs/fungi/coccidioides_pangenome/immitis_in_posadasii_out_v5/run.json`
  (2026-10-11, n_families 47780, n_strains 529); `study.json` current = immitis_in_posadasii_out_v5.
- v4 and v5 `frequency_table.tsv` are byte-identical (`cmp`). v5 differs only downstream
  (stratified co-occurrence null).
- **The only posadasii-in reciprocal run with output** is `rescue_freqpol_posadasii_in_immitis_out`
  (`CP/results/rescue_freqpol_posadasii_in_immitis_out/output/pangenome/`), pipeline 91e3157
  (`CP/.nf_launch/reciprocal_posadasii_in_immitis_out/slurm-*.out`: "pipeline commit: 91e3157..."),
  launched before `ni`; report re-rendered at f87fd1e for #212/#214 (`docs/.../run.json`
  `report_rerender`). Published 2026-09-29. No species tree, no Dollo columns, no island
  locus view / clinker.
- Superseded / retired: `mmseqs_immitis`, `mmseqs_posadasii`, `mmseqs_*_in_*_out`, `mmseqs_wholeset`
  (Sep 15-17, before the 2026-09-21 strain corrections in `CP/CORRECTIONS.md`); `*ureesii*`,
  `tier1_sweep` (U. reesii excluded, README "Decision: U. reesii excluded (2026-09-28)");
  `cocci_locusview_preview` (hand-built preview).

### Bd_pangenome (listed only)

Same pipeline type. 344 IN + 2 OUT (`config.csv`). `full_v1` (05fcf3a) published at
docs/fungi/Bd_pangenome/full_v1; its slurm .out ends `[FAILED] completed=1437 failed=7`
while `ni_run.json` says exit_status 0 (not investigated). `full_v2` is marked `current: true`
but has no results/.nf_launch directory (not launched).

## 2. Output files of the current runs

Paths relative to `AF/results/full_v072/config/pangenome/` (AF) and
`CP/results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome/` (v5).

| File | AF size | Header / content |
|---|---|---|
| frequency_table.tsv | 2.5 MB, 47955 lines | `family frequency strain_count bin frequency_out strain_count_out bin_out` |
| presence_matrix.tsv / presence_matrix.rescued.tsv | 99 MB / 100 MB | families x strains |
| cluster/tier1_cluster.tsv, tier1_rep_seq.fasta | 149 MB / 23 MB | mmseqs clusters |
| cooccurring_pairs.tsv.zst | 21 MB | `family_a family_b jaccard fisher_p fdr_q permutation_p direction_a clade_composition asymmetry_a direction_a_freq permutation_q` (+ v5 only: `direction_a_tree gain_node n_loss_events loss_clades asr_method`) |
| pair_classification.tsv.zst | 22 MB | pair classes |
| significant_islands.tsv, island_regions.tsv, island_loci.json | 2.6 MB / 0.45 MB / 2.5 MB | islands |
| report_tables/per_strain_summary.tsv | 18 KB | `Short group is_representative bins_from n_families core soft_core shell cloud singleton nonrep_only outgroup_only singleton_z is_outlier` |
| report_tables/group_class_overlap.tsv | AF 0 bytes; v5 1010 B | `ingroup_class outgroup_class n_families` |
| report_tables/classification_counts.tsv, island_size_distribution.tsv | small | counts |
| report/pangenome_openness.tsv | 90 B | `kappa gamma r_squared is_open core_inf tau` |
| report/report.md, report/figures/, report/figures_pdf/ | | rendered report |
| rescue_funnel.tsv, diagnostics/diagnostics.tsv | small | rescue and QC diagnostics |
| clinker/ | 105 MB | 100 loci, 25 published |

**Gain/loss output:** the pipeline writes polarity only per co-occurring pair, for
`family_a` (`direction_a` strict, `direction_a_freq` with `pangenome_polarity_loss_min_frac` 0.9 /
`gain_max_frac` 0.1, and with a species tree `direction_a_tree` by Dollo). There is **no
family-level gain/loss table for all families**, and the diagnostic `gain_loss_skew` is
"not computed" in both runs (`diagnostics/diagnostics.tsv`).

## 3. Headline numbers

### Composition (pipeline categories; bins over representative strains of each group)

Source: `report/report.md` "Pangenome composition", and `tail -n +2 frequency_table.tsv | cut -f4 | sort | uniq -c` (col 7 for `bin_out`).

| Class | A. fumigatus (AF v072, IN, 121 reps) | C. immitis (v5 IN, 133 reps) | C. posadasii (v5 OUT `bin_out`, 278 reps) | C. posadasii (posadasii-in run, IN) | C. immitis (posadasii-in run, OUT `bin_out`) |
|---|---|---|---|---|---|
| total families in table | 47954 | 47780 | 47780 | 47719 | 47719 |
| families in the group (sum of 5 classes) | 36302 | 26047 | 32609 | 32576 | 26027 |
| core (>=0.95) | 7165 (19.7%) | 6846 (26.3%) | 6357 (19.5%) | 6361 | 6848 |
| soft_core (>=0.90) | 667 (1.8%) | 610 (2.3%) | 636 (2.0%) | 640 | 604 |
| shell (>=0.15) | 4612 (12.7%) | 4025 (15.5%) | 4628 (14.2%) | 4624 | 3998 |
| cloud (<0.15) | 10561 (29.1%) | 6764 (26.0%) | 11894 (36.5%) | 11907 | 6772 |
| singleton (1 rep) | 13297 (36.6%) | 7802 (30.0%) | 9094 (27.9%) | 9044 | 7805 |
| nonrep_only | 5814 | 1406 | 3591 | 3580 | 1413 |
| outgroup_only / ingroup_only | 5838 (outgroup_only) | 20327 (outgroup_only) | 11580 (ingroup_only) | 11563 (outgroup_only) | 20279 (ingroup_only) |

Note: in the posadasii-in run the outgroup (immitis) bins are core 6848, soft_core 604, shell 3998,
cloud 6772, singleton 7805 (sum 26027). These agree with v5's immitis IN bins within ~1%,
but the runs used different pipeline commits (91e3157 vs a2836ad).

Afumigatus outgroup is not binned: 2 strains < `pangenome_outgroup_min_bin_strains = 3`
(report.md line 59: "the outgroup is not binned (too few representative genomes)").

### Immitis vs posadasii cross-table (v5 `report_tables/group_class_overlap.tsv`)

Examples: core/core 5424; immitis core but absent in posadasii 296; immitis singleton absent in
posadasii 5386; absent in immitis but posadasii cloud 7663, singleton 7551, core 199.
Full 49-row table in the file.

### Openness (Heaps' law, `report/pangenome_openness.tsv`)

| Run | kappa | gamma | R2 | open | core_inf |
|---|---|---|---|---|---|
| AF full_v072 | 7492.08 | 0.3171 | 0.980 | True | 2423 |
| CP v5 | 9019.95 | 0.2627 | 0.9974 | True | 295 |
| CP posadasii-in | 9019.05 | 0.2626 | 0.9975 | True | 295 |

Report text: the curve "counts every strain in the presence matrix, including near-identical
strains and the outgroup". So the Cocci curve mixes both species. A per-species curve is not found.

### Per-strain singleton counts (`per_strain_summary.tsv` col 10, my awk min/median/max)

| Run | group / rep | n | min | median | max |
|---|---|---|---|---|---|
| AF v072 | IN rep | 121 | 113 | 272 | 2203 |
| AF v072 | IN non-rep | 172 | 120 | 218 | 1442 |
| AF v072 | OUT | 2 | 471 | - | 742 |
| CP v5 | IN (immitis) rep | 133 | 178 | 270 | 779 |
| CP v5 | OUT (posadasii) rep | 278 | 293 | 400 | 1094 |
| CP posadasii-in | IN (posadasii) rep | 278 | 292 | 398 | 1090 |
| CP posadasii-in | OUT (immitis) rep | 133 | 183 | 268 | 785 |

Families per strain: AF min 9525, median 10365, max 11911 (n=295); CP v5 min 8972,
median 9586, max 10324 (n=529) (report.md "Per-strain summary").

### Gain / loss (per unique `family_a` in `cooccurring_pairs.tsv.zst`)

Command: `zstd -dc f | tail -n +2 | awk -F'\t' -v c=COL '{print $1"\t"$c}' | sort -u | cut -f2 | sort | uniq -c`.
Scope: only families that appear as `family_a` in an FDR-significant co-occurring pair, not all families.

| Run | pairs | unique family_a | direction_a (strict) | direction_a_freq (0.9/0.1) | direction_a_tree (Dollo) |
|---|---|---|---|---|---|
| AF v072 | 1,363,573 | 8232 | gain 6575, loss 874, ambiguous 783 | identical to strict | not computed (no species tree) |
| CP v5 (immitis in) | 77,059 | 4161 | gain 1556, loss 12, ambiguous 2593 | gain 2726, loss 357, ambiguous 1078 | loss 4161 (all) |
| CP posadasii-in (91e3157) | 157,454 | 7628 | gain 3623, loss 86, ambiguous 3919 | gain 5561, loss 660, ambiguous 1407 | not computed |

Warning on the v5 Dollo column: all 4161 families are "loss"; `gain_node` is
`MRCA(1M0..h5384)` for 4113 of them. Why is not checked. Do not use the Dollo column on a
slide until it is understood. For AF, strict and freq polarity are the same because the
outgroup has only 2 strains.

Per-branch / per-clade gain/loss counts: not found as a table. The pair table has
`clade_composition` (carriers per clade) and, in v5, `loss_clades`/`n_loss_events`.

### Other

- Islands: AF v072 11676 accessory islands (6758 single-strain), report.md line 63-70.
- Pair classes: AF trans_unconfirmed 616750, trans 543757, insufficient_data 188361,
  unexplained_physical 10264, ambiguous_linkage 4441; CP v5 trans 65147, unexplained_physical 5048,
  trans_unconfirmed 3316, ambiguous_linkage 2009, insufficient_data 1539 (`classification_counts.tsv`).
- Rescue funnel: AF applied 288996 cells; CP v5 applied 571866 (`rescue_funnel.tsv`).

## 4. Figures that exist and gaps

Published per run in `docs/fungi/<study>/<run>/figures/` (PNG) and `figures_pdf/` (PDF):

- AF full_v072: accumulation_curve, frequency_distribution, group_composition,
  island_domain_enrichment, island_size_distribution, pair_classification_summary,
  per_genome_class_composition, presence_absence_matrix; plus `locus_summary/`
  (locus_overview.png, locus_strips_selected.png, locus_strips.pdf, starship_overlay.png).
- CP v5, v4, rescue_freqpol_* : the same 8 plus group_class_overlap_heatmap,
  group_class_overlap_shared, group_class_overlap_upset; v5 has `locus_summary/` (no starship overlay).

| Needed for talk | Status |
|---|---|
| Cluster counts + size distribution | exists per run: `frequency_distribution.png` (ingroup only, fraction of strains, coloured by class) |
| Core/accessory bar | exists: `group_composition.png`; labels say "ingroup"/"outgroup", not species names |
| Accumulation curve | exists: `accumulation_curve.png`; mixes all strains incl. outgroup; no per-species curve |
| Gain vs loss summary by direction | **missing** (numbers only in pair table) |
| Gain vs loss by clade | **missing** |
| Per-strain / per-species singleton counts | partial: `per_genome_class_composition.png`; no per-species singleton box plot |
| Cross-study comparison (3 analyses side by side) | **missing** |
| Summary table section in report | **missing**: composition is a bullet list in report.md; report.html has 3 tables (top islands, trans modules, Pfam), none summarise composition |

## 5. Method comparability

| Item | AF full_v072 | CP v5 (immitis-in) | CP posadasii-in (rescue_freqpol) |
|---|---|---|---|
| Pipeline commit | 05fcf3a | a2836ad (PR #235 branch) | 91e3157 (+ report steps re-rendered at 8c55750, f87fd1e) |
| Clustering | mmseqs 0.9 id / 0.8 cov (defaults in nextflow.config lines 238-240, same at 91e3157) | same | same |
| Rescue | on, per-strain DB, max_target_seqs 200 | same | same (91e3157 nextflow.config line 249) |
| Polarity thresholds | 0.9 / 0.1 (default) | same | same |
| Species tree / Dollo | none | coccidioides_coccidioides_only.full.nwk | none |
| Stratum | taxon_group (default) | mash_clade | taxon_group |
| Outgroup | 2 other species (A. lentulus, A. fischeri) | 360 C. posadasii | 169 C. immitis |
| Genomes | 293 + 2 | 169 + 360 | 360 + 169 |
| Annotation | NCBI proteins/GFF3 (292) + UniProt (3) | funannotate freeze 20260112 | same as v5 |
| Diagnostics | rescue_redundancy triggered (87.8%); assembly confound OK (rho 0.08) | rescue_redundancy triggered (77.5%); assembly_quality_confound **triggered** (rho -0.44, n=169) | not checked |

Known issues from memory:

- Rescue-pass bugs (memory `afumigatus_pangenome_rescue_pass_bug.md`): "Both bugs below were
  found and fixed on 2026-09-15". The fix was in the old bin/ chain. full_v072 runs the
  pipeline's per-strain DB design (`pangenome_rescue_max_target_seqs = 200 // per-strain DB is tiny`),
  so the truncation bug does not apply. Its "Final corrected numbers" (core 53.0%, singleton 22.6%,
  `AF/results/SUMMARY.md`) bin over all 293 strains and are **not comparable** with v072
  (bins over 121 representatives; core 19.7%, singleton 36.6%). Do not mix them on slides.
- #212 singleton bin (memory `pangenome_singleton_bin_issue.md`): "Relaunches under new names
  at a commit >= f87fd1e carry #212 and #214 natively". full_v072 (05fcf3a) and v5 (a2836ad)
  were launched after that, and their tables have `nonrep_only`/`outgroup_only` bins and
  `bin_out` columns, so the fix is in. The posadasii-in run got it by offline re-render
  (run.json `report_rerender`, 8c55750).
- Cocci primary comparison (memory): "reciprocal immitis vs posadasii is the main comparison";
  "U. reesii ... excluded".
- Island reference strains (memory `island_reference_strain_defaults.md`): Afum `Asfu_Af293,Asfu_A1163`;
  immitis-in `RS,WA_211` (confirmed); posadasii-in `B3245,COCPO_Galgiani`.
- Test sets (memory): in10_v1 0 islands; isl45_v1 2,859 islands; c1c7 not launched.
- v5 report: `assembly_quality_confound` triggered: "Between-strain accessory-content comparisons
  are confounded by assembly fragmentation on this dataset" (`diagnostics/diagnostics.tsv`).

## 6. Open questions for the user

1. Cocci "posadasii pangenome": use the posadasii OUT bins of v5 (same run, same commit as immitis),
   or the older posadasii-in run (91e3157), or launch `posadasii_in_immitis_out_v2` (registered,
   never run) with reference `B3245,COCPO_Galgiani`? Per-species gain/loss polarity for posadasii
   needs a posadasii-in run.
2. Do you want species-only runs (`immitis_only_v2`, `posadasii_only_v2`, registered, not run) for
   per-species accumulation curves, or a per-group curve computed from the v5 matrix?
3. Polarity labels: strict (`direction_a`), frequency (`direction_a_freq`, 0.9/0.1) or Dollo?
   The v5 Dollo column calls every family "loss"; needs a check before use.
4. Gain/loss scope: only co-occurring `family_a` families exist now. Do you want polarity for all
   accessory families (new computation from presence matrix + outgroup)?
5. Afumigatus: is a 2-species outgroup polarity acceptable on a slide, given strict = freq here?
6. Which numbers count as "clusters": total families in table (47954 / 47780) or families present
   in the group (36302 / 26047 / 32609)? With or without nonrep_only?
7. Afumigatus SUMMARY.md (Sep 15, bin/ chain) vs full_v072: confirm the old numbers are retired.

## 7. Proposed slide-ready outputs (not built)

| # | Output | Data source | Work needed | Effort |
|---|---|---|---|---|
| 1 | Summary table: genomes, reps, families, core/soft/shell/cloud/singleton per species (A. fum, C. imm, C. pos) | 3 `frequency_table.tsv` (+`bin_out`), `per_strain_summary.tsv`, `pangenome_openness.tsv` | new `bin/pangenome_talk_summary.py` (reads N run dirs, writes TSV + PNG table) | ~2 h |
| 2 | Side-by-side composition bars (3 species, % of families) | same as 1 | same script, one matplotlib figure | ~1 h |
| 3 | Frequency histograms per species (immitis, posadasii from `frequency_out`) | `frequency_table.tsv` cols 2, 5 | same script | ~1 h |
| 4 | Per-species accumulation curves | presence_matrix.rescued.tsv split by group | script with permutations; or launch species-only runs | ~3 h script; or 1-2 runs |
| 5 | Gain vs loss bar per study (strict and freq) | `cooccurring_pairs.tsv.zst` per unique family_a | same script | ~1 h |
| 6 | Gain/loss by clade | `clade_composition` column (and v5 `loss_clades`) | script + decision on clade labels; check Dollo first | ~3-4 h |
| 7 | Per-strain singleton box plot by species | `per_strain_summary.tsv` col 10, reps only | same script | ~1 h |
| 8 | Immitis x posadasii shared/specific families | existing `group_class_overlap_shared.png` / upset in v5 | none (relabel ingroup/outgroup to species if wanted) | 0-1 h |

Placement: a general multi-run summary script fits `bin/` (takes run dirs as args); study-specific
labelling can live in `studies/<d>/<set>/bin/`. Outputs from 1-7 should go to an analysis dir
(e.g. `studies/fungi/<set>/analysis/talk_2026/`), with provenance per CLAUDE.md before commit.
