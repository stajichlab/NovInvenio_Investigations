# Phase 2 status: validation battery and cross-method support

Date: 2026-10-07. Source todo: `NovInvenio/todo/cross-method-support-column.md`. Source ADR: `NovInvenio/docs/adr/0002-family-profile-search-pathway.md` (Q8).

This note records what is done and what remains. A box is checked only where I found evidence in code, an issue, a note or a result file. "Not verified" means I did not check it. I did not rerun any tool for this note.

Revised 2026-10-07 after an independent Opus review. Corrections I checked myself are marked (verified). Findings I took from the review without checking the code are marked (review only).

Where the raw sweep results are: `/bigdata/stajichlab/jstajich/projects/NI_Sweep/` (`pezizo5_coverage/`, `pezizo5_coverage_ext/`, `sordariales_shallow/`, `deep_broad_1kfg/`). My first draft looked only under `NovInvenio/results` and wrongly said the raw results were missing.

Phase 2 has three parts, matching the closed issues #5 (controls scorer), #6 (sweep) and #7 (support column).

## 2a. Controls scorer (issue #5, closed)

- [x] `bin/score_controls.py` exists (436 lines). It scores family mode (`--cluster_tool mmseqs`) and pairwise mode.
- [x] It was used for real: per-control tables for three clades are in `studies/fungi/pezizo_set1_cluster/tier_comparison/` and `studies/fungi/sordariales_shallow_cluster/tier_comparison/`. The summary is in `README.md` of this folder.
- [x] Control sets exist for six clades in `NovInvenio/configs/controls/`:

  | Clade | Positives | Negatives | Notes |
  |---|---|---|---|
  | pezizo_set1 | 6 | 10 | plus 1 `EXAMPLE_` placeholder |
  | pezizo5 | 6 | 11 | plus 1 placeholder; the same genes as pezizo_set1, ID-translated |
  | sordariales_shallow | 1 | 8 | the 8 are BUSCO |
  | Agaricales | 2 | 5 | |
  | Mucorales | 2 | 0 | |
  | Chaetothyriales | 0 | 0 | all 4 rows are `EXAMPLE_` placeholders |

  Counts (verified) come from parsing each CSV with Python's csv module and dropping `EXAMPLE_` rows. pezizo_set1 and pezizo5 share their positives, so the 9 resolved positives in the tier comparison are about 7 independent genes or fewer (review only).
- [x] Result so far (from `README.md`): Tier P (pairwise diamond) recovered 6/6 resolved pezizo_set1 positives and had 0/8 false positives on sordariales_shallow BUSCO negatives. Tier C+H (mmseqs plus family HMM) recovered 2/5 pezizo_set1 positives.
- [ ] Pooled recall estimate. The three clades give only 9 resolved positives in total, which is too few. The existing note states this and does not pool.
- [ ] Negatives resolved for pezizo_set1 (0/10 resolved in the tier comparison). Agaricales has 5 negatives in `configs/controls`, but the tier comparison reports no negatives for agaricomycetes, so they were not scored there.
- [ ] Loss-direction positive controls. None exist in any clade.
- [ ] `fasta` anchors in pairwise mode. The scorer reports them unresolved.
- [ ] Circularity in the positives (review only): the positives were screened with `--very-sensitive` diamond ("zero hit in any outgroup"), so Tier P recovering them is partly true by construction. "Tier P is the reference" is therefore not established by the controls alone.
- [ ] Leakage (review only): `hmm_presence_domain_evalue = 1e-5` was chosen to fix ada-1 and ham-5, which are also scoring positives. Tuning on a control and then reporting recall on it is not a held-out test.
- [ ] Fix stale docs: `configs/controls/README.md` still says the scorer is "not yet built".

## 2b. Parameter sweep and BUSCO battery (issue #6, closed)

Issue #6 was closed, but its checklist in the issue body is unchecked. Its scope was larger than what was run (see below).

Built:
- [x] `bin/run_param_sweep.sh` (now in `NovInvenio_Investigations/legacy/novinvenio_scripts/`), `bin/collate_sweep.py` (225 lines), BUSCO clustering recovery (`bin/busco_family_recovery.py`) and BUSCO presence recovery (`bin/busco_presence_recovery.py`).
- [x] `collate_sweep.py` has a documented knee rule: quality gates, a composite score, ties broken to fewer novelties. It skips a metric that was never measured instead of treating blank as worst.
- [x] The metrics include number of families, number of novelties, BUSCO recovery, presence recovery, recall, FP rate, and the TBLASTN-removal count (`tblastn_removed`).

Run (from `.living/decisions.md` #12, #13, #17 and `todo/validate-hmm-presence-coverage-broader-sweep.md`):
- [x] pezizo5 coverage grid: `hmm_presence_cov {0.5,0.3,0.2,0.1} x hmm_presence_min_residues {0,100}`, 8 points (verified in `NI_Sweep/pezizo5_coverage_ext/.../sweep_metrics_merged_all8.tsv`). My first draft said 4 points. The 0.2 / 0 point produced 0 novelties and 0 presence recovery but `run_ok=1`, so it is a failed point that the flag missed. Among the other points, 0.3 / 100 has presence recovery 0.995 (shipped default then 0.5 / 0: 0.960) and TBLASTN-removed 21 (then 554).
- [x] The same kind of grid on sordariales_shallow. The review reports 2 of the 4 rows have `run_ok=0` (0.5 / 0 and 0.5 / 100), so only 2 rows are valid. Decision #17 says 1 failed. I did not recheck this file. The claim that the ranking agrees with pezizo5 rests on very few valid points.
- [x] A `min_domain_evalue` sweep on pezizo_set1_cluster. Commits `80df78a`, `c603953`, `d077621`.
- [x] Defaults changed: `hmm_presence_cov` 0.5 -> 0.3 (`cd297ba`, 2026-09-22) and `hmm_presence_domain_evalue` set to 1e-5. `nextflow.config` shows `hmm_presence_cov = 0.3`, `hmm_presence_min_residues = 100`, `hmm_presence_domain_evalue = 1e-5`.
- [x] deep_broad_1kfg (130 taxa) attempted. The review says `sweep_metrics.tsv` has only 1 row and it has `run_ok=0`. The `busco_recovery` of 0.011 therefore comes from a failed run and is weak evidence that clustering fails at scale (review only). Decision #17 says not to finish its remaining points without a clustering sweep.

Not done:
- [ ] **The family-definition grid was never run**: `min-seq-id {0.2,0.3,0.4,0.5} x cov {0.5,0.7} x E {1e-3,1e-5}`. Decision #17 states clustering was held at `--min_seq_id 0.3 -c 0.8` for every sweep run. `family_min_seq_id = 0.3` and `family_cov = 0.8` are still the ADR guess. The only evidence is `busco_recovery` 0.753 on pezizo5 and 0.011 on deep_broad_1kfg. This is the main gap.
- [x] **Correction (verified): recall and fp_rate were combined into a sweep table once.** In `NI_Sweep/pezizo5_coverage_ext/.../sweep_metrics_merged_all8.tsv`, 3 of 8 rows have real values: 0.2 / 100 gives recall 1.0, FP 0.0; 0.1 / 0 and 0.1 / 100 give recall 0.75, FP 0.0. `sweep_scores_all8.tsv` is a `collate_sweep.py` output with 0.2 / 100 marked `chosen=1`. The todo says the user overruled it and the shipped default stayed at 0.3 / 100. Decision #17 is stale on this point.
- [ ] **That selection is not trustworthy as it stands.** The other 5 rows have blank recall and FP (controls were not scored for them). In a mixed column `collate_sweep.py` treats a blank as recall 0 and FP 1.0 (review only; the data show this: 0.3 / 100 has recall 0.0, FP 1.0 and `admissible=0`). So the "chosen" point won against rows that were never measured. Controls must be scored for every grid point before any knee is meaningful.
- [ ] A knee-chosen default for the family-definition parameters in `nextflow.config` and the ADR.
- [ ] The residue axis of the broader grid. The todo suggests `min_residues {0,50,100,150,200}`. Only {0,100} was run. The cov axis was extended to {0.5,0.3,0.2,0.1}.
- [ ] A per-length-bucket breakdown of presence recovery in the sweep table. Not verified.
- [ ] Known defects in `collate_sweep.py` and `score_controls.py` (review only; I did not read the lines): fallback pool includes `run_ok=0` rows; the default BUSCO gate (0.9) rejects pezizo5 (0.753) and sordariales (0.849), so the all8 "admissible" flags needed an unrecorded override; the docstring contradicts itself on blanks; unresolved controls are left out of the recall denominator, so a setting that loses a positive's family can score better; predicate thresholds (0.75 / 0.0) are CLI defaults, not read from the run; `run_ok` records only the Nextflow exit code. These need checking and fixing before any rerun.
- [ ] A third clade for the `hmm_presence_cov` result. deep_broad_1kfg was unusable.

Related, done outside Phase 2:
- [x] Diamond sensitivity benchmark (`notes/diamond-sensitivity/`). `--diamond_sensitivity very-sensitive` is the default (issue #171).
- [x] Cross-method comparison of Tier P and Tier C+H on three clades (`README.md` of this folder): genome-wide agreement is poor (Jaccard 0.08-0.21), and Tier C+H cost 465x to 1215x the CPU-hours of Tier P at these sizes.

## 2c. Cross-method support column (issue #7, closed)

- [x] `lib/report_data.py` has a `support` field per row (`'pairwise'`, `'mmseqs'`, `'pairwise+mmseqs'`). `bin/make_report.py` takes `--support_matrix` and `--support_method`. The novelties page has a "Concordant (both methods)" filter. A test exists (`tests/test_report_data.py`, support tests).
- [ ] **Not wired into the pipeline.** No workflow or `main.nf` passes `--support_matrix`. The column works only if someone runs `make_report.py` by hand with a second pathway's matrix. A workflow runs one `--cluster_tool` at a time, so a design is needed: combine two finished runs after the fact, or run both pathways in one invocation.
- [ ] Support column for the **core** and **losses** reports. The #7 closing comment defers it. The losses payload and template have no `support` field.
- [ ] Support for a third method (`orthofinder`, ADR Phase 4). Not started and not expected yet.
- [ ] The support side uses `derive_novelties`, which skips TBLASTN, while the primary side may be TBLASTN-filtered, so "concordant" may compare two different predicates (review only; `report_data.py` near line 370).
- [ ] Reconsider what the column means. The todo says multi-method candidates are "high confidence". The measured agreement is low (Jaccard 0.08-0.21), and 82-98% of Tier C+H extras are candidates where default diamond found few or no seed-group homologs. A concordant filter would keep a small subset. The data do not say which method is right when they disagree.

## Cross-cutting open items

- [ ] `todo/cross-method-support-column.md` still says "open". Its three acceptance boxes are unchecked. Update it to match this note.
- [ ] A third method to decide the HMM-outgroup misses that have no raw diamond hit (`README.md`, "Open items").
- [ ] The pezizo_set1 Tier P search cost is unmeasured.
- [ ] Cost scaling beyond about 12 proteomes is unmeasured, so no crossover between Tier P and Tier C+H can be stated.
- [ ] NovInvenio issue #137 (three-way comparison of `pairwise`, `mmseqs` and `novelty_discovery`) is open and overlaps this work. It should absorb the unfinished comparison items.
- [ ] `mmseqs` and `famsa` crash on this cluster's non-AVX2 nodes (issues #148, #181). Any family-definition sweep must run on SLURM nodes pinned to supported CPUs, which the existing slurm config already does.

## Proposed order (my suggestion, not decided)

0. Before any new compute: fix the `collate_sweep.py` and `score_controls.py` defects listed above, and decide how controls are split into selection and held-out sets. With about 7 independent positives, selection should rest on the curation-free metrics (BUSCO clustering recovery, BUSCO presence recovery, TBLASTN-removed). Controls would then be used once, as a held-out check. Spiked or simulated controls (ingroup genes with known outgroup orthologs removed) would give a larger test set (review's suggestion).
1. Housekeeping: update the todo, the controls README, decision #17 and issue #6's checklist so they match this note. No compute needed.
2. Decide whether the family-definition sweep is still wanted. If the pairwise pathway stays the reference (the evidence so far favours it), the sweep protects only the `mmseqs` pathway. #137 may settle that first.
3. If wanted: run the 16-point grid on pezizo_set1 and sordariales_shallow, with controls and BUSCO scored in one table for every point, then run `collate_sweep.py`. At 2408 and 3305 CPU-hours for one Tier C+H run, 16 points per clade would be on the order of 40,000-100,000 CPU-hours (my extrapolation from single runs, not a measurement). Make #137 a hard gate: do not launch this until #137 says the `mmseqs` pathway is worth keeping. This needs your decision before launch.
4. Wire `--support_matrix` into the workflow only after the meaning of "concordant" is settled. Settle a third adjudicator first (TBLASTN or genomic evidence, or OrthoFinder), because the two methods differ in sensitivity and agreement alone does not show which call is right.
