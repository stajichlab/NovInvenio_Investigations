# Issue #132 tier-1 sweep: downstream metrics for all 8 grid points

Addendum to `2026-09-20-tier1-clustering-sweep.md`. That note ran islands, trans
pairs and Leiden modules for 3 settings only. This note fills the gap: all 8
`min_id x cov` grid points now have the full downstream metric set.

## What ran

- Dataset: `config_genus_vs_ureesii.csv` (530 proteomes, 4,512,100 proteins), the
  same input as the 2026-09-20 sweep.
- Pipeline: `pangenome.nf` at NovInvenio commit `af6fd68b69311f75417e06a962c906cfe736af41`.
  The runs used a `git archive` snapshot of that commit
  (`studies/fungi/coccidioides_pangenome/.nf_launch/tier1_sweep_pipeline_af6fd68/`).
  Live commits in the shared checkout could not change the runs.
- Launcher: `studies/fungi/coccidioides_pangenome/run_tier1_sweep.sh <min_id> <cov>`.
  One launch dir per point (`.nf_launch/tier1_sweep_id<min_id>_cov<cov>/`).
  Outputs are in `results/tier1_sweep/coccidioides_sweep_id<min_id>_cov<cov>/pangenome/`.
- Every point, including the 0.9/0.8 baseline, ran with rescue off.
  So all 8 points are like-for-like.
- `cov-mode` was 0 for every point (pipeline default).
- Collector: `studies/fungi/coccidioides_pangenome/bin/collect_tier1_sweep.py`.
  Output: `results/tier1_sweep/tier1_sweep_metrics.tsv`.
  Test: `bin/test_collect_tier1_sweep.py` (5 tests pass).

## What was reused, and what was not

- Nothing was reused from the 2026-09-20 sweep. All 8 points were re-run
  pipeline-native, including clustering.
- Reason for re-clustering: `CLUSTER_TIER1` has no `storeDir`. Injecting the old
  cluster TSVs would need pipeline edits. Clustering took 4-6 min per point
  (32 CPUs), so re-running was cheaper.
- Reason for re-running 0.70/0.50: the 2026-09-20 tier-2 numbers came from a
  standalone script chain, not the pipeline. That chain used Leiden seed 42 (pipeline:
  0) and the scripts' own defaults for pair classification. It also did not produce
  `assembly_quality_correlations.tsv`. Its numbers are therefore not comparable with
  these runs. Example: the standalone baseline gave 966,865 trans pairs; the
  pipeline baseline here gives 543,867.
- Re-clustering does not give exactly the same families as the 2026-09-20 runs.
  Baseline 0.9/0.8: 54,412 families here vs 54,421 then. 0.70/0.50: 30,435 vs
  30,406. I did not investigate the cause.

## Deviations from a default pipeline run

1. Rescue off. It is set as `params.pangenome_rescue_enable = false` in
   `.nf_launch/tier1_sweep_resources.config`, not on the command line.
   Nextflow 26.04 keeps a CLI `--pangenome_rescue_enable false` as the String
   `"false"`. Groovy treats that String as true. The first launch did run the rescue
   branch. I cancelled it before any `TBLASTN_PER_STRAIN` task started. This affects
   every boolean CLI flag of the pipeline (see "Pipeline issues found" below).
2. `--pangenome_island_pfam_hmm` points at a one-profile stub HMM
   (`results/tier1_sweep/_inputs/single_profile_stub.hmm`, PF26733.1).
   `BUILD_ISLANDS` only runs when this parameter is set. Islands do not use Pfam
   output. The full Pfam-A scan costs about 80 CPU-h per point at this scale.
   Consequence: `island_pfam_enrichment.tsv`, `module_domains.tsv` and
   `pfam.domtblout` in these runs have no meaning. Do not use them.
3. `COOCCURRENCE` and `PAIR_CLASSIFICATION` got a 16 h / 12 h first-attempt limit.
   Six of eight `COOCCURRENCE` jobs ran on the `stajichlab` partition (c01, AMD
   Opteron 6376). The preempt account was at its CPU cap. Two ran on preempt.
   This changes only wall time, not results.
4. Both `--pangenome_project` and `--project` are passed. With only
   `--pangenome_project`, the outputs go to `<outdir>/output/`. All 8 points would
   have overwritten each other.

## Results (all 8 points, rescue off, cov-mode 0)

Baseline for AMI/ARI = 0.9/0.8.

| min_id | cov | families | core | soft_core | shell | cloud | singleton | islands | trans pairs | Leiden modules | largest module |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.70 | 0.5 | 30,435 | 6,548 | 264 | 3,444 | 9,541 | 10,638 | 13,887 | 279,728 | 6 | 4,404 |
| 0.70 | 0.8 | 46,518 | 5,476 | 462 | 4,931 | 15,821 | 19,828 | 27,005 | 565,298 | 6 | 6,818 |
| 0.80 | 0.5 | 33,051 | 6,517 | 274 | 3,467 | 9,978 | 12,815 | 13,850 | 285,110 | 7 | 4,418 |
| 0.80 | 0.8 | 48,754 | 5,502 | 470 | 4,875 | 15,995 | 21,912 | 26,454 | 548,136 | 6 | 6,694 |
| 0.90 | 0.5 | 38,812 | 6,364 | 271 | 3,770 | 11,293 | 17,114 | 16,107 | 299,733 | 7 | 4,996 |
| **0.90** | **0.8** | **54,412** | 5,439 | 457 | 5,003 | 17,021 | 26,492 | **27,729** | **543,867** | 7 | 7,045 |
| 0.95 | 0.5 | 44,666 | 5,991 | 300 | 4,448 | 13,266 | 20,661 | 20,950 | 332,949 | 5 | 5,836 |
| 0.95 | 0.8 | 60,079 | 5,214 | 461 | 5,468 | 18,622 | 30,314 | 32,270 | 560,707 | 7 | 7,563 |

| min_id | cov | AMI (all proteins) | ARI (all) | AMI (module proteins only) | ARI (module only) | n proteins in a module in both runs | partial rho accessory~N50 | partial rho accessory~contigs |
|---|---|---|---|---|---|---|---|---|
| 0.70 | 0.5 | 0.394 | 0.572 | 0.480 | 0.627 | 830,606 | -0.570 | +0.466 |
| 0.70 | 0.8 | 0.704 | 0.866 | 0.667 | 0.739 | 1,274,737 | -0.577 | +0.493 |
| 0.80 | 0.5 | 0.415 | 0.595 | 0.480 | 0.606 | 854,924 | -0.564 | +0.462 |
| 0.80 | 0.8 | 0.713 | 0.878 | 0.661 | 0.751 | 1,273,997 | -0.575 | +0.492 |
| 0.90 | 0.5 | 0.466 | 0.648 | 0.518 | 0.638 | 929,444 | -0.563 | +0.466 |
| **0.90** | **0.8** | 1 | 1 | 1 | 1 | 1,365,947 | **-0.572** | **+0.490** |
| 0.95 | 0.5 | 0.390 | 0.596 | 0.467 | 0.601 | 967,979 | -0.562 | +0.471 |
| 0.95 | 0.8 | 0.669 | 0.828 | 0.659 | 0.741 | 1,306,499 | -0.572 | +0.492 |

Definitions:

- Islands = data rows of `significant_islands.tsv`.
- Trans pairs = rows of `pair_classification.tsv.zst` with `classification == trans`
  (permutation-confirmed). `trans_unconfirmed` counts are in the TSV.
- AMI/ARI compare protein-level module labels. Each protein is mapped to its family
  (`cluster/tier1_cluster.tsv`), then to its Leiden module (`family_modules.tsv`).
  Family IDs are never joined across runs. A protein whose family is in no module
  gets the label `none`. "All proteins" = 4,512,100 proteins, `none` is one label.
  "Module proteins only" = proteins with a real module in both runs.
- AMI uses arithmetic normalisation (scikit-learn default). The implementation was
  checked against scikit-learn 1.9.1 on random and hand cases. At n = 4.5 M the
  expected-MI correction is negligible: AMI equals NMI to 4 decimals here.
- Partial rho = the pipeline's `ASSEMBLY_QUALITY_QC` value
  (`rho_partial_length`, controls for total assembly length, n = 529 strains).
  This is not the same number as the 2026-09-20 note's rho (-0.519), which used
  the 411 dereplicated ingroup strains and a different accessory definition.
  Compare rho values within this table only.

## What the numbers show

- Coverage drives the downstream counts more than identity does.
  - cov 0.5: 13,850-20,950 islands and 279,728-332,949 trans pairs.
  - cov 0.8: 26,454-32,270 islands and 543,867-565,298 trans pairs.
- At cov 0.8, trans pairs are almost flat across min_id 0.70-0.95
  (543,867-565,298, a 4% range). Family count changes 29% over the same range
  (46,518-60,079). Islands change 22% (26,454-32,270).
- At cov 0.5, trans pairs rise with min_id (279,728 to 332,949).
- The fragmentation confound does not move. Partial rho(accessory ~ N50) spans
  -0.562 to -0.577 over all 8 points, a range of 0.015. At the same min_id, the
  cov 0.5 value is 0.007-0.011 smaller in magnitude than the cov 0.8 value. This agrees with the
  2026-09-20 tier-1 result (range 0.017 there, on its own rho definition).
- Leiden at resolution 1.0 gives 5-7 modules at every point. The largest module
  holds 4,404-7,563 families. This is the known too-coarse resolution problem
  (#132 priority 4). It is still present at every grid point.
- Module structure agreement with the baseline splits by coverage.
  - cov 0.8 points: AMI 0.67-0.71 (all proteins), ARI 0.83-0.88.
  - cov 0.5 points: AMI 0.39-0.47, ARI 0.57-0.65.
  - Changing identity at fixed cov 0.8 keeps AMI near 0.7. Changing cov to 0.5 at
    fixed identity 0.9 drops AMI to 0.47.
- No grid point gives a reason to change the 2026-09-20 recommendation. None
  reduces the confound. The downstream subject (islands, trans pairs, module
  partition) depends mostly on `cov`.

## Limits

- One taxon, one clustering tool (mmseqs2), cov-mode 0 only.
- Leiden was run once per point with seed 0. Seed-to-seed variation of the module
  partition was not measured here. Part of the AMI < 1 at non-baseline points may be
  Leiden variation, not clustering. A seed-replicate of the baseline would measure
  this. I did not run one.
- Rescue was off. With rescue on, family presence and all downstream counts change.
  These numbers describe the unrescued matrix only.

## Cost

- Allocated CPU-h per point (cpus x realtime, all Nextflow sessions in the launch dir,
  including the cancelled first launch): 8.6 (0.70/0.5), 23.1 (0.70/0.8), 8.7 (0.80/0.5),
  22.0 (0.80/0.8), 15.6 (0.90/0.5), 24.2 (0.90/0.8), 17.2 (0.95/0.5), 28.2 (0.95/0.8).
  Total 147.6 allocated CPU-h. `PAIR_CLASSIFICATION` reserves 8 CPUs, so real
  CPU use is lower than this figure.
- Wall time of the two long steps per point: `COOCCURRENCE` 1.0-7.0 h (single CPU),
  `PAIR_CLASSIFICATION` 0.5-2.0 h. The runs on the Opteron node took longer than the
  two on preempt nodes.
- Disk: `results/tier1_sweep/` holds 14 GB. The launch dirs' `work/` directories
  hold more. Neither has been cleaned.

## Pipeline issues found (not fixed here)

1. Boolean CLI flags do not work under Nextflow 26.04. `--flag false` arrives as the
   String `"false"`, which is truthy. Verified with a one-line test script. This
   affects `--pangenome_rescue_enable` (observed) and, by the same mechanism, any
   other boolean param tested with `if (params.x)` (not checked one by one). Workaround: set the value in a `-c` config file.
2. `pangenome.nf` sets `params.project = params.pangenome_project` inside the
   workflow. The publishDir closures do not see it. Earlier runs, for example
   `results/rescue_structural_genus_vs_ureesii/`, published to `<outdir>/output/`.
3. `BUILD_ISLANDS` only runs when `--pangenome_island_pfam_hmm` is set, although it
   does not use Pfam output.
4. A second `withName: '.*COOCCURRENCE'` block in a later `-c` file does not
   override `stajichlab_queue.config`'s combined `'.*TBLASTN_PER_STRAIN|.*COOCCURRENCE'`
   selector. The override must use the same selector string.

## Draft comment for issue #132 (not posted)

```markdown
**Gap (a) closed: downstream metrics for all 8 tier-1 grid points**

All 8 `min_id x cov` points (cov-mode 0) were re-run pipeline-native on the
530-proteome genus_vs_ureesii set, rescue off for every point including the
0.9/0.8 baseline. NovInvenio commit af6fd68. Nothing was reused from the
2026-09-20 standalone chain (its trans-pair counts and Leiden seed differ from the
pipeline, so it is not comparable). Full write-up:
`notes/pangenome-method-investigations/2026-09-25-tier1-sweep-downstream.md`.

| min_id | cov | families | islands | trans pairs | modules | AMI vs 0.9/0.8 | ARI | partial rho acc~N50 |
|---|---|---|---|---|---|---|---|---|
| 0.70 | 0.5 | 30,435 | 13,887 | 279,728 | 6 | 0.394 | 0.572 | -0.570 |
| 0.70 | 0.8 | 46,518 | 27,005 | 565,298 | 6 | 0.704 | 0.866 | -0.577 |
| 0.80 | 0.5 | 33,051 | 13,850 | 285,110 | 7 | 0.415 | 0.595 | -0.564 |
| 0.80 | 0.8 | 48,754 | 26,454 | 548,136 | 6 | 0.713 | 0.878 | -0.575 |
| 0.90 | 0.5 | 38,812 | 16,107 | 299,733 | 7 | 0.466 | 0.648 | -0.563 |
| 0.90 | 0.8 | 54,412 | 27,729 | 543,867 | 7 | 1 | 1 | -0.572 |
| 0.95 | 0.5 | 44,666 | 20,950 | 332,949 | 5 | 0.390 | 0.596 | -0.562 |
| 0.95 | 0.8 | 60,079 | 32,270 | 560,707 | 7 | 0.669 | 0.828 | -0.572 |

AMI/ARI are on protein-level module labels (protein -> family -> module, proteins in
no module labelled `none`); family IDs are never joined across runs. Partial rho is
the pipeline's ASSEMBLY_QUALITY_QC value (n = 529), not the 411-strain value in the
2026-09-20 note.

Findings:
- The fragmentation confound does not move: partial rho(accessory ~ N50) spans
  -0.562 to -0.577 across all 8 points.
- `cov` drives the downstream subject. At the same min_id, cov 0.5 gives 51-65% of
  the islands and 49-59% of the trans pairs of cov 0.8. At cov 0.8, trans pairs vary only 4% across min_id
  0.70-0.95.
- Module agreement with the baseline: AMI 0.67-0.71 at cov 0.8, 0.39-0.47 at cov 0.5.
  Leiden seed variation was not measured, so part of AMI < 1 may be Leiden, not
  clustering.
- Leiden r=1.0 gives 5-7 modules at every point (priority 4 still open).
- No change to the recommendation: keep 0.9/0.8, cov-mode 0.

Caveat: islands were produced with a one-profile stub Pfam HMM (BUILD_ISLANDS is
gated on `--pangenome_island_pfam_hmm`). Island counts are valid; the Pfam
enrichment outputs of these runs are not.

Side finding: under Nextflow 26.04, `--pangenome_rescue_enable false` on the command
line is the String "false" and is truthy, so rescue still runs. Set booleans in a
`-c` config file instead.
```
