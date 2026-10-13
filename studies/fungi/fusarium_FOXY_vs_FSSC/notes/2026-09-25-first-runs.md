# 2026-09-25 — first reciprocal runs, FOXY vs FSSC

Status: both runs stopped at COOCCURRENCE. No islands, trans pairs or Leiden
modules were produced. The numbers below are the steps that finished.

## Setup

- Pipeline: NovInvenio commit af6fd68b69311f75417e06a962c906cfe736af41
  (detached worktree `NovInvenio-worktrees/fusarium-runs-af6fd68`).
- Nextflow 26.04.6. Backend mmseqs, tier-1 0.9 id / 0.8 cov (defaults).
- Rescue disabled (`params_norescue.yaml`). Island Pfam HMM set.
- Strains: 328 total. 209 FOXY, 119 FSSC. Selection rule: `../README.md`.
- Pipeline dedup (mash < 0.001) kept 327 of 328. It dropped virg_Mont-1.
- Mash clades (ingroup only): FOXY-in 5 clades (47/19/44/71/28).
  FSSC-in 8 clades (16/8/25/11/35/6/7/11).

## Jobs

| attempt | FOXY_in_FSSC_out | FSSC_in_FOXY_out | result |
|---|---|---|---|
| 1 | 29084981 | 29084982 | cancelled: rescue ran (CLI boolean bug, see README) |
| 2 | 29086546 | 29086547 | cancelled at 6 h 05 min: COOCCURRENCE cannot finish |

## Compute

- Head-job wall time: 0 h 48 min (attempt 1) + 6 h 05 min (attempt 2).
- Used cpu-h, completed tasks, from `trace.txt` (realtime x %cpu):
  35.6 (FOXY-in) and 49.7 (FSSC-in). FAMILY_PFAM_SCAN is 29.2 and 44.1 of it.
- Allocated cpu-h, all our SLURM jobs, all attempts (sacct CPUTimeRAW): 856.4.
  FAMILY_PFAM_SCAN is 793 of that (8 cpus per task, 2 h first-attempt limit;
  17 tasks hit the limit, exit 140, and reran on the 8 h retry).
- Full rescue, estimated from the measured matrix: 13,347 cpu-h per direction.

## Families

Tier-1 families: 399,749 (FOXY-in run) and 399,851 (FSSC-in run). Same
input proteins; mmseqs does not give the identical count on the two runs.
Coccidioides: 47,745 families for 529 strains.

Frequency bins (fraction of ingroup representative strains). The pipeline
puts families with 0 ingroup carriers in the `singleton` bin; they are split
out here.

| bin | FOXY-in | FSSC-in |
|---|---|---|
| core (>= 0.95) | 4,600 | 879 |
| soft_core | 3,605 | 1,516 |
| shell | 9,838 | 15,633 |
| cloud | 66,664 | 91,726 |
| singleton, 1 ingroup strain | 76,212 | 133,263 |
| 0 ingroup strains (outgroup only) | 238,830 | 156,834 |

Families per strain (median): FOXY 14,876 (12,898-16,531).
FSSC 14,713 (414-18,401).

## How the divergence shows up

All numbers from the unrescued `presence_matrix.tsv` of each run. The
Coccidioides values are from the unrescued matrix of its reciprocal runs
(NovInvenio 91e3157), so rescue does not differ between the two studies.

| measure | FOXY-in | FSSC-in | Cocci immitis-in | Cocci posadasii-in |
|---|---|---|---|---|
| families with >= 1 ingroup carrier | 160,919 | 243,512 | 25,729 | 34,795 |
| ingroup core (>= 95%) | 4,600 | 879 | 6,357 | 5,914 |
| core / ingroup-present families | 2.9% | 0.36% | 24.7% | 17.0% |
| families in both groups | 4,572 (1.1%) | 4,574 (1.1%) | 12,783 (26.8%) | 12,787 (26.8%) |
| families in >= 95% of both groups | 253 | 254 | 5,178 | 5,185 |
| families in all strains | 3 | 3 | n/a | n/a |

- At 90% identity, only 1.1% of families contain both FOXY and FSSC proteins.
  In Coccidioides the value is 26.8%.
- 253 families are in >= 95% of both complexes. A single FOXY strain has
  about 14,900 families. So the cross-complex "core" at this clustering
  identity is under 2% of one genome's gene set.
- FSSC core is 879 families. FSSC here spans 42 names (about 40 species,
  mash 0-0.133 to *F. vanettenii*). At 90% identity most orthologs
  from different FSSC species fall into different families. This is a
  clustering-identity effect, not only biology. I did not test other identities.
- FOXY core is 4,600 families. FOXY is one species complex (mash 0-0.029
  to Fol4287).

## Assembly quality vs content (`assembly_quality_correlations.tsv`)

accessory_present vs N50, partial rho controlling assembly length:

| run | n | rho_partial |
|---|---|---|
| FOXY-in | 209 | -0.295 |
| FSSC-in | 119 | +0.132 |
| Cocci immitis-in (rescued) | 169 | -0.436 |
| Cocci posadasii-in (rescued) | 360 | -0.374 |
| Cocci 529-strain tier-1 sweep baseline (#132) | 529 | -0.519 |

- FOXY-in rho is weaker than all Coccidioides values.
- FSSC-in rho is positive and small. The FSSC ingroup mixes about 40 species,
  so accessory content per strain likely tracks species identity more than
  assembly quality. I did not test this.
- Coccidioides rows are rescued matrices; Fusarium rows are not rescued.

## Why COOCCURRENCE did not finish

`pangenome_cooccurrence.py` tests every pair of shell+cloud families with
>= 5 ingroup carriers.

| | eligible families | candidate pairs |
|---|---|---|
| Cocci immitis-in | 7,081 | 25.1 M |
| FOXY-in | 39,805 | 792.2 M |
| FSSC-in | 45,370 | 1,029.2 M |

- FOXY-in attempt 2 (64 GB) screened 752.6 M of 792.2 M pairs in 3 h 51 min.
  121.6 M pairs had passed the prefilter (screen_alpha 0.2), about 16%.
  It was then killed for memory (MaxRSS 67.1 GB).
- FSSC-in: 51.7 M prefilter survivors after 205.8 M pairs (25%) at 2 h 50 min,
  31 GB. By rate it would reach about 260 M survivors.
- After BH, every FDR survivor gets a Python per-pair exact stratified test.
  I timed that test at 0.58 ms per pair (n = 209, 5 clades). 10 M pairs is
  about 1.6 h; 100 M pairs is about 16 h, plus one result dict per pair.
- So more memory alone does not make this finish in a sensible time. The
  retry ladder (64 -> 96 GB) would have failed. I stopped both runs.

Eligible families by carrier floor (`pangenome_min_strain_count`):

| floor | FOXY-in eligible / pairs | FSSC-in eligible / pairs |
|---|---|---|
| 5 (default) | 39,805 / 792 M | 45,370 / 1,029 M |
| 10 | 24,945 / 311 M | 24,296 / 295 M |
| 20 | 13,513 / 91 M | 14,505 / 105 M |
| 40 | 8,439 / 36 M | 8,612 / 37 M |

A floor of 40 gives about the Coccidioides pair count. It changes what is
tested, so it is a user decision.

## Data problems found

- tuar_NRRL_46518 (GCA_013363205.1) has 416 proteins in its BFD funannotate
  output (genome BUSCO 95.7%). It is a failed annotation. The selection QC
  used genome BUSCO only and did not catch it.
- 20 of 119 FSSC strains have < 75% of the FSSC median protein count
  (14,730): 7,233-10,177 proteins. 4 are EtdFoc/EthFoc isolates, 5 are UCR
  isolates. I did not check whether this is biology or annotation.
  FOXY has none (12,900-17,367). Table: `selection/protein_counts.tsv`.
- F. tuaranense NRRL 46518 is in BFD twice (two assemblies). The two are
  at mash >= 0.001, so dedup kept both.

## Options for the next run (not done)

1. Drop failed annotations: add a protein-count floor to `bin/select_strains.py`.
2. Make COOCCURRENCE feasible: raise `pangenome_min_strain_count`, or split
   the study so FSSC is not ~40 species in one ingroup.
3. Pfam scan: 2 h first-attempt limit is too short at this scale
   (17 of 62 chunks timed out). A study-config override to `time = 8.h`
   for FAMILY_PFAM_SCAN avoids the wasted first attempt.
4. Rescue (13,347 cpu-h per direction) needs a user decision.
5. Issue #132 caveat (b): these runs are the 0.9/0.8 baseline only. The
   cov-mode-1 paralog test has not been run on Fusarium.
