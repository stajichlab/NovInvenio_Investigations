# 2026-09-25 — first runs, F. oxysporum f. sp. lycopersici (FOL)

Both runs finished (SUCCESS). Run 1 had rescue off; run 2 had rescue on.

## Setup

- Pipeline: NovInvenio af6fd68b69311f75417e06a962c906cfe736af41
  (detached worktree `NovInvenio-worktrees/fusarium-runs-af6fd68`).
  Nextflow 26.04.6. mmseqs tier-1 0.9 id / 0.8 cov.
- Params files: `params_norescue.yaml`, `params_rescue.yaml`
  (rescue flag, `pangenome_pfam_chunk_size: 2000`,
  `pangenome_mash_threshold: 0.0005`). Island Pfam HMM set.
- 25 strains: 21 FOL (IN), 4 FOSC (OUT). The pipeline dedup kept all 21
  ingroup strains as representatives.

## Strains and QC

28 BFD strains are labelled exactly "f. sp. lycopersici". All 28 pass
BUSCO genome >= 90% (98.6-99.6%) and the 11,000-protein floor
(13,800-16,388 proteins; FOXY median 14,956, MAD 328). The 3-MAD floor
(13,972) was not applied; it would drop Fol59 (13,800, N50 3.0 Mb).

Clone dedup at mash < 0.0005 dropped 7 (kept strain in brackets):
- Fol074 (Fol038), mash 0.000498.
- VPRI44099, VPRI44107, VPRI44272, VPRI44275, VPRI44276 (VPRI44100).
  VPRI44100 = VPRI44276 and VPRI44107 = VPRI44272 at mash 0, same N50 and
  length: duplicate assemblies.
- VPRI44273 (VPRI44274), mash 0.000315.

Ingroup (21): 4287, Fol002, Fol004, Fol007, Fol014, Fol016, Fol018, Fol026,
Fol029, Fol038, Fol069, Fol072, Fol073, Fol075, Fol59, MN25, race 3,
RR23I-1.2, VPRI44100, VPRI44274, VPRI44277.

Outgroup (4): GL1080 (f. sp. fragariae), AJ275 (f. sp. narcissi), LA1E
(f. sp. vasinfectum), *F. odoratissimum* 36102. Rule in `../README.md`.
Minimum mash to any FOL strain: 0.0107, 0.0110, 0.0113, 0.0243.

## Mash structure of the 28 FOL strains (BFD sketches)

- Maximum pairwise distance: 0.0252.
- Single-linkage groups: 21 at < 0.0005, 7 at < 0.002, 5 at < 0.005,
  3 at < 0.01 and < 0.015.
- At < 0.005: one group of 22 (4287, race 3, Fol029, Fol59, Fol002-038,
  Fol073, Fol074, all VPRI); {Fol016, Fol072, Fol075}; Fol069; MN25;
  RR23I-1.2.
- At < 0.01, Fol069 and MN25 join {Fol016, Fol072, Fol075}. RR23I-1.2 stays
  alone: 0.024 from every other FOL strain, about as far as
  *F. odoratissimum* (0.0243).
- The pipeline's own mash clades on the 21 ingroup strains give the same
  5 groups: clade_0 = 15, clade_1 = {Fol016, Fol072, Fol075}, and
  RR23I-1.2, Fol069, MN25 alone.
- So f. sp. lycopersici here is at least 3 lineages (one clonal-rich
  lineage of 22, one of 5, and RR23I-1.2). This fits a polyphyletic forma specialis.

## Pre-launch checks

- Manual mmseqs with the same parameters: 42,465 families; 6,601 eligible
  co-occurrence families, 21.8 M pairs.
- Rescue estimate: 72.5 cpu-h (model 2.0e-6 cpu-h per query x genome Mb).

## Runtime and compute

| | run 1 (no rescue) | run 2 (rescue) |
|---|---|---|
| head job | 29099753 | 29099764 |
| wall time | 0 h 39 min | 2 h 01 min |
| used cpu-h (trace, realtime x %cpu) | 4.6 | 72.8 |
| TBLASTN_PER_STRAIN used cpu-h | — | 67.1 (25 tasks, 32-55 min) |
| FAMILY_PFAM_SCAN used cpu-h | 4.0 (8 chunks, 9-31 min) | 5.2 (9 chunks) |
| allocated cpu-h (sacct) | 27.6 | 107.4 |
| failed tasks | 0 | 0 |

The rescue estimate (72.5) against the measured TBLASTN (67.1): 8% high.
No Pfam chunk came near the 4 h limit.

## Families and bins

| | run 1 | run 2 |
|---|---|---|
| families | 43,364 | 43,367 |
| core (>= 0.95 of 21) | 6,171 | 9,354 |
| soft_core | 3,872 | 2,382 |
| shell | 7,833 | 8,715 |
| cloud | 6,927 | 7,691 |
| singleton, 1 ingroup strain | 10,327 | 9,044 |
| outgroup only (0 ingroup) | 8,234 | 6,181 |
| in all 25 strains | 4,141 | 6,302 |
| families per ingroup strain, median (range) | 14,929 (13,652-15,437) | 16,721 (16,159-19,820) |

Rescue: 7,791,855 tblastn rows; 582,724 passed the thresholds; 450,830
rejected (overlap with another family's gene), 33,537 (short rep), 28,504
(hotspot); 50,737 cells changed absent -> present.

Openness fit: gamma 0.336 (run 1), 0.293 (run 2); both is_open = True.

## Co-occurrence, pairs, islands, modules

| | run 1 | run 2 |
|---|---|---|
| eligible families / candidate pairs | 6,865 / 23.6 M | 7,607 / 28.9 M |
| passed prefilter | 7,238,601 | 6,941,662 |
| FDR-significant pairs | 820,493 | 749,978 |
| trans_unconfirmed | 765,438 | 727,048 |
| trans | 49,880 | 17,778 |
| unexplained_physical | 2,896 | 2,809 |
| ambiguous_linkage | 2,279 | 2,343 |
| significant islands | 2,865 | 3,205 |
| Leiden modules (r = 1.0) | 132 (largest 698, 590, 341) | 267 (largest 358, 306, 266) |

## Assembly quality vs content (`assembly_quality_correlations.tsv`, n = 21)

| comparison | run 1 partial rho | run 2 partial rho |
|---|---|---|
| accessory_present vs N50 | -0.269 | -0.159 |
| accessory_present vs n_contigs | 0.242 | 0.313 |
| core_missing vs N50 | 0.705 | 0.708 |
| core_missing vs n_contigs | -0.557 | -0.573 |

- Coccidioides accessory vs N50: -0.436 / -0.374 (reciprocal runs), -0.519
  (529-strain baseline). FOL's -0.269 / -0.159 are weaker, with n = 21.
- core_missing rises with N50 (+0.71). The direction is the opposite of a
  fragmentation artifact. The 5 strains with N50 > 2 Mb (4287, race 3,
  Fol029, Fol59, MN25) include the more divergent ones (MN25). I did not test
  the cause.

## Top island Pfam enrichments (FDR < 0.05)

Run 1: 27 domains. Top by q: Epimerase (67/4,526 island families vs
128/14,760 background, q 1.8e-4), bZIP_1 (24 vs 32, q 1.8e-4), adh_short
(99 vs 215, q 5.3e-4), adh_short_C2, NAD_binding_10, ADH_N, HLH,
Methyltransf_23, PDDEXK_12, NmrA, Methyltransf_11, Glyco_hydro_43, NACHT.

Run 2: 33 domains. Top by q: OrsD (31/5,386 vs 37/16,406, q 2.8e-7),
NPHP3_N (105 vs 201, q 5.4e-6), bZIP_1 (q 5.4e-6), NACHT, WHD_GPIID,
Ank_2, Ank, Helitron_like_N, Glyco_hydro_43, ADH_zinc_N, Ank_KRIT1,
Methyltransf_23.

## No rescue vs rescue

- Rescue moved 3,183 families into core (6,171 -> 9,354) and cut
  outgroup-only families by 2,053.
- Trans pairs fell 64% (49,880 -> 17,778). Total FDR pairs fell 9%.
- Leiden module count doubled (132 -> 267), and the largest module shrank
  (698 -> 358).
- The accessory-vs-N50 partial rho weakened (-0.269 -> -0.159).
- The island Pfam enrichment kept bZIP_1, NACHT, Glyco_hydro_43 and
  Methyltransf_23. It gained NACHT/ankyrin/WHD domains (NPHP3_N, WHD_GPIID,
  Ank, Ank_2) and OrsD.

## Files

`results/mmseqs_norescue/output/pangenome/`, `results/mmseqs_rescue/output/pangenome/`
(report in `report/`, `island_synteny.html`). Traces in `results/*/trace.txt`.
