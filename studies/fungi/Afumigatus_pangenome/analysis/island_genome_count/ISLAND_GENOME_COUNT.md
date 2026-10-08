# Island genome count: how many genomes give island predictions?

Started 2026-09-29. Question from the PI: what is the smallest, fastest test
set that still yields accessory-island predictions, so island drawing can be
tested without an 18-21 h full run? And do known Starships show up as islands?

## Why

- `Afumigatus_test10/in10_v1` (9 representative genomes, NovInvenio f87fd1e)
  completed but COOCCURRENCE found 0 FDR-significant pairs (636,681 passed the
  prefilter), so 0 islands. The island, locus and clinker steps ran only on
  empty input.
- In `full_v070` (293 strains, 121 representatives) the island steps are cheap:
  COOCCURRENCE 16 min 53 s, PAIR_CLASSIFICATION 5 min 51 s, BUILD_ISLANDS 20 s.
  The 18 h is upstream (clustering, rescue).

## Design

- Parent data: `results/full_v070/config/pangenome/` (presence_matrix.rescued,
  family_positions, gene_positions, tier1_cluster, strain_inventory). This run
  completed through BUILD_ISLANDS (11,775 islands) before it failed at
  ISLAND_LOCI (NovInvenio #213).
- Sampling unit: representative ingroup strains (`is_representative == 1`),
  because COOCCURRENCE uses only representatives. N = 10, 15, 20, 30, 45, 60,
  80, 100 with seeds 1-3, plus N = 121 (all). Af293 is kept in every subset;
  both outgroup references are kept (as in the full run).
- Steps rerun per subset, with the commands and parameters copied from the
  pipeline's `.command.sh` (NovInvenio f87fd1e): FREQUENCY_BINS, COOCCURRENCE,
  PAIR_CLASSIFICATION, BUILD_ISLANDS. Environment: NovInvenio pixi env.
- Starship truth (Gluck-Thaler et al. 2025, mBio, doi:10.1128/mbio.01092-25;
  `mbio.01092-25-s0002.xlsx`, in the study DATA_MANIFEST):
  - Table S7: 8 Af293 elements (7 high confidence) ->
    `inputs/af293_starships_table_s7.tsv` via `scripts/extract_af293_starships.py`.
    RefSeq NC_007194.1-NC_007201.1 map to GenBank CM000169.1-CM000176.1; all 8
    lengths match exactly (NCBI esummary, 2026-09-29).
  - Table S21 (20 high-confidence Starships, presence per strain) ->
    `results/id_crosswalk/ground_truth_starships_by_short.tsv` (earlier session).
- Scoring (`scripts/score_starships.py`): see its docstring. Coordinate check
  in Af293; presence-pattern best Jaccard vs a 200-permutation null (95th pct).

## Limits

- Subsets reuse families clustered from all 293 genomes. A real N-genome run
  clusters its own families, so family content and rescue calls would differ.
- Only Af293 has Starship coordinates in our assembly naming; our A1163
  assembly (DS499594.1...) is not the one used in the paper.
- Rescued (genome_only) cells have no gene coordinates and are not placed in Af293.

## Reference: full_v070 (293 strains, pipeline a394ace)

Scored with `score_starships.py` against all 261 strains with a genotype:
- 11,775 islands; 7,192 placed in Af293; 985 of those hit an Af293 Starship.
- All 7 high-confidence Af293 Starships hit by >= 1 island (21 to 307 islands each).
- 21.3% of placed island genes lie in a Starship, vs 5.2% of non-core Af293 genes.
- Presence pattern: 19 of 19 scorable Starships above the permutation null
  (e.g. Gnosis-h1 best Jaccard 0.98 vs null 0.60; Tardis-h1 0.97 vs 0.57).

## Results (SLURM 29207128-31, all COMPLETED 2026-09-29, 25 runs)

Per run: `outputs/sweep_summary.tsv`; per N: `outputs/sweep_by_n.tsv`.

| N reps | co-occurring pairs (3 seeds) | islands | Af293 high-conf Starships hit (of 7) | presence matches > null / scored | island steps (s) |
|---|---|---|---|---|---|
| 10 | 0, 0, 0 | 0 | 0 | 0 | 35-66 |
| 15 | 0, 0, 0 | 0 | 0 | 0 | 91-100 |
| 20 | 0, 0, 0 | 0 | 0 | 0 | 117-143 |
| 30 | 154,865 / 2,273 / 5,246 | 2,306 / 1,549 / 1,546 | 7 / 6 / 7 | 4/5, 7/8, 5/7 | 170-316 |
| 45 | 115,209 / 269,098 / 75,194 | 3,152-3,570 | 7 all | 8/13, 7/11, 10/13 | 324-471 |
| 60 | 47,191-530,806 | 3,682-4,792 | 7 all | 10/12, 12/14, 12/15 | 354-616 |
| 80 | 195,399-405,503 | 4,989-5,701 | 7 all | 14/15, 14/15, 12/13 | 498-639 |
| 100 | 214,515-999,062 | 5,957-7,045 | 7 all | 15/15, 15/16, 15/16 | 591-922 |
| 121 | 1,363,573 | 8,226 | 7 | 16/17 | 1,418 |

- The threshold lies between 20 and 30 representatives: every run at <= 20 had 0
  FDR-significant pairs; every run at >= 30 had > 1,500 islands.
- Island genes in Starships: 14.8-20.6% at every N >= 30, vs a non-core
  background of 4.4-5.7%.
- Explained (2026-10-08, `outputs/gap_check/`): N=121 here gave 8,226 islands vs 11,775 in
  full_v070 because `scripts/make_subset.py` keeps only the kept strains' rows of
  `family_positions.tsv.zst`. full_v070's PAIR_CLASSIFICATION saw gene positions for all 295
  strains (123 representatives + 172 near-duplicates and outgroups). With the same
  f87fd1e scripts and the same full_v070 co-occurrence output (1,363,573 pairs):
  all 295 strains' positions give 11,775 islands (exactly full_v070); the 123 sweep strains'
  positions give 8,226 (exactly this sweep). `insufficient_data` pairs: 71,662 vs 188,361;
  `trans`: 564,245 vs 543,704. So script version (a394ace vs f87fd1e) and #212 are not the cause.
  Which step causes the drop (variant C, job 29630183, same scripts): pair classification
  with the 123 representatives' positions but islands built from all 295 strains' positions
  gives 11,676 islands, 99 fewer (0.8%) than 11,775. So the classification evidence
  (`n_co_carrying`, `linkage_fraction`, which use every strain in family_positions) accounts for
  about 100 islands. The other 3,450 (11,676 to 8,226) come from BUILDING islands with only the
  representatives' positions: they are island entries found only in the 172 non-representative
  strains. Whether those are islands that moved or were gained in near-identical strains, or
  one-gene differences that split one island into several entries (islands are deduplicated on
  the exact member set), was not tested. A real run on N dereplicated genomes would resemble
  the sweep, not full_v070.

## Reproduce

`bash run.sh`, then `python scripts/summarize.py` after the jobs finish.
