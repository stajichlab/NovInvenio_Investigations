# Why the first run found 0 co-occurring pairs (2026-09-25)

Read-only analysis of `results/mmseqs_norescue/output/pangenome/`. Script:
`bin/diag_cooccurrence_power.py` (exact one-sided Fisher p for every pair of
families, via a lookup on the 2x2 counts; the same test as
`bin/pangenome_cooccurrence.py`). Family filter: 5 <= carriers < 24, bin
shell/cloud/soft_core. This gives 12,163 families and 74.0 M pairs, a little more
than the pipeline's 56.6 M, so the counts below differ slightly from the run.

## Result: a near miss, not a structural limit

- 24 representative ingroup strains from 16 species (1-3 strains per species).
- Smallest possible pair p = 3.7e-7, and pairs reach it.
- BH over all pairs: a block of N pairs with p <= t passes when N >= t * m / 0.05.

| p <= | pairs (all 74.0 M) | N needed | pairs (carriers 8-16, 8.4 M) | N needed |
|---|---|---|---|---|
| 1e-4 | 45,086 | 147,926 | 24,496 | 16,757 |
| 1e-5 | 6,881 | 14,793 | 4,828 | 1,676 |
| 1e-6 | 573 | 1,479 | 573 | 168 |

- Over all pairs, the strong tail is 2-3x too small at every threshold.
- Restricting to families carried by 8-16 of the 24 strains (4,094 families) keeps
  the same strongest pairs but cuts m by 9x, and all thresholds pass. Carriers
  10-14 (2,166 families) also passes.
- Most eligible families span many species (median about 9 of 16), so a pair
  rarely co-occurs perfectly across 24 strains.

## Caveat before using a carrier window

With 16 species and 1-3 strains each, pairs that pass are likely dominated by
shared ancestry (families gained or lost together on one branch). The pipeline
assigned only 2 mash clades here, so its within-clade stratified test controls
little of that. A tree-aware analysis (`--pangenome_species_tree`, #185) is the
better fit for this set.
