# #132 priority 2: core / soft-core / shell cutoffs (0.95 / 0.90 / 0.15)

2026-09-28. Script: `p2_frequency_cutoffs.py`. Outputs: `p2_summary.tsv`, `p2_sensitivity.tsv`,
`p2_mixture.tsv`, `p2_spectra.png`. Runtime 17 s, no jobs.

## Question

Do the cutoffs sit on natural breaks in the family frequency spectrum, or is the spectrum a
continuum where any cutoff is a convention?

## Data

The pipeline's own `frequency_table.tsv` (per-group bins after #212 where available). Families
with `strain_count >= 1` in the ingroup representatives; `n` = ingroup representatives after
dereplication. Old-format tables put `strain_count = 0` in `singleton`; those are excluded here.

| dataset | table | n | families |
|---|---|---|---|
| Cocci_Ci | `coccidioides_pangenome/results/rescue_freqpol_immitis_in_posadasii_out` | 133 | 26,051 |
| Cocci_Cp | `coccidioides_pangenome/results/rescue_freqpol_posadasii_in_immitis_out` | 278 | 32,576 |
| Afumigatus | `Afumigatus_pangenome/results/full_v070/config` (v0.7.0; tables before the ISLAND_LOCI failure) | 121 | 36,302 |
| FOXY | `fusarium_FOXY_vs_FSSC/results/mmseqs_norescue_FOXY_in_FSSC_out` | 209 | 160,919 |
| FSSC_ambrosia | `FSSC_ambrosia/results/mmseqs_norescue` | 24 | 61,791 |
| FOL | `fusarium_FOL/results/mmseqs_rescue` | 21 | 37,186 |
| FSSC_solani | `FSSC_solani/results/mmseqs_norescue` | 13 | 34,286 |

All paths under `studies/fungi/`, `output/pangenome/frequency_table.tsv` unless noted. Class
counts computed here match the pipeline's `bin` column (checked on Cocci_Ci: 7,808 / 6,786 /
4,012 / 609 / 6,836).

## Results

**1. No trough at any cutoff.** Every spectrum is U-shaped (figure). The minimum of the
smoothed spectrum lies at frequency 0.48-0.66, in the middle of the shell class. Over
0.2-0.8 the smoothed density varies only 2.3x-4.7x. 0.95 and 0.90 sit on the rising flank
of the core peak; 0.15 sits on the smooth decay of the cloud arm.

**2. Class sizes move steadily with the cutoff** (`p2_sensitivity.tsv`; change against the
default):

| dataset | core at 0.90 | 0.97 | 0.99 | 1.00 | shell at 0.10 | 0.20 |
|---|---|---|---|---|---|---|
| Cocci_Ci | +9% | -10% | -34% | -64% | +17% | -13% |
| Cocci_Cp | +10% | -6% | -39% | -86% | +24% | -14% |
| Afumigatus | +9% | -11% | -29% | -52% | +19% | -12% |
| FOXY | +78% | -20% | -64% | -94% | +33% | -17% |

So a strict 100% core is not usable: 52-94% of the default core misses at least one
representative (annotation or assembly gaps).

**3. Binomial mixture** (zero-truncated, Snipen et al. 2009; EM, K = 2-10, 5 starts each;
`p2_mixture.tsv`).
- Large panels (n = 121-278): BIC keeps falling up to K = 10. The mixture finds no fixed
  number of classes, which is a continuum result.
- The mixture's core boundary (smallest k whose most likely component is the top one) rises
  with K. For K = 6-10 it lies at 0.947-0.977 (Ci), 0.942-0.975 (Cp), 0.959-0.983
  (Afumigatus) and 0.871-0.962 (FOXY). 0.95 is inside or next to these ranges.
- At K = 3 (a core / shell / cloud picture) the low-component boundary is 0.150 (Ci),
  0.133 (Cp), 0.149 (Afumigatus), 0.077 (FOXY). So 0.15 matches a three-class fit for
  Cocci and Afumigatus. But BIC rejects K = 3 by a wide margin on every large panel, so
  this is not evidence for 0.15.
- Small panels (n = 13-24): BIC picks K = 5-6. The core boundary is k = n (FOL, FSSC_solani)
  or k = 22 of 24 (FSSC_ambrosia).

**4. Small n breaks the class definitions.** The cutoffs map to integer counts
(`k_core>=`, `k_soft>=`, `k_shell>=` in `p2_summary.tsv`):
- n = 13: core = 13/13, soft-core = 12/13, shell = k >= 2. The **cloud class is empty**:
  ceil(0.15 x 13) = 2. Cloud exists only for n >= 14.
- n = 21 and 24: core = n-1 or more, soft-core = exactly n-2. "Soft core" is one count.

## Recommendation and evidence class

| cutoff | recommendation | evidence class |
|---|---|---|
| core 0.95 | keep | empirical: within the mixture core boundary range (0.94-0.98) on the four large panels; a strict 1.0 loses 52-94% of core |
| soft-core 0.90 | keep, report as convention | convention: no break; at n < 25 it is one strain count |
| shell 0.15 | keep, report as convention; report shell/cloud counts with their sensitivity (about +20% / -13% for 0.10 / 0.20) | convention: the spectrum is a continuum there |

Pipeline changes this suggests (not made):
- Warn in the report when n < 14 (cloud cannot exist) and when n < 25 (soft-core is one count).
- For small panels, report the spectrum itself (families per k), which carries all the
  information the classes do.

## Limits

- The mixture treats families as independent draws; lineage structure (the SNP-tree result
  in `studies/fungi/coccidioides_pangenome/analysis/sans/`) makes presence correlated.
  The mixture is a description of the spectrum, not a generative model of it.
- One table per taxon; the tier-1 clustering (0.9 / 0.8) is fixed across all of them.
