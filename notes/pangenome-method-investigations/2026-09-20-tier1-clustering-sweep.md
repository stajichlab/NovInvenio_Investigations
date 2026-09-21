# Issue #132 tier-1 clustering sensitivity analysis

Dataset: Coccidioides 529-strain pangenome (`mmseqs_genus_vs_ureesii`), 530 proteomes,
4,512,100 proteins, 2.2 GB `all_strains.pep.fa`. Every number below is measured on this
dataset in this run; nothing is carried over from a prior report except where stated.

Method: 12 independent `mmseqs easy-cluster` runs over the real input, matching
`bin/pangenome_cluster_backend.py`'s invocation (`--cluster-reassign`, `--threads 32`)
and varying only the swept parameters. Frequency bins use the pipeline's own cutoffs
(0.95 / 0.90 / 0.15) over the 411 dereplicated ingroup strains, exactly as
`bin/pangenome_frequency_bins.py` does in the pipeline.

Reproduction check: the re-run at `min_id 0.9 / cov 0.8 / cov-mode 0` gives **54,421
families**, identical to the shipped pipeline result, and a median family-representative
length of **201 aa**, identical to the value reported in issue #133. The measurement
chain is therefore reproducing the pipeline, not a different computation.

## Definitions used here

- **accessory** = shell + cloud + singleton families a strain carries (the
  fragmentation artefact lands mostly in singleton, so excluding it would hide the thing
  under test).
- **rho** = **partial** Spearman correlation of a per-strain count against assembly
  contig count or N50, **controlling for total assembly length** (rank-linear residuals),
  over the 411 dereplicated ingroup strains.
- Baseline confound, measured here: **rho(accessory ~ N50) = -0.519**,
  **rho(accessory ~ contigs) = +0.435**.

## Cost actually incurred

11 full clusterings ran as one SLURM array (6 concurrent on one node), wall clock
**15 minutes total**; per-run `/usr/bin/time -v` shows 4:41-5:47 elapsed and 3.9-4.1 GB
max RSS at 32 threads. The post-cluster merge variant ran separately in under 5 minutes.
Compressed cluster TSVs + representative FASTAs for all 12 clusterings: **405 MB** on
`/rhome`. All mmseqs tmp directories were written to node-local `$SCRATCH` and deleted
in-job. Tier-2 (co-occurrence -> pair classification -> islands -> Leiden) was run for
3 settings only (1h14m, 1h29m and 2h57m wall clock). Total disk footprint of the whole
investigation: **425 MB** on `/rhome`, none on `/bigdata`.

## The mechanism is real and directly confirmed

Before sweeping anything: using the run's own `gene_positions.tsv`, a protein was scored
as **contig-terminal** if it is the first or last gene model on its contig - the position
where a gene model is most likely truncated by the assembly.

| protein population (baseline clustering) | n | fraction contig-terminal |
|---|---|---|
| all clustered proteins | 4,512,100 | **9.38%** |
| members of single-strain (private) families | 24,105 | **31.38%** |
| members of families in >=400 strains | 3,345,373 | **7.78%** |

Private-family proteins are **3.3x** enriched for sitting at a contig end relative to the
proteome as a whole, and **4.0x** relative to broadly-conserved families. The hypothesis
in #132/#133 - that contig-broken gene models fail `-c 0.8` and become singleton families -
is confirmed on this dataset.

## Tier 1 - full grid

| setting | families | core | soft_core | shell | cloud | singleton | strain-private | copies/cell | multi-copy fams |
|---|---|---|---|---|---|---|---|---|---|
| `SHIPPED_BASELINE` | 54,421 | 5,437 | 455 | 5,019 | 17,044 | 26,466 | 24,080 | 1.0006 | 1,535 |
| `id0.70_cov0.50_cm0` | 30,406 | 6,552 | 267 | 3,442 | 9,524 | 10,621 | 9,456 | 1.0027 | 1,877 |
| `id0.70_cov0.80_cm0` | 46,439 | 5,480 | 459 | 4,927 | 15,797 | 19,776 | 17,654 | 1.0021 | 1,706 |
| `id0.80_cov0.50_cm0` | 33,042 | 6,520 | 271 | 3,478 | 9,948 | 12,825 | 11,624 | 1.0013 | 1,799 |
| `id0.80_cov0.80_cm0` | 48,792 | 5,499 | 468 | 4,884 | 16,005 | 21,936 | 19,789 | 1.0009 | 1,649 |
| `id0.90_cov0.50_cm0` | 38,909 | 6,359 | 274 | 3,770 | 11,327 | 17,179 | 15,743 | 1.0009 | 1,674 |
| `id0.90_cov0.80_cm0` | 54,421 | 5,437 | 455 | 5,019 | 17,044 | 26,466 | 24,080 | 1.0006 | 1,535 |
| `id0.95_cov0.50_cm0` | 44,608 | 5,998 | 295 | 4,427 | 13,296 | 20,592 | 18,826 | 1.0006 | 1,563 |
| `id0.95_cov0.80_cm0` | 60,138 | 5,216 | 463 | 5,451 | 18,659 | 30,349 | 27,637 | 1.0005 | 1,478 |
| `id0.90_cov0.80_cm1` | 33,703 | 6,500 | 212 | 3,049 | 9,109 | 14,833 | 13,774 | 1.0370 | 4,331 |
| `id0.90_cov0.80_cm2` | 46,473 | 5,959 | 352 | 4,284 | 15,632 | 20,246 | 17,924 | 1.0007 | 1,610 |
| `id0.90_cov0.80_cm5` | 48,338 | 5,507 | 437 | 4,966 | 15,507 | 21,921 | 19,845 | 1.0011 | 1,653 |
| `postmerge_cm5_id0.90_cov0.90` | 49,528 | 5,482 | 469 | 4,908 | 15,860 | 22,809 | 20,612 | 1.0008 | 1,589 |

| setting | rho acc~contigs | rho acc~N50 | rho priv~contigs | rho priv~N50 | acc med best40 | acc med worst40 | priv med best40 | priv med worst40 | rep len median | % reps <150aa |
|---|---|---|---|---|---|---|---|---|---|---|
| `SHIPPED_BASELINE` | +0.435 | -0.519 | +0.408 | -0.401 | 2613 | 2929 | 21 | 61 | 201 | 37.2 |
| `id0.70_cov0.50_cm0` | +0.418 | -0.524 | +0.374 | -0.388 | 1666 | 1932 | 6 | 23 | 172 | 43.5 |
| `id0.70_cov0.80_cm0` | +0.440 | -0.528 | +0.399 | -0.406 | 2556 | 2876 | 17 | 53 | 190 | 39.3 |
| `id0.80_cov0.50_cm0` | +0.412 | -0.516 | +0.371 | -0.377 | 1704 | 1978 | 7 | 24 | 182 | 41.6 |
| `id0.80_cov0.80_cm0` | +0.437 | -0.525 | +0.401 | -0.403 | 2540 | 2862 | 18 | 54 | 196 | 38.1 |
| `id0.90_cov0.50_cm0` | +0.411 | -0.512 | +0.362 | -0.339 | 1860 | 2140 | 12 | 29 | 189 | 40.4 |
| `id0.90_cov0.80_cm0` | +0.435 | -0.519 | +0.408 | -0.401 | 2613 | 2929 | 21 | 61 | 201 | 37.2 |
| `id0.95_cov0.50_cm0` | +0.417 | -0.511 | +0.339 | -0.347 | 2202 | 2486 | 15 | 34 | 196 | 39.2 |
| `id0.95_cov0.80_cm0` | +0.437 | -0.519 | +0.388 | -0.395 | 2825 | 3138 | 28 | 69 | 207 | 36.1 |
| `id0.90_cov0.80_cm1` | +0.290 | -0.350 | +0.314 | -0.328 | 1492 | 1686 | 12 | 20 | 220 | 35.9 |
| `id0.90_cov0.80_cm2` | +0.427 | -0.519 | +0.399 | -0.377 | 2188 | 2482 | 16 | 44 | 188 | 40.2 |
| `id0.90_cov0.80_cm5` | +0.431 | -0.518 | +0.394 | -0.378 | 2558 | 2871 | 14 | 48 | 201 | 37.8 |
| `postmerge_cm5_id0.90_cov0.90` | +0.432 | -0.519 | +0.399 | -0.395 | 2554 | 2870 | 16 | 49 | 199 | 37.7 |

## Headline result: the fragmentation confound does NOT shrink across the min_id x cov grid

This is the measure the issue called the interesting one, and the answer is negative.

Across all eight `min_id x cov` combinations, `rho(accessory ~ N50)` moves within
**-0.511 to -0.528** - a range of 0.017, against a baseline of -0.519. Loosening identity
from 0.95 to 0.70 does not touch it. Dropping coverage from 0.8 to 0.5 does not touch it.
`--cov-mode 2` (-0.519), `--cov-mode 5` (-0.518) and the post-cluster merge at 90% identity
over 90% of the shorter sequence (-0.519) are likewise indistinguishable from baseline.

Eleven of the twelve distinct clusterings tested (the table's 13 rows include the shipped
baseline and its byte-identical re-run) leave the confound exactly where it was. Family counts move
by a factor of two (30,406 to 60,138) and the core/accessory split moves substantially,
but the fragmentation signal rides straight through all of it.

**`--cov-mode 1` is the single exception**: rho(accessory ~ N50) **-0.519 -> -0.350**
(a 33% reduction) and rho(accessory ~ contigs) **+0.435 -> +0.290**. It also does exactly
what the hypothesis predicts a fix should do: singleton families drop 26,466 -> 14,833
(-44%), strain-private families 24,080 -> 13,774 (-43%), and the contig-terminal enrichment
of private-family proteins falls from **31.4% to 18.1%** - the truncation-derived
singletons really are being reabsorbed. The median private-family count in the 40 most
fragmented assemblies drops from **61 to 20**, i.e. down to the level the 40 best-assembled
assemblies show at baseline (21).

Note that even at its best the fix is partial: -0.350 is still a substantial confound, and
the private-family correlation barely moves (+0.408 -> +0.314).

## The paralog-merging cost, quantified - and it is disqualifying

`--cov-mode 1` requires the alignment to cover the target only, so a short sequence can be
absorbed into a long family regardless of how little of that family it explains. Measured
consequences on this dataset:

| measure | baseline (cov-mode 0) | cov-mode 1 |
|---|---|---|
| mean copies per present cell | 1.0006 | **1.0370** |
| families with any multi-copy strain | 1,535 | **4,331** |
| same-strain copy pairs in the whole matrix | 2,557 | **171,477 (67x)** |

A same-strain copy pair is ambiguous on its own: two halves of one contig-broken gene also
land in the same strain. So each of the 168,920 copy pairs that cov-mode 1 creates and the
baseline does not was classified by the contig position of its two members:

| class of newly-created same-strain copy pair (cov-mode 1) | n | share |
|---|---|---|
| **both members contig-terminal** (consistent with a broken gene) | 17,586 | **10.4%** |
| one member contig-terminal | 22,982 | 13.6% |
| **neither member contig-terminal** (two interior, intact genes) | 128,352 | **76.0%** |
| members on different contigs | 46,919 | 27.8% |

Three quarters of what `--cov-mode 1` merges is two interior, intact gene models being
collapsed into one family. That is not fragmentation repair; that is family collapse.

Pfam corroborates it. Of the 3,055 cov-mode-1 families that absorb two or more baseline
families where at least two of those baseline representatives carry a Pfam annotation
(existing `pfam.domtblout`, 9,163 representatives covered), **28.9% merge representatives
with completely disjoint domain architectures** and only 23.2% merge representatives with
identical architectures. The comparable figures for the alternatives that do not reduce the
confound are far lower: `--cov-mode 5` 10.2% disjoint, `min_id 0.70 / cov 0.50` 12.1%
disjoint. `--cov-mode 1` is roughly **3x** more likely to merge proteins with no shared
domain than any other setting tested.

## Tier 2 - downstream effect for the shortlist

Shortlist justification from tier 1: **`--cov-mode 1`** because it is the only setting that
reduces the confound at all; **`min_id 0.70 / cov 0.50`** as the loosest grid endpoint (the
largest family-count change available, 30,406 vs 54,421, to show what a big clustering
change does downstream *without* the confound benefit); and **`0.9/0.8/cov-mode 0`** re-run
through the same standalone chain so the comparison is like-for-like rather than against
the shipped pipeline run.

Chain run per setting: presence matrix -> frequency bins -> family positions ->
co-occurrence (Fisher + BH-FDR + exact stratified test) -> pair classification (k=10) ->
islands -> Leiden (resolution 1.0, seed 42). Wall clock: 1h14m (cov-mode 1), 1h29m
(0.70/0.50), 2h57m (baseline - it has the most families and therefore the most candidate
pairs, 110.9M vs 36.5M).

| tier-2 measure | 0.9/0.8 cov-mode 0 (baseline) | 0.9/0.8 cov-mode 1 | 0.70/0.50 cov-mode 0 |
|---|---|---|---|
| families | 54,421 | 33,703 | 30,406 |
| accessory islands | **27,836** | 11,079 (-60%) | 13,909 (-50%) |
| islands per 1,000 families | 511 | 329 | 457 |
| `trans` pairs (permutation-confirmed) | **966,865** | 397,699 (-59%) | 516,908 (-47%) |
| `trans_unconfirmed` pairs | 11,773,437 | 4,686,911 | 5,369,677 |
| `unexplained_physical` | 34,124 | 16,168 | 17,876 |
| Leiden modules (r=1.0) | 5 | 5 | 5 |
| largest module | 7,503 | 4,371 | 4,623 |
| NMI vs baseline modules (shared families) | - | 0.695 | 0.679 |
| ARI vs baseline modules | - | 0.779 | 0.773 |

Reading: island count and trans-pair count track the family count almost linearly - both
alternatives cut them roughly in proportion to the families they remove. Nothing here
argues that one partition is more correct than another; it says the downstream analysis
subject scales with the clustering, which is what makes the clustering parameter
consequential in the first place and is why it needed an independent criterion (the
fragmentation confound) rather than a downstream count to choose on.

The Leiden result is uninformative as a discriminator: all three settings give **5 modules
at resolution 1.0** with a giant first module, reproducing the known too-coarse-resolution
problem already recorded as #132 priority 4. NMI 0.68-0.70 and ARI 0.77-0.78 between each
alternative and the baseline mean the coarse module structure is broadly conserved under a
2x change in family count - so module identity is not sensitive enough at r=1.0 to
distinguish clustering settings either. Resolve the resolution question (priority 4) before
using module structure as evidence for anything.

## Recommendation

**Keep `pangenome_tier1_min_id = 0.9` and `pangenome_tier1_cov = 0.8`, with cov-mode 0.**

Not because those values were shown to be optimal - they were not, and nothing in this
sweep distinguishes 0.9/0.8 from 0.8/0.8 or the post-merge variant on the fragmentation
measure - but because **no setting tested buys a meaningful reduction in the fragmentation
confound at an acceptable cost**, and the one that reduces it does so by destroying
biological signal.

**Evidence class: empirical, negative.** The recommendation is "the parameter is not the
lever" rather than "this value is right". State it that way in
`docs/pangenome-assumptions.md`: 0.9/0.8 is a convention that this sweep failed to improve
on, not a validated optimum. What the sweep does establish positively is that the
fragmentation confound is **not a clustering-parameter artefact** - it survives a 2x change
in family count - so it has to be addressed elsewhere.

### What to do instead

The confound is an assembly-quality property of the input, not of the clustering. Three
routes, none of them a `-c` value:

1. **Filter at the input, not the cluster step.** A private-family protein sitting at a
   contig end is identifiable directly from `gene_positions.tsv` (31.4% of them do, vs
   9.4% background) with no re-clustering. Flagging or excluding contig-terminal singleton
   families attacks the artefact without touching what a family means for every other
   protein. This is the same structural criterion issue #133 proposes for rescue, applied
   one step earlier.
2. **Control for assembly quality in the analysis** rather than removing it from the data -
   report accessory counts as residuals against N50/contig count, or restrict
   fragmentation-sensitive claims to an assembly-quality-matched subset.
3. **If a coverage change is wanted anyway**, `--cov-mode 5` (coverage of the shorter
   sequence) or the post-cluster merge is the least harmful: 10.2% disjoint-domain merges
   vs cov-mode 1's 28.9%, copies/cell 1.0011 vs 1.0370. Both reduce the family count
   ~10% and both leave the confound at -0.518/-0.519, so neither buys anything. They are
   safe, not useful.

### Risk statement (asked for explicitly)

Looser coverage handling does merge genuine paralogs, and this dataset can measure it:
`--cov-mode 1` raises within-strain copy pairs 67x, and 76% of those merges join two
interior gene models rather than two contig-end fragments, with 28.9% joining proteins of
disjoint Pfam architecture. Adopting it to improve the N50 correlation would trade a
known, measurable confound for an unmeasurable deflation of real copy-number and
gene-family variation. Do not adopt it.

### Caveats

- Single taxon. Coccidioides is clonal and low-diversity; the paralog-merge cost of loose
  clustering could be quite different in a recombining or paralog-rich genus. The note in
  #132 about confirming on other taxa stands.
- The Pfam concordance test uses the existing island-scan `pfam.domtblout`, which covers
  9,163 of 54,421 baseline representatives (17%). The disjoint-domain fractions are
  relative comparisons between settings on the same covered subset, not absolute rates.
- "Contig-terminal" is gene-index terminality (first/last gene model on the contig), not a
  measured distance to the contig end in bp. It is a conservative proxy: it will miss a
  truncated gene that happens to have a neighbour beyond it.
- `--cluster-reassign` was on for every run, matching the pipeline.

## Artefacts

All under `/rhome/jstajich/.claude/jobs/3354f7d2/tmp/tier1sweep/`:
`runs/<setting>/tier1_cluster.tsv.zst` (12 clusterings, 405 MB total),
`metrics/<setting>.json` (tier-1 metrics), `tier2/<setting>/` (tier-2 outputs),
`logs/` (SLURM logs and `/usr/bin/time -v` resource records),
`metrics.py`, `terminal.py`, `discriminate.py`, `pfam_merge.py`, `tier2_job.sh`
(the analysis code). Nothing was written into either repository.
