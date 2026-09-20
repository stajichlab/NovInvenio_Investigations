# Genome relatedness and phylogrouping for the NovInvenio pangenome pipeline

**Status**: investigation / method recommendation. No pipeline code was changed.
**Date**: 2026-09-20
**Evidence base**: the real 529-strain *Coccidioides* study
(`NovInvenio_Investigations/studies/fungi/coccidioides_pangenome`, run
`mmseqs_genus_vs_ureesii`), plus new read-only mash runs done for this investigation.
**Scratch / reproducibility**: all new data and scripts under
`/rhome/jstajich/.claude/jobs/3354f7d2/tmp/phylo/` (`a1`–`a17` scripts, `sketch.sbatch`,
`splithalf.sbatch`). Nothing was written into either repo.

This document is written to **generalize to other pangenome datasets**. Section 1
describes the current strategy, Section 2 how well it works and how it fails, Sections
3–6 answer the four scientific questions, Section 7 gives parameter guidance and a
validation checklist for a new dataset, Section 8 ranks the recommendations.

---

## 1. Current strategy (what the pipeline does today)

Three mash-based steps, all in `workflows/pangenome_profile.nf` and
`modules/pangenome/mash.nf`.

### 1.1 `MASH_SKETCH`
Run **twice** per study, over two different genome sets:

- `MASH_SKETCH_ALL` — every strain (ingroup + outgroup). Feeds `DEREPLICATE`.
- `MASH_SKETCH_INGROUP` — ingroup only. Feeds `ASSIGN_CLADES`.

Both call `mash sketch -p <cpus> -o <prefix> <genomes>` then
`mash dist -p <cpus> -t <prefix>.msh <prefix>.msh`. **No `-s` or `-k` is passed**, so
mash defaults apply: sketch size `s = 1000`, k-mer `k = 21`. The split into two sketches
is deliberate: including an outgroup species swamps the within-ingroup signal (the
between-species distance is ~5x the within-species distance — measured below).

### 1.2 `DEREPLICATE` → `strain_inventory.tsv`
`bin/pangenome_dereplicate_strains.py` thresholds the all-strain mash matrix at
`params.pangenome_mash_threshold` (default **0.001**), forms transitive-closure
(single-linkage) dedup groups, and picks one representative per group by an
assembly-quality proxy (contig count / N50 / length). Output columns:
`Short, n_contigs, n50, total_length, dedup_group, is_representative`.

**Important**: the inventory does **not** gate any compute. Its only consumers are
`bin/pangenome_frequency_bins.py` and `bin/pangenome_cooccurrence.py`, which use it as an
*optional denominator* (count only `is_representative == 1` strains). Every per-strain
process still runs on every strain.

### 1.3 `ASSIGN_CLADES` → `clade_assignments.tsv` → `FILL_TAXON_GROUP`
`bin/pangenome_assign_clades.py`: ingroup mash distance matrix → **classical PCoA**
(double-centred Gower transform + eigendecomposition, negative eigenvalues dropped) →
keep `params.pangenome_clades_n_pcoa` (default 4) components → **k-means** (`scipy
kmeans2`, `minit='++'`) over `k` in `params.pangenome_clades_k_range` (default `2,10`) →
pick the `k` with the best **mean silhouette coefficient**. `params.pangenome_clades_k_fixed`
skips the sweep. `FILL_TAXON_GROUP` writes the labels into the samplesheet's *empty*
`TaxonGroup` cells only — a curated value is never overwritten.

### 1.4 Downstream consumption of the clade labels
Exactly two places, both discrete:

- `bin/pangenome_cooccurrence.py::exact_stratified_pvalue()` — an **exact** within-clade
  stratified association test. Shuffling family B's presence within each clade fixes every
  margin, so the null overlap count is a sum of independent `Hypergeom(n_c, |a_c|, |b_c|)`
  draws, one per clade, convolved. No simulation. This is a Cochran–Mantel–Haenszel-style
  test stratified by clade.
- `bin/pangenome_pair_classification.py` — a pair is called `trans` only if
  `linkage_fraction <= trans_threshold` **and** `permutation_p < perm_alpha` **and** the
  pair's `clade_composition` spans `>= params.pangenome_pair_class_min_clades` (default 2)
  distinct clades. Otherwise `trans_unconfirmed`.

### 1.5 Where it is configured (`nextflow.config`, ~lines 218–232, 259–264)
```
pangenome_dereplicate     = true
pangenome_mash_threshold  = 0.001
pangenome_assign_clades   = true
pangenome_clades_k_range  = '2,10'
pangenome_clades_k_fixed  = null
pangenome_clades_n_pcoa   = 4
pangenome_pair_class_min_clades = 2
```
There is **no tree anywhere in the pipeline** and no continuous relatedness metric is
carried past `clade_assignments.tsv`'s PCoA coordinate columns (which nothing reads).

---

## 2. How well it works, and how it fails

### 2.1 It works at the species level, decisively
On the 529-strain *Coccidioides* ingroup the silhouette sweep chose **k = 2**, giving
527/529 = **99.6%** concordance with the curated `Species` column (the 2 exceptions are
metadata errors — Section 6).

New measurements confirming the species split is not marginal:

| quantity (mash `s=100000`, this investigation) | value |
|---|---|
| PCoA axis 1, whole ingroup, variance explained | **80.8%** |
| PERMANOVA R² of the species label on mash distance | **0.907** |
| median within-species distance (immitis / posadasii) | 0.0033 / 0.0048 |
| median between-species distance | 0.0217 (**~5x**) |
| split-half reliability of the k=2 partition (Section 3.2) | **ARI = 1.000** |

So the **minimum bar the user named is met with a very large margin**, and it is
reproducible, not a lucky cut.

### 2.2 Where it will work worse on another dataset — the failure modes

1. **Single-species study.** With one species there is no dominant axis. The silhouette
   sweep will still return some `k` (silhouette is defined for any partition and always
   picks *something*) and `FILL_TAXON_GROUP` will write it into `TaxonGroup` with no
   warning. Measured here: within either *Coccidioides* species alone, every k-means /
   Ward partition at k = 2–10 has **split-half ARI 0.02–0.26** — i.e. essentially
   irreproducible. A single-species study gets a confident-looking but meaningless
   `TaxonGroup` column, and `min_clades >= 2` then filters `trans` pairs on noise.
   This is the most dangerous failure mode and it is silent.
2. **A genus with less divergence between the groups you care about.** The success above
   rests on a 5x distance ratio. Two *recently* diverged species, or two host-associated
   populations, will not produce that; the split will be a cut through a continuum.
3. **Small N.** k-means + silhouette over k=2..10 on 20 genomes is unreliable, and the
   exact stratified test loses all power once strata get small (a stratum with
   `|a_c| == 0` or `|b_c| == 0` contributes nothing; a 2-strain stratum can barely
   contribute). Below roughly 40–50 ingroup genomes, stratification beyond 2–3 groups is
   counterproductive.
4. **Highly fragmented / uneven assemblies.** Mash distance itself is fairly robust
   (measured split-half `r = 0.98` overall), but **gene presence/absence is not**: within
   dedup groups, gene-content Hamming distance correlates with `|Δ n_contigs|` at
   Spearman **0.55** and 48% of the differing families are unique to the *smaller*
   assembly (Section 5). A study mixing long-read and short-read assemblies will get
   assembly-quality-driven "population structure" in the gene-content data even though
   the mash grouping is fine.
5. **More than 2 real groups with very different sizes.** Silhouette is biased toward
   balanced, well-separated, globular clusters and toward small `k`. It will under-report
   a 3rd group that is small or elongated.
6. **`s = 1000` sketch noise.** See 2.3 — this limits fine structure, though it is not the
   binding constraint here.

### 2.3 The sketch size is 100x smaller than it needs to be, and fixing it is free
With `mash sketch` defaults, a distance is a count of shared hashes out of 1000, so it is
**quantized in steps of ~2.7e-5**. Across all 139,656 ingroup pairs there are only **323
distinct distance values**. Within-species distances sit at 0.003–0.005, i.e. resolved to
roughly 110–180 levels, each carrying the binomial sampling noise of 1000 draws.

Re-sketching all 529 genomes at higher `s` (job `28951656`, `short` queue, 24 cores):

| sketch size | distinct distance values | sketch+dist wall time (24 cores) |
|---|---|---|
| `s=1000` (current) | 323 | ~3 min (from the study's own trace) |
| `s=10000` | 2,322 | ~2 min |
| `s=100000` | **16,092** | **~2.5 min (2 min sketch + 34 s dist)** |

`MASH_SKETCH_INGROUP` took 2.7 min and `MASH_SKETCH_ALL` 3.1 min in the real run, against
**18,357 CPU-minutes (306 CPU-hours) for `TBLASTN_PER_STRAIN`**. Raising `-s` to 100000
changes the study's total cost by well under 1%. There is no cost argument for leaving it
at the default.

What it buys, measured: the s=1000 clustering of within-species structure is **not
reproducible against the s=100000 clustering** (ARI 0.30–0.61 at k=3..8), i.e. a
meaningful share of the apparent sub-structure at s=1000 is sketch noise. What it does
**not** buy: real sub-species clusters (Section 3) — within-species PCoA1 rises only from
12.5% to 14.5% variance explained. Bump `-s` because it removes an avoidable noise source
and makes the distances usable as a continuous metric, not because it will find clusters.

---

## 3. Q1 — Is there finer structure that k=2 is hiding?

### 3.1 Honest answer: no discrete sub-species clusters exist in this dataset
Splitting by species first and clustering within each species (the hierarchical approach)
*does* produce partitions, and their silhouette scores look respectable (immitis best
k=5, sil 0.42; posadasii best k=2, sil 0.49). Those numbers are misleading. Three
independent lines of evidence say the within-species partitions are not real structure:

**(a) The within-species pairwise distance distributions are unimodal.** Both species show
a single dominant mode (immitis ~0.0032, posadasii ~0.0048) with a long *left* tail of
near-clonal pairs. There is no second mode, i.e. no set of pairs that are "the same
distance as each other but clearly further than everyone else" — the signature a genuine
sub-population split would leave.

**(b) Any within-species grouping explains little variance and barely separates.**
PERMANOVA on the mash matrix (`s=100000`):

| grouping | R² | median within / between distance | ratio |
|---|---|---|---|
| species (whole ingroup) | **0.907** | 0.0040 / 0.0217 | **5.4** |
| Ward k=2 within immitis | 0.108 | 0.00326 / 0.00335 | 1.03 |
| Ward k=8 within immitis | 0.381 | 0.00275 / 0.00334 | 1.21 |
| Ward k=2 within posadasii | 0.131 | 0.00456 / 0.00491 | 1.08 |
| Ward k=8 within posadasii | 0.279 | 0.00422 / 0.00480 | 1.14 |

Two strains in *different* within-species "clusters" are only 3–21% further apart than two
in the same one. (The PERMANOVA p-values are all 0.005 and are **meaningless** — the
groupings were derived from the same distance matrix, so the test is circular. R² and the
distance ratio are the informative quantities and they are reported above.)

**(c) Split-half reliability kills it.** This is the decisive test. Each genome's contigs
were split alternately into two halves (odd contigs → half A, even → half B), each half
sketched independently at `s=50000` and a distance matrix built from each (job `28951717`,
~6 min wall). Two near-independent estimates of the same relatedness; a real cluster must
appear in both.

| partition | split-half ARI |
|---|---|
| **whole ingroup, k=2 (the species split)** | **1.000** |
| within immitis, Ward k=2 / 3 / 5 / 10 | 0.056 / 0.021 / 0.123 / 0.164 |
| within posadasii, Ward k=2 / 3 / 5 / 10 | −0.007 / 0.229 / 0.261 / 0.256 |

The distance *values* are reliable (within-species distance r between halves = 0.88; all
pairs r = 0.98). The *partitions* built on them are not, because there are no gaps to cut
at. **The 4–8 clusters are not recoverable from mash on these genomes.**

### 3.2 Is silhouette the wrong criterion? Partly, but it is not the problem
- Gap statistic on the whole ingroup also chooses k=2; within immitis it chooses k=1
  ("no structure"), within posadasii k=5. So the criteria disagree within species —
  itself a symptom of no real structure.
- UPGMA (average linkage) on the within-species distance matrix does not split at all: at
  every k from 2 to 8 it just peels off singleton outliers (largest cluster 158–168 of
  169). That is the classic topology of a near-star phylogeny.
- Changing the criterion would only change *which* unreproducible partition you report.

### 3.3 Verdict on DAPC
**Not worth adopting here, and it would not have given the answer the user remembers.**
- DAPC (Jombart et al. 2010) is *discriminant* analysis of principal components: it finds
  axes that maximise **between-group** relative to within-group variance, *given* groups.
  The groups normally come from `find.clusters` (k-means over retained PCs, selected by
  **BIC**, not silhouette). So DAPC is not a different way of discovering structure — the
  discovery step is the same k-means-on-PCs the pipeline already does. What DAPC adds is a
  *projection* that makes whatever groups you fed it look maximally separated. Run on
  unreproducible clusters it will draw beautiful, well-separated ellipses around noise.
  This is a well-known DAPC failure mode, and it is exactly the trap here.
- DAPC's normal input is **genome-wide SNP allele frequencies**, not a k-mer sketch
  distance. The resolution difference is the whole story: *A. fumigatus* clade structure
  is resolved from ~hundreds of thousands of SNPs; mash `s=1000` gives 323 distinct
  distance values over the entire dataset. If 4–8 *A. fumigatus* clusters were recovered,
  it was from variant data (or from a species with genuinely deeper, more discrete
  population structure — *A. fumigatus* clades A/B are well separated), not from a
  minhash sketch of this resolution. I cannot reproduce a minhash route to 4–8 clusters,
  and the split-half test says one does not exist for *Coccidioides*.
- **If sub-species clusters are genuinely wanted**, the route is variant calling (map all
  strains to a reference, call SNPs, then DAPC/ADMIXTURE/fineSTRUCTURE on the SNP matrix),
  or a real marker-gene phylogeny (`nf_phyling` → BUSCO markers → IQ-TREE). That is a
  separate study, not a pangenome pipeline step. See Section 8c.

---

## 4. Q2 — A continuous relatedness metric

### 4.1 Recommended metric
**Mash distance at `-s 100000 -k 21`, published as a full strain × strain matrix
(`mash_distance.tsv`) alongside `clade_assignments.tsv`.**

Justification against the alternatives:

| candidate | verdict |
|---|---|
| **mash, `s=100000`** | **Recommended.** Already in `pixi.toml` and in the container. Measured cost for 529 genomes: 2.5 min on 24 cores, <1% of study runtime. Split-half `r = 0.98`. Only a parameter change to a process that already runs. |
| mash, `s=1000` (current) | Adequate for the species split, too coarse to publish as a continuous metric (323 distinct values, step 2.7e-5). |
| `skani` / `fastANI` (true ANI) | Better-calibrated ANI, and `skani` is fast. But **not in `pixi.toml` and not installed** — a new dependency and a new container build. The measured limitation here is not distance accuracy (r=0.98 split-half) but the absence of real structure; ANI would not change that. Worth it only if a future study needs a defensible species-boundary ANI number (the 95% ANI convention). |
| `sourmash` | Same family as mash, extra dependency, no benefit over `mash -s 100000`. |
| `mashtree` | A neighbour-joining tree from mash distances. It would produce a tree, but the tree's within-species topology inherits exactly the instability measured in Section 3.2. **Do not use it to justify gain/loss polarisation** — a tree with unreproducible internal branches is worse than no tree, because it looks authoritative. |
| **real marker phylogeny** (`nf_phyling` BUSCO → PhyKIT concat → IQ-TREE) | The only thing that would genuinely support branch-level ancestral state reconstruction. Costly at 529 genomes (BUSCO per genome is the bottleneck; realistically hours to a day of cluster time even parallelised). Correct for a *dedicated* phylogenomics study; too slow as an unconditional per-study pipeline step. Make it opt-in. |
| kinship / covariance matrix (GRM) | Requires SNPs. Same blocker as DAPC. |

### 4.2 (a) Does a better structure correction improve co-occurrence statistics?
**Measured, and the answer is no — the current 2-clade stratification is already
sufficient for this dataset.** This is the most important negative result in this report.

First, the confounding is real and large. Accessory gene content (12,582 families with
2% < frequency < 95%) tracks relatedness strongly:

| pair set | Spearman(mash distance, gene-content Jaccard) |
|---|---|
| all pairs | 0.886 |
| within-species pairs | **0.726** (linear R² = 0.57) |
| within immitis / within posadasii | 0.466 / 0.591 |

Within species, a linear model of Jaccard on mash distance + `|Δ n_contigs|` reaches
R² = 0.706, of which mash alone gives 0.570 and assembly contiguity alone 0.164.

Second, I measured the *calibration* of the actual test the pipeline uses. 300 random
pairs of real accessory families (5% < f < 50%), which carry real phylogenetic
autocorrelation but have no reason to be biologically associated. Fraction with p < 0.05
(should be ~0.05):

| test | false-positive rate at α = 0.05 |
|---|---|
| plain Fisher, no stratification | **0.327** |
| exact stratified, k=2 (species — what the pipeline does) | **0.050** |
| exact stratified, k=4 hierarchical | 0.040 |
| exact stratified, k=10 | 0.037 |
| exact stratified, k=20 | 0.037 |
| exact stratified, k=40 | 0.040 |
| exact stratified, k=80 | 0.033 |

Unstratified, **one in three random family pairs is "significantly co-occurring"**. The
existing k=2 stratification fixes it exactly. Going finer buys ~0.01 and starts to become
conservative (losing real power) by k=80.

Third, restricting to one species at a time — which is what a *single-species* study looks
like — shows why:

| within immitis (n=169) | FPR | within posadasii (n=360) | FPR |
|---|---|---|---|
| plain Fisher | 0.048 | plain Fisher | 0.072 |
| stratified k=2 | 0.044 | stratified k=2 | 0.052 |
| stratified k=10 | 0.036 | stratified k=10 | 0.048 |
| stratified k=20 | 0.036 | stratified k=20 | 0.040 |

(250 pairs each; binomial SE ≈ 0.014.) **All of the confounding lives on the species
axis.** Within a species, the unstratified test is already close to nominal.

Consequences, stated plainly:
- **Do not** replace `exact_stratified_pvalue` with a phylogenetically-aware permutation,
  a mixed model, or relatedness-weighted permutation. Measured benefit on this data:
  approximately zero. All three are substantially more code, slower, and harder to
  defend. This is the clearest "sounds sophisticated, not worth it" item in this report.
- **Do** keep the stratification. Removing it would be catastrophic (0.327 FPR).
- **Do** add the guardrail in Section 7.3: a study with only one real group gets a
  `TaxonGroup` column of noise, and `min_clades >= 2` then filters `trans` on noise.
  Detecting that is worth more than any refinement of the statistic.

### 4.3 (b) Polarising gain vs loss
`polarize_direction()` calls direction from raw outgroup presence counts: present in
*every* outgroup strain → "loss"; absent from *all* → "gain"; mixed → "ambiguous". No
phylogenetic context at all (`todo/pangenome-phylogeny-aware-gain-loss.md`).

**A rough tree does not fix this; it needs a properly rooted one.** Direction is
inherently a rooting question: "gain on this branch" and "loss on the complementary
branch" are the *same* unrooted pattern. What makes the current rule work at all is that
it is implicitly rooted on the outgroup. So:

- A **mashtree** would add nothing — it gives an unrooted topology whose within-species
  branches are, per Section 3.2, not reproducible.
- What *would* add real information is a **rooted marker-gene tree** (`nf_phyling` →
  IQ-TREE with the outgroup as the root) plus **Dollo or Fitch/Sankoff parsimony** over
  the presence/absence characters, counting independent gain and loss *events* per family
  instead of one binary call. Tools: `Count`, `BadiRate`, or a direct Sankoff
  implementation.
- **But**: at the species level, which is the only level this dataset resolves
  reproducibly, the tree is `((immitis, posadasii), ureesii)` — which the existing
  outgroup rule already encodes. A real tree would only change the answer where
  *within-species* branch placement matters, and that is precisely the resolution mash
  cannot deliver. So the tree must come from markers/SNPs or not at all.
- **Prerequisite before believing any gain/loss number**: this run's presence matrix is
  the *unrescued* one (the `RESCUE_PASS` zstd-symlink bug — `presence_matrix.tsv` and
  `presence_matrix.rescued.tsv` are byte-identical, zero `genome_only` cells). At 15.6%
  density, an "absence" is currently "no protein-level cluster member", not "not in the
  genome". The reported 99.3% gain / 0.7% loss split is what you get when absences are
  systematically over-called: a family that is genuinely present but unannotated in a
  strain looks like an absence, which inflates apparent "gain" and suppresses "loss".
  **Re-run with the rescue fix before any phylogenetic polarisation work.** Investing in a
  tree to polarise a matrix with a known systematic absence bias would be spending the
  expensive resource on the unreliable input.

---

## 5. Q3 — "Clusters that move between differently related clades"

### 5.1 What the pipeline asks today
`trans` requires `clade_composition` to span `>= min_clades` (default 2) clades. With
k=2 clades that is simply "carried by at least one strain of each species". It is a
presence test, not a statistic: it has no null model, no notion of *how* the carriers are
distributed, and no way to distinguish one ancient gain shared by both species from many
independent recent acquisitions.

### 5.2 The better question, and the statistic
The right question is **phylogenetic incongruence**: is a family's presence pattern
*more scattered* across the relatedness structure than a vertically inherited family of
the same frequency would be? Scattered presence across distant genomes, with close
relatives lacking it, is the classic HGT / mobile-element signature.

Concrete statistic, computable from the mash matrix alone, no tree needed:

> For family *f* with carrier set *C*, |C| = m ≥ 5:
> `D_obs(f) = mean pairwise mash distance among carriers`.
> Null: draw m strains uniformly at random, recompute; z-score against that null
> (60+ draws per carrier-set size).
> `z >> 0` → over-dispersed (carriers unusually distant — HGT / mobile element candidate).
> `z << 0` → under-dispersed (clonally inherited, lineage-restricted).

I implemented and ran this (`a15_disp.py`, 16,783 testable families). It runs in seconds
(one `A @ M` matrix product). Result distribution: median z = −5.3, 71% of families at
z < −3 (strongly lineage-restricted — as expected for vertical inheritance), 5.3% at
z > +2.

**However, run as-is on this two-species dataset it is dominated by the species axis and
is not yet useful.** The most under-dispersed families are simply posadasii-specific
(immitis_frac = 0.00); the most over-dispersed are families present in ~40% immitis /
~60% posadasii, i.e. shared genes with patchy calls. It is re-discovering "which species
is this gene in", plus presence/absence noise from the unrescued matrix.

**Fixes needed before this is a deliverable statistic** (this is why it is a ticket, not a
do-now):
1. **Compute the null within species, not across.** Stratify the random carrier draw by
   species so the species axis cannot drive z — exactly the same stratification logic
   already in `exact_stratified_pvalue`.
2. **Re-run on the rescued matrix.** At 15.6% density with over-called absences, a
   patchy pattern is as likely to be annotation dropout as biology. This is the single
   biggest blocker.
3. **Condition on assembly quality.** Gene-content Jaccard correlates with
   `|Δ n_contigs|` at Spearman 0.31 within species; a family absent only from fragmented
   assemblies will look over-dispersed.
4. **Report it per family** as a new `phylo_dispersion_z` column on the frequency table,
   and **replace the binary `min_clades` gate in `pair_classification` with it** (or add
   it alongside): a `trans` pair is far more interesting when both families are
   over-dispersed than when they merely both happen to occur in both species.

---

## 6. Q4 — Collapsing near-identical genomes

### 6.1 Current state, measured
On the real run: **530 strains → 412 dedup groups**, i.e. **118 strains (22%) collapsed**
at mash 0.001. 56 multi-member groups; the largest have 19 and 17 members. And again:
**this saves no compute today**, because `DEREPLICATE` runs at step 4 of the workflow,
*after* the per-strain `TBLASTN_PER_STRAIN` fan-out at step 3, and only feeds the
frequency/co-occurrence denominators.

### 6.2 The threshold is about right — and it sits on a plateau
Single-linkage (the same transitive-closure semantics the current dedup uses):

| threshold | groups | strains collapsed | notional TBLASTN saving (of 306 CPU-h) |
|---|---|---|---|
| 0.0005 | 437 | 92 (17%) | ~53 CPU-h |
| **0.001 (current)** | **409** | **120 (23%)** | **~69 CPU-h** |
| 0.002 | 346 | 183 (35%) | ~106 CPU-h |
| 0.003 | 225 | 304 (57%) | ~176 CPU-h |
| 0.005 | **2** | 527 (100%) | — |

(mash `s=100000`; `s=1000` gives essentially the same picture.) The collapse to 2 groups
at 0.005 is **single-linkage chaining** — at that threshold the transitive closure eats
each entire species. 0.001 is a safe operating point; 0.002 is defensible; **0.003 and
above is dangerous with single linkage** and should be blocked or warned about. If the
threshold is ever to be raised, switch to complete or average linkage first so chaining
cannot silently merge a whole species.

### 6.3 The user's second clause — "and no variation in clusters present" — does not hold
This is the key finding for Q4. Strains inside a dedup group (mash < 0.001, i.e. genomes
that are ~indistinguishable) **still differ substantially in gene-family content**:

| comparison | mean gene-family Hamming distance (of 54,421 families) |
|---|---|
| pairs within a dedup group (n = 477 pairs) | **1,569** (median 1,500; min 1,009; max 3,913) |
| random strain pairs | 3,218 (median 2,720) |

Even the closest pair in the dataset (mash 2.4e-5) differs in 1,437 family calls.
**49.6% of all families vary within at least one dedup group** — including **61.7% of
"core" (f ≥ 95%) families** and **97.9% of softcore**, which is diagnostic: a core gene
cannot genuinely be absent from a near-clone.

This variation is **largely technical, not biological**:
- Hamming vs `|Δ n_contigs|`: Spearman **0.545** (p = 2.6e-38)
- Hamming vs `|Δ total_length|`: Spearman **0.406**
- **48%** of the differing families are present only in the *smaller* assembly of the pair
  — i.e. the more fragmented/shorter assembly is missing them.

And this is exactly what `RESCUE_PASS` exists to fix — and exactly what silently did not
happen on this run (zero cells rescued, zstd-symlink bug). So:

**Answer to "is a mash-distance threshold the right criterion, versus one based on actual
pangenome content variation?" — mash is the right criterion, and content is the wrong
one, *for now*.** Gene-content variation between near-clones is currently dominated by
assembly and annotation artefacts, so a content-based collapse rule would collapse
strains according to how similar their *assemblies* are, not their genomes. Revisit only
after a rescued matrix is in hand; the same measurement (mean within-dedup-group Hamming,
and the fraction of *core* families that vary within a group) is the natural acceptance
test for whether the rescue worked.

### 6.4 What collapsing more would actually save
Runtime profile of the real run (from `trace.txt`, 530 strains):

| process | tasks | total realtime |
|---|---|---|
| `TBLASTN_PER_STRAIN` | 530 | **18,357 min ≈ 306 CPU-h (94% of total)** |
| `PAIR_CLASSIFICATION` | 1 | 93 min |
| `COOCCURRENCE` | 1 | 84 min |
| `PREFIX_GENOME` / `PREFIX_PROTEOME` / `MAKE_STRAIN_GENOME_DB` | 530 each | 39 / 34 / 28 min |
| `MASH_SKETCH_ALL` / `_INGROUP` / `DEREPLICATE` | 1 each | 3.1 / 2.7 / 3.1 min |

Cost is ~linear in strain count and concentrated in one per-strain process.

- **Already-available saving, currently unrealized**: gating the per-strain fan-out on
  `is_representative` at the *existing* 0.001 threshold would drop 118 of 530 tasks —
  **~69 CPU-hours, 22%**, for a DAG reordering (move `MASH_SKETCH_ALL`/`DEREPLICATE` ahead
  of the per-strain block).
- Raising to 0.002 adds another ~11%.
- **Caveat that must be respected**: collapsing at the *TBLASTN* stage means the collapsed
  strains lose their presence-matrix columns entirely. You cannot both run the analysis on
  530 strains and skip TBLASTN for 118 of them. So this is a `pangenome_dereplicate_analysis`
  decision (run the whole analysis on representatives only), not a free optimisation —
  it must be an explicit, separate opt-in flag, and the report must state which strain
  set the results describe.
- `COOCCURRENCE`/`PAIR_CLASSIFICATION` scale with families², not strains, so dereplication
  barely touches them.

---

## 7. Parameter guidance and a validation check for a new dataset

### 7.1 Choosing the parameters
| parameter | guidance |
|---|---|
| `mash sketch -s` | Set **100000** (currently the mash default 1000 is used). Costs ~2 min on 24 cores at 500 genomes. `-k 21` is fine for fungal genomes; raise to 23–25 only for genomes > 200 Mb. |
| `pangenome_clades_k_range` | Leave `'2,10'`. The sweep is cheap and the range is not the problem; the *criterion* is. Always inspect the per-k silhouette curve, not just the winner. |
| `pangenome_clades_k_fixed` | **Use it whenever you already know the grouping** (e.g. a two-species study → `k_fixed 2`; a study with curated DAPC/MLST clades → fill `TaxonGroup` by hand and let `FILL_TAXON_GROUP` leave it alone). A known grouping always beats a discovered one. |
| `pangenome_clades_n_pcoa` | Leave 4. Increasing it does not help: within-species PCoA axes carry 5–15% variance each with no elbow, so extra axes add noise. Only raise it if PCoA1 explains < ~40% *and* you have independent evidence of > 4 real groups. |
| `pangenome_mash_threshold` | Leave **0.001**. 0.002 is defensible. **Do not exceed 0.002** while dedup uses single linkage — 0.005 collapsed this entire dataset to 2 groups by chaining. |
| `pangenome_pair_class_min_clades` | Only meaningful if the clade labels are real. With 2 real clades, 2 is right. **With one real group, set `pangenome_assign_clades = false` and do not use `min_clades` at all** rather than accept a noise partition. |

### 7.2 One-line summary for a new dataset
> Sketch at `-s 100000`; trust the clade grouping **only** if PCoA1 explains > ~50% of the
> distance variance, the k=2..3 partition survives a split-half ARI > 0.8, and the
> between-group / within-group median distance ratio is > ~2 — otherwise set
> `pangenome_assign_clades = false`, treat the study as unstratified, and keep
> `pangenome_mash_threshold` at 0.001.

### 7.3 A concrete worked check to run before trusting downstream results
Five diagnostics, all cheap, all runnable from the mash matrix the pipeline already
produces. Scripts `a1`, `a7`, `a10`, `a14`, `a16` in the scratch dir are working
implementations.

1. **PCoA axis-1 variance.** From `mash dist -t`, double-centre and eigendecompose. If
   axis 1 carries **> 50%**, you have a dominant real split. Here: 80.8% → trustworthy.
   If the top axes are all < 20% with no elbow, there is no discrete structure — stop and
   set `pangenome_assign_clades = false`.
2. **Distance histogram.** Plot the within-group pairwise distances. **Bimodal → real
   split. Unimodal with a left tail → a continuum plus clones**, and any k-means cut is
   arbitrary. Both *Coccidioides* species are unimodal within species.
3. **Between/within ratio.** Median between-group distance ÷ median within-group distance
   for the chosen partition. **> 2 is a real grouping**; 1.0–1.2 is a cut through a
   continuum. Here: species = 5.4 (real); every within-species partition = 1.03–1.21
   (not real).
4. **Split-half ARI — the decisive test.** Split each genome's contigs odd/even, sketch
   each half independently, cluster each, compare with adjusted Rand index. **ARI > 0.8 →
   the partition is reproducible.** Cost measured here: ~6 minutes on 24 cores for 529
   genomes. Here: species k=2 → 1.000; every within-species partition → 0.02–0.26.
5. **Label concordance QC.** For each strain, compare its mean mash distance to its own
   labelled group against its mean distance to every other group; flag any strain that is
   closer to a different group. On this dataset it flagged **exactly 2 of 529** (Section
   8), with ratios 0.16 and 0.22 against a 1st-percentile background of 3.83 — no false
   positives and a huge margin.

A sixth check, for the *gene content* rather than the grouping: compute mean gene-family
Hamming distance within dedup groups, split out by frequency bin. **If > ~50% of `core`
families vary between near-clones, your presence/absence matrix has an annotation/assembly
dropout problem** (here: 61.7%, caused by the no-op rescue pass) and every downstream
gain/loss and co-occurrence number should be treated as provisional.

---

## 8. The two discordant strains

### 8.1 Verdicts

**`485B-1_L_OLD_CPA0023` — labelled *C. immitis*. Verdict: MISLABELLED. It is
*C. posadasii*. Confidence: very high.**
- Mean distance to immitis 0.02068, to posadasii 0.00458 (**4.5x closer**).
- Decisive evidence: its three nearest genomes are **its own sibling isolates** —
  `485B-4_S_NEW_CPA0029` at **0.000072**, `SOIL_485B-0_L_OLD_CPA0020` and
  `485B-0_S_OLD_CPA0021` both at 0.000192 — and all three are labelled *C. posadasii*.
  These are the same 485B isolate series (an OLD/NEW passage or soil/clinical set). A
  distance of 7.2e-5 is at the dataset's resolution floor: this is the same strain
  background as genomes everyone agrees are posadasii.
- **Not a hybrid**: a hybrid would sit intermediate between the two clouds. It is not
  intermediate at all — it is *inside* the posadasii cloud (0.00458 vs a posadasii
  within-species median of 0.0048), i.e. a perfectly typical, central posadasii genome.
- **Not contamination**: assembly stats are unremarkable and match its siblings
  (581 contigs, 27.60 Mb, vs 572 contigs / 27.62 Mb for `485B-4` and 524 / 27.54 Mb for
  `485B-0_S`). A mixed immitis/posadasii assembly would be inflated in length and would
  show *reduced* distance to immitis; neither is true.

**`B3476` — labelled *C. posadasii*. Verdict: MISLABELLED. It is *C. immitis*.
Confidence: very high.**
- Mean distance to posadasii 0.02111, to immitis 0.00346 (**6.1x closer**), and 0.00346 is
  slightly *below* the immitis within-species median of 0.0033–0.0035 — i.e. it is a
  central, not peripheral, member of the immitis cloud.
- Nearest genome: `Michoacan_2` (immitis) at **0.00112**; the next several are all immitis
  at 0.0023–0.0027. Note `Michoacan_2`'s own nearest neighbour is `B3476` — they are a
  closely related immitis pair (but *not* duplicates: 0.00112 is above the 0.001 dedup
  threshold, so the pipeline correctly kept both).
- **Not a hybrid, not contamination**: assembly is the *best* of the pair discussed here
  (257 contigs, N50 274 kb, 27.02 Mb) — a clean, contiguous assembly sitting cleanly
  inside one species' cloud.

**No hybrids anywhere in the dataset.** Across all 529 strains, the ratio
(mean distance to nearest other species ÷ mean distance to own labelled species) has
median 4.54 and 1st percentile 3.83; **zero strains** fall in the 0.9–1.5 window an
intermediate/hybrid genome would occupy. The distribution is strictly bimodal: 527 strains
at ratio ≈ 4–5 and the two mislabels at 0.16 and 0.22.

### 8.2 What should actually change
1. **The error is upstream of the samplesheet, in the source assembly freeze.** The
   study's `species.csv` records the source file as
   `.../annotation_freeze/20260112/pep/Coccidioides_immitis_485B-1_L_OLD_CPA0023.proteins.fa`
   — the species is baked into the **filename** in
   `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation_freeze/20260112/`,
   and likewise `Coccidioides_posadasii_B3476.*`. Every derived config
   (`config.csv`, `config_genus_vs_ureesii.csv`, `config_immitis.csv`,
   `config_posadasii_in_immitis_out.csv`, …) inherits it, so **six config files and the
   `DATA_MANIFEST.yaml` all carry the same wrong label**. Fixing only the samplesheet
   would leave the error to reappear in the next study built from that freeze.
   **Recommend: correct the annotation freeze's metadata (and ideally the filenames), then
   regenerate the study configs.** Flag it to whoever owns the 2025_All_Cocci freeze.
2. **Impact on the results already produced**: the two `IN`-group strains were used in a
   dataset where both species are ingroup, so the presence matrix is unaffected. But they
   were *also* used in `config_immitis_in_posadasii_out.csv` and
   `config_posadasii_in_immitis_out.csv`, where each sits on the **wrong side of the
   ingroup/outgroup boundary**. Any novelty/loss call from those two runs involving the
   species boundary is contaminated by a single misplaced genome — enough to turn a
   genuinely species-specific family into "present in the outgroup". **Those two runs
   should be re-run after the labels are fixed.**
3. **Do not confirm by mash alone if it matters for publication** — mash placement this
   clean is strong evidence, but a confirmatory marker (the standard *Coccidioides*
   diagnostic loci, or a quick BUSCO/marker tree placing the two strains) is cheap
   insurance before amending a public record.
4. **Yes, the pipeline should detect this automatically.** The check in Section 7.3 item 5
   is a few lines, uses the mash matrix `ASSIGN_CLADES` already has in hand, and flagged
   exactly these two with zero false positives. **Where it belongs**: a new
   `label_concordance.tsv` output from `ASSIGN_CLADES` (or a small sibling process in
   `modules/pangenome/mash.nf`), emitting one row per strain with
   `Short, label, nearest_group, d_own, d_nearest, ratio, flagged`, and a loud stderr
   warning per flagged strain. It should **warn, never fail** — a genuinely admixed or
   recently diverged population could trip it legitimately, and the pipeline must not
   refuse to run on a biologically interesting dataset. Grouping key: the samplesheet's
   `Species` column when populated, else `TaxonGroup`.

---

## 9. Recommendations, ranked

### (a) Do now — cheap, high value
1. **Fix the two species labels at the annotation-freeze source**, regenerate the study
   configs, and re-run `config_immitis_in_posadasii_out` / `config_posadasii_in_immitis_out`
   (Section 8.2). Zero compute, direct correctness impact.
2. **Set `mash sketch -s 100000`** in `modules/pangenome/mash.nf` (ideally as a new
   `pangenome_mash_sketch_size` param defaulting to 100000). Measured cost at 529 genomes:
   ~2.5 min on 24 cores, < 1% of study runtime. 323 → 16,092 distinct distance values.
3. **Publish the full mash distance matrix** as `mash_distance.tsv` next to
   `clade_assignments.tsv`. It is already computed and thrown away. It is the continuous
   relatedness metric, and every diagnostic in Section 7.3 needs it.
4. **Add the label-concordance QC** (Section 8.2 item 4) and the **structure diagnostics**
   (PCoA1 variance, between/within ratio) as stderr warnings + a small TSV from
   `ASSIGN_CLADES`. This is what turns a silent failure mode into a visible one.
5. **Guard `ASSIGN_CLADES` against the no-structure case**: if PCoA1 variance < ~40% or the
   winning partition's between/within median ratio < ~1.5, emit a loud warning that the
   clade labels are unreliable and that `min_clades` should not be trusted. Consider
   refusing to write `TaxonGroup` in that case.
6. **Re-run this study with the `RESCUE_PASS` fix.** Everything in Sections 4.3, 5 and 6.3
   is gated on having a presence matrix whose absences are real.

### (b) Worth a ticket
7. **Reorder the DAG so `DEREPLICATE` precedes the per-strain fan-out**, and add an opt-in
   `pangenome_dereplicate_analysis` flag that runs the whole analysis on representatives
   only. Saves ~22% (~69 CPU-h here) at the current threshold. Must be explicit, because
   collapsed strains lose their matrix columns.
8. **Switch dedup from single to average/complete linkage** (or hard-cap the threshold at
   0.002), so raising `pangenome_mash_threshold` cannot chain an entire species into one
   group as 0.005 does today.
9. **`phylo_dispersion_z` per family** (Section 5.2), computed within-species, on a rescued
   matrix, as a new frequency-table column — and use it to replace or supplement the
   binary `min_clades` gate in `pair_classification`. Seconds of compute; the prerequisite
   is the rescue fix, not the statistic.
10. **Make `nf_phyling` an opt-in per-study step** producing a rooted BUSCO marker tree,
    and only then implement Dollo/Sankoff ancestral state reconstruction for gain/loss
    polarisation (`todo/pangenome-phylogeny-aware-gain-loss.md`). Hours-to-a-day of cluster
    time; justified for a study that will publish gain/loss claims, not for every run.
11. **Add the within-dedup-group gene-content Hamming diagnostic** (Section 7.3 item 6) as
    a standing QC metric. It is the cheapest existing detector of a rescue/annotation
    dropout problem, and it would have caught the no-op `RESCUE_PASS` immediately.

### (c) Sounds sophisticated, not worth it here — and why
12. **DAPC.** Its discovery step is the same k-means-on-PCs the pipeline already runs; its
    discriminant projection would make unreproducible clusters *look* well separated. Its
    natural input is SNPs, not a k-mer sketch. See Section 3.3.
13. **Phylogenetically-aware permutation, mixed models, or relatedness-weighted permutation
    for the co-occurrence null.** Measured benefit: ~0.01 in false-positive rate over the
    existing exact 2-strata test (0.050 → 0.037 at 20 strata), and it becomes conservative
    by 80 strata. The existing exact stratified test already removes the 0.327 → 0.050
    inflation. Substantially more code and runtime for no measurable gain. Section 4.2.
14. **Finer sub-species clade labels for this dataset.** They are not reproducible
    (split-half ARI 0.02–0.26) and feeding them to `min_clades` would filter `trans` pairs
    on noise. Section 3.
15. **`mashtree` / any distance-based tree as the basis for gain/loss polarisation.** It
    inherits the same unreproducible within-species topology but presents it as a tree,
    which is worse than having no tree because it invites branch-level claims the data
    cannot support. Section 4.1 / 4.3.
16. **Switching to `skani` / `fastANI` / `sourmash`.** New dependencies and a container
    rebuild to improve a distance estimate that is already at split-half `r = 0.98`. The
    binding constraint is the absence of real structure, not distance accuracy. Reconsider
    only if a study needs a defensible ANI species-boundary number.
17. **A gene-content-based (rather than mash-based) dereplication criterion.** Within-dedup
    -group gene-content variation is currently ~50% explained by assembly contiguity
    differences (Spearman 0.55 with `|Δ n_contigs|`; 48% of differing families unique to the
    smaller assembly). A content-based rule would dereplicate by assembly quality. Revisit
    only after a rescued matrix exists. Section 6.3.
