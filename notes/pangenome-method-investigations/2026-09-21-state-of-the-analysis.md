# Where the pangenome methods analysis stands — 2026-09-21

A summary of what four investigations established, what changed as a result, and what is
still open. Companion to the detailed notes in this directory.

## The short version

The pangenome subworkflow works. Its **inputs** are the problem, and the two largest
input defects were both invisible: a silently-failed rescue pass, and an assembly-quality
confound running through the accessory genome. Neither is fixable by tuning the parameter
that causes them.

## 1. What was confirmed working

- **`LEIDEN_MODULES`** reproduces the study's hand-run **byte-identically**
  (`md5 7bb3903422df1779fb6d3be3ee88c366`); a full pipeline re-run gives ARI 0.926 and
  97.1% identical module labels.
- **Species-level phylogrouping** is essentially perfect: mash -> PCoA -> k-means gives
  **527/529 concordance** with species labels, and the two exceptions are mislabels, not
  biology.
- **The co-occurrence null is correctly calibrated.** Measured FPR at alpha=0.05:
  unstratified Fisher **0.327**, the shipped k=2 stratified test **0.050**, k=10/20/40/80
  all 0.033-0.040. The existing stratification is right and does not need replacing.

## 2. What was found broken

### RESCUE_PASS applied zero cells on the flagship run

`presence_matrix.tsv` and `presence_matrix.rescued.tsv` were byte-identical; zero
`genome_only` cells. Cause: 530 x `zstd -dc` refusing Nextflow's staged symlinks, with the
exit code ignored. **301.6 task-hours of TBLASTN produced nothing**, and every downstream
table was computed on unrescued data — and the run was reported, reviewed and treated as
valid.

Fixed (`e086c78`), plus a guard (#127) so it can never be silent again.

### But the rescue should not simply be switched on

Exhaustive analysis of all 12,877,646 rescuable cells:

| where the hit lands | share |
|---|---|
| on a predicted gene of a **different family** | **77.6%** |
| intergenic | 22.4% |
| same family | 0% |

Split: redundant-locus **74.2%** | true paralog 3.3% | dubious short/repeat 21.2% |
**plausible genuine annotation dropout ~1.3%**.

Enabling it as configured would take density 15.63% -> 60.28%, three-quarters of it the
same locus counted twice. **Thresholds are the wrong lever**: the class mix is flat at
19-23% intergenic from (90,80) all the way to (99.9,99), because these hits are
100%-identical *real DNA*.

### Assembly fragmentation inflates the accessory genome

`rho(accessory ~ N50)` = **-0.53**, `rho(accessory ~ contigs)` = **+0.53**, after
controlling for assembly size. Within *C. posadasii*, the 15 most fragmented assemblies
carry **283 more accessory families (+11%)** than the 15 best, on slightly **less** DNA.

Mechanism, confirmed directly: proteins in strain-private families sit at a contig
terminus **31.4%** of the time, against **9.38%** of all 4.5M proteins — **3.3x
enrichment**. A gene broken by a contig boundary yields two partial models that fail
clustering and become private families.

### Gain/loss direction is currently an artifact

Coccidioides: **99.3% gain, 0.7% loss, 0% ambiguous** (143:1).
A. fumigatus, where rescue worked (5,871,769 cells applied): **77.5% / 10.9% / 11.6%**
(7.1:1).

The method is sound. Coccidioides is the broken case, and the difference is the rescue.

## 3. What the tier-1 sweep settled — a negative result

12 clusterings across `min_id` x `cov` plus coverage-mode variants and a post-cluster merge.

**Eleven of twelve leave `rho(accessory~N50)` in a 0.017 range (-0.511 to -0.528)** while
family count moves 2x (30,406 to 60,138). The confound rides straight through clustering.

`--cov-mode 1` was the sole exception (-0.519 -> -0.350, singletons -44%) and is
**poisoned**: of the 168,920 pairs it newly merges, **76.0% have neither member at a
contig end** and 28.9% join representatives with disjoint Pfam architectures. It improves
the metric by collapsing intact genes, not by repairing broken ones.

**Conclusion: keep 0.9 / 0.8 / cov-mode 0. Evidence class: empirical, negative.** Not
"this value is optimal" — nothing distinguishes it from 0.8/0.8 — but **"this parameter
is not the lever."**

This reorders the plan: **#133 (structural rescue criterion) becomes the primary fix**,
having been sequenced as a fallback behind #132.

## 4. What generalizes to other pangenome datasets

- **A single-species study gets a confident-looking `TaxonGroup` of pure noise.**
  Split-half ARI is 1.000 for species and **0.02-0.26** for every within-species
  partition. `min_clades >= 2` then filters on that noise silently.
- **Trust a clade grouping only if** PCoA1 > ~50% of distance variance, k=2-3 survives
  split-half ARI > 0.8, and between/within median distance ratio > ~2. Otherwise set
  `pangenome_assign_clades = false` and treat the study as unstratified.
- **Sketch at `-s 100000`.** The default s=1000 yields only **323 distinct distance
  values** across 139,656 pairs. Costs 2.5 min at 529 genomes.
- **Sub-species structure does not exist here** — two real groups plus a continuum. DAPC
  would not help; its discovery step is the same k-means-on-PCs and its discriminant
  projection would make noise look separated. The A. fumigatus 4-8 clusters almost
  certainly came from variant data, not minhash.
- **Do not build** a phylogenetically-aware or mixed-model co-occurrence null. All the
  confounding is on the species axis and the k=2 stratification already removes it.

## 5. Where it leaves the flagship study

**Every number from `genus_vs_ureesii` is computed on an unrescued matrix at 15.6%
density.** The DNA-only findings (relatedness, clade structure, the two mislabels) are
unaffected. Everything that depends on presence/absence — islands, Leiden modules,
gain/loss, the report tables — needs re-measuring after #133 lands and the study is re-run.

That re-run is the single highest-value action outstanding, and it is also the experiment
that tests this whole analysis: **if gain:loss moves from 143:1 toward 7:1 with a non-zero
ambiguous fraction, the absence-calling explanation is confirmed.** If it does not,
something else is wrong.

## 6. Two strains to fix

| Short | labelled | clusters with | evidence |
|---|---|---|---|
| `485B-1_L_OLD_CPA0023` | *C. immitis* | posadasii | 0.000072 from a posadasii sibling in the same isolate series |
| `B3476` | *C. posadasii* | immitis | sits centrally in the immitis cloud, clean assembly |

Across all 529 strains **zero** fall in the 0.9-1.5 hybrid window; only these two, at 0.16
and 0.22. Mislabels, not hybrids or contamination. **The error is upstream** — the species
is baked into the annotation-freeze filenames and propagates into six config CSVs. Both
sit on the wrong side of the IN/OUT boundary in `config_immitis_in_posadasii_out.csv` and
`config_posadasii_in_immitis_out.csv`, so those two runs are invalid as labelled.

## 7. Method note worth carrying forward

`--cov-mode 1` was recommended, mid-session, on the strength of a 33% confound reduction.
The paralog check — requested only because the issue template demanded a stated risk —
is what revealed that 76% of its merges were intact interior genes. **The headline metric
improved for the wrong reason.**

Every remaining sweep in #132 (core/shell cutoffs, `pair_class_k`, Leiden resolution)
needs an equivalent "what could this break" measure alongside the metric being optimised.

## Open issues

#130 assembly-quality QC · #131 assumptions register · #132 sensitivity analysis
(priority 1 complete and negative; 2-4 open) · #133 structural rescue criterion ·
#134 diagnostics surfacing with `--pangenome_strict`
