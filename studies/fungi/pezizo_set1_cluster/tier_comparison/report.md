# Cluster-vs-pairwise sensitivity: pezizo_set1

## Analysis 1: controls recall/FP by tier

| tier | n_controls | n_positive | n_positive_resolved | positive_hits | recall | n_negative | n_negative_resolved | negative_fp | fp_rate | n_unresolved | n_placeholder_skipped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 16 | 6 | 6 | 6 | 1.0 | 10 | 0 | 0 |  | 10 | 1 |
| C+H | 16 | 6 | 5 | 2 | 0.4 | 10 | 0 | 0 |  | 11 | 1 |
| C | 16 | 6 | 5 | 4 | 0.8 | 10 | 0 | 0 |  | 11 | 1 |
| R | 16 | 6 | 5 | 4 | 0.8 | 10 | 0 | 0 |  | 11 | 1 |


### Per-control outcomes

| control_id | class | expected_call | anchor_type | resolved_family | actual_call | outcome | note | tier |
|---|---|---|---|---|---|---|---|---|
| POS_HEX1 | positive | novel | protein_id | sp|P87252|HEX1_NEUCR | novel | hit |  | P |
| POS_LAH | positive | novel | protein_id | tr|V5IRA6|V5IRA6_NEUCR | novel | hit |  | P |
| POS_ADA1 | positive | novel | protein_id | tr|Q7SE74|Q7SE74_NEUCR | novel | hit |  | P |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | novel | hit |  | P |
| POS_HAM8 | positive | novel | protein_id | tr|Q7SEK3|Q7SEK3_NEUCR | novel | hit |  | P |
| POS_SPA1 | positive | novel | protein_id | tr|Q7SI25|Q7SI25_NEUCR | novel | hit |  | P |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| POS_HEX1 | positive | novel | protein_id | tr|A0ACF5BXR3|A0ACF5BXR3_9PEZI | not_novel | miss |  | C+H |
| POS_LAH | positive | novel | protein_id |  | unresolved | unresolved | protein not in any profiled family | C+H |
| POS_ADA1 | positive | novel | protein_id | tr|A0ACF5CBF7|A0ACF5CBF7_9PEZI | not_novel | miss |  | C+H |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | not_novel | miss |  | C+H |
| POS_HAM8 | positive | novel | protein_id | tr|J3K3I3|J3K3I3_COCIM | novel | hit |  | C+H |
| POS_SPA1 | positive | novel | protein_id | tr|J3K8R1|J3K8R1_COCIM | novel | hit |  | C+H |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| POS_HEX1 | positive | novel | protein_id | tr|A0ACF5BXR3|A0ACF5BXR3_9PEZI | novel | hit |  | C |
| POS_LAH | positive | novel | protein_id |  | unresolved | unresolved | protein not in any profiled family | C |
| POS_ADA1 | positive | novel | protein_id | tr|A0ACF5CBF7|A0ACF5CBF7_9PEZI | novel | hit |  | C |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | novel | hit |  | C |
| POS_HAM8 | positive | novel | protein_id | tr|J3K3I3|J3K3I3_COCIM | novel | hit |  | C |
| POS_SPA1 | positive | novel | protein_id | tr|J3K8R1|J3K8R1_COCIM | not_novel | miss |  | C |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| POS_HEX1 | positive | novel | protein_id | tr|A0ACF5BXR3|A0ACF5BXR3_9PEZI | novel | hit |  | R |
| POS_LAH | positive | novel | protein_id |  | unresolved | unresolved | protein not in any profiled family | R |
| POS_ADA1 | positive | novel | protein_id | tr|A0ACF5CBF7|A0ACF5CBF7_9PEZI | novel | hit |  | R |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | novel | hit |  | R |
| POS_HAM8 | positive | novel | protein_id | tr|J3K3I3|J3K3I3_COCIM | novel | hit |  | R |
| POS_SPA1 | positive | novel | protein_id | tr|J3K8R1|J3K8R1_COCIM | not_novel | miss |  | R |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |


### Second clade: Agaricales (spc14/spc33)

Independent ingroup/outgroup split (`agaricomycetes_pairwise`/`agaricomycetes_mmseqs`), 2 positive controls, 0 negatives — run to check the pezizo_set1 findings above aren't a one-clade artifact.

| tier | n_controls | n_positive | n_positive_resolved | positive_hits | recall | n_negative | n_unresolved |
|---|---|---|---|---|---|---|---|
| P | 2 | 2 | 0 | 0 | n/a | 0 | 2 |
| C+H | 2 | 2 | 2 | 2 | 1.0 | 0 | 0 |
| C | 2 | 2 | 2 | 2 | 1.0 | 0 | 0 |
| R | 2 | 2 | 2 | 2 | 1.0 | 0 | 0 |

Both controls are `fasta`-anchored. Tier P's `n/a` here is not a biological miss — `score_controls.py`'s pairwise mode has no family HMM database to `hmmsearch` a `fasta` anchor against by design (its own module docstring states this), so it structurally cannot score a `fasta`-anchored control at all. With only 2 positive controls and 0 negatives, this set cannot differentiate C+H/C/R from each other (all three tie at 1.0) — its value here is confirming there's no second-clade surprise, not adding statistical power.

## Analysis 2: genome-wide concordance (gains)

| baseline_count_vs_C | delta_pct_vs_C | delta_vs_C | jaccard_vs_P | n_candidates | precision_vs_P | recall_vs_P | refined_count_vs_C | tier |
|---|---|---|---|---|---|---|---|---|
| 21258 | -95.1 | -20214 | 0.125 | 1044 | 0.4722 | 0.1453 | 1044 | C+H |
|  |  |  | 0.06 | 21258 | 0.0656 | 0.4111 |  | C |
| 21258 | -0.0 | -7 | 0.06 | 21251 | 0.0656 | 0.4111 | 21251 | R |


## Analysis 2: genome-wide concordance (losses)

| baseline_count_vs_C | delta_pct_vs_C | delta_vs_C | jaccard_vs_P | n_candidates | precision_vs_P | recall_vs_P | refined_count_vs_C | tier |
|---|---|---|---|---|---|---|---|---|
| 9901 | -99.8 | -9877 | 0.0709 | 24 | 0.4167 | 0.0787 | 24 | C+H |
|  |  |  | 0.0029 | 9901 | 0.0029 | 0.2283 |  | C |
| 9901 | -0.2 | -15 | 0.0029 | 9886 | 0.0029 | 0.2283 | 9886 | R |


## Analysis 3: real compute cost (CPU-hours)

| tier | cpu_hours | note |
|---|---|---|
| P | 0.116 | PARTIAL/lower-bound only -- true all-vs-all DIAMOND_SEARCH cost is unmeasurable: storeDir cache hit means Nextflow never logs a trace row for it (unlike ordinary -resume CACHED rows). This total covers only DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, TBLASTN(+MAKEDB), BUILD_PRESENCE_MATRIX -- do NOT read this as Tier P's full cost. |
| C+H | 42.707 | PROFILE_SEARCH:* (gain-side family-profile pathway), measured directly. |
| C+H_loss | 35.394 | PROFILE_LOSS_SEARCH:* (loss-side equivalent), measured directly. |
| C | 0.0 | Free byproduct of Tier C+H's own clustering step. |
| R | 0.032 | Measured directly with `time` around refine_ambiguous_families.py's diamond step -- not read from any trace file (plain subprocess call, no trace row). |


## Analysis 4: decision framework

**Bottom line:** raw mmseqs cluster membership (Tier C) is a real, near-zero
-cost proxy that *beats* the family-HMM pipeline's own recall on curated
controls (0.8 vs 0.4), but at genome scale it trades that recall gain for a
severe loss of precision (0.066 vs 0.47 on gains; 0.0029 vs 0.42 on losses).
It is not a drop-in replacement for either Tier P or Tier C+H — but it is a
strong, essentially-free first-pass triage filter ahead of a slower,
specific pass. Tier R adds negligible cost (0.032 CPU-h) on top of Tier C
but produced no measurable improvement over Tier C on this dataset, either
on the curated controls or genome-wide — its case rests entirely on future
data where its target failure mode is more common than it was here.

### What the controls tell us (Analysis 1)

Cluster+HMM's headline weakness (recall 0.4) is not one uniform quality
problem — it is two opposite-direction failure modes that happen to net out
close to even. Dropping the HMM step (Tier C) recovers ADA1 and HAM5 (both
correctly ingroup-only clustered families the HMM's permissive threshold
falsely called present in Mcir/Spom) but *loses* SPA1, a family mmseqs only
partially clustered (2 of 5 ingroup orthologs merged) that the HMM step's
own sequence search compensated for by finding the 3 ingroup members
clustering missed. Net, Tier C nets +2 controls (4/5 vs 2/5) — a real gain,
but a trade, not a strict improvement: an HMM presence step can add true
signal by rescuing under-clustered orthologs, not only introduce false
signal by over-calling conserved ones.

HEX1 reads `hit` under Tier C and Tier R alike, but Tier R's `hit` is
inherited unmodified from Tier C, not a genuine paralog-split rescue — its
family didn't actually split (the HEX1/eIF-5A pair shares independent
qualifying diamond hits beyond the one forced-cut edge, so the connected
component stayed whole; see Task 5's investigated finding). Tier R getting
credit for HEX1 here is coincidental, not evidence the refinement mechanism
works as designed.

The second clade (Agaricales) confirms C+H/C/R all correctly resolve both
positive controls — but with only 2 positive/0 negative controls it cannot
differentiate the tiers from each other. Its real finding is methodological:
Tier P's own controls-scoring path structurally cannot score `fasta`
-anchored controls at all, so its recall there is permanently `n/a`, not a
biological 0.

### What genome-wide concordance tells us (Analysis 2)

Tier C+H is conservative relative to Tier P (1,044 gain candidates,
precision 0.47, recall 0.15). Tier C is the opposite: 21,258 gain candidates
(20x more), recall triples to 0.41, but precision collapses to 0.066 —
roughly 14 of every 15 Tier-C "novel" calls are not confirmed by Tier P.
Losses show the same pattern, more extreme (precision 0.42 → 0.003, a ~140x
drop). This is the genome-scale expression of the same HMM-detection
-overreach mechanism found in Analysis 1: without any step that actually
screens outgroups, raw clustering cannot distinguish a truly lineage
-restricted family from one that merely failed to cluster with its real,
conserved outgroup orthologs.

Tier R is statistically indistinguishable from Tier C at genome scale
(Δ = -7 candidates / -0.0% gains, -15 / -0.2% losses) — refinement
essentially never fires with consequence on this dataset, consistent with
Analysis 1: only 1,227 of 7,281 families (~17%) were ever flagged ambiguous,
and few enough of those actually split (or split without changing their
novelty call) to move the genome-wide count by more than a fraction of a
percent. No candidate-count *inflation* was observed (the over-splitting
risk flagged as failure mode 4 in the spec) — if anything, refinement very
slightly shrinks the candidate count in both directions. This dataset gives
no evidence the refinement approach is unsafe; it also gives very little
evidence it does much of anything.

### What the real cost data tells us (Analysis 3)

Tier C+H's fully-measured cost — 42.7 CPU-h gain-side + 35.4 CPU-h loss-side,
~78 CPU-h total — is the reliable baseline here. Tier C is a free byproduct
of that same run (0 CPU-h marginal cost) and Tier R adds essentially nothing
on top (0.032 CPU-h, under two minutes). Tier P's true all-vs-all diamond
cost could not be measured at all (its results were `storeDir`-cached before
any trace history this investigation has access to begins) — the 0.116
CPU-h reported is a small fragment (self-search, hit-parsing, tblastn) that
excludes the O(N²) search this whole investigation exists to avoid, so the
central "how much compute do we save" question only has an indirect answer.
Tier C+H's own 42.7/35.4 CPU-h totals also bundle mmseqs clustering together
with the HMM search step (`PROFILE_SEARCH:*` was measured as one prefix,
not broken into its `MMSEQS_FAMILY_CLUSTER` vs. `FAMILY_HMMSEARCH:*`
sub-costs) — Tier C's *own* marginal cost, if run without ever building the
family HMMs at all, is very likely much smaller than "free once you've
already paid for C+H," but this investigation did not isolate that number.

### Recommendation

1. **Do not deploy Tier C (or R) as a replacement for Tier P or Tier C+H's
   final novelty/loss calls.** Genome-wide precision (0.03-0.07) is too low
   to trust directly — most of Tier C's extra candidates are not confirmed
   by the pairwise reference.
2. **Do use Tier C as a cheap, high-recall first-pass triage filter feeding
   into Tier P**, per the framework the spec proposed: run mmseqs clustering,
   call Tier C's raw-membership presence genome-wide, and escalate every
   Tier-C-flagged candidate — plus any Tier C+H near-miss sitting close to
   the `other_max_frac` boundary — to full pairwise confirmation. Tier C's
   higher recall than Tier C+H on the curated controls (0.8 vs 0.4) makes
   this ordering less likely to silently drop a true positive than using
   Tier C+H alone as the triage step; the risk this framework accepts is
   Tier-C-style false positives, not Tier-C+H-style false negatives, and
   every Tier C candidate still gets a full pairwise check downstream.
3. **Tier R is not yet earning its keep on this dataset** — it neither
   changes genome-wide precision/recall relative to Tier C, nor
   demonstrably fixes the one control case (HEX1) it was designed around.
   Hold it back from the promotion path (per the spec's own gate) until a
   dataset with more, and more clearly ambiguous, paralog-merge cases can
   test its graph-splitting mechanism under real load — or until the
   MCL-based fallback (noted but not built here) is tried against the same
   HEX1 case, since the simple forced-single-edge-cut approach demonstrably
   failed to separate a family with more than one qualifying inter-member
   edge.
4. **What would invalidate this framework**: if a future dataset shows
   HEX1-style clustering-inherent paralog inflation is common (not the
   apparent near-singleton case here), Tier C's recall advantage over Tier
   C+H would shrink or reverse, since neither tier's presence rule fixes
   that failure mode — only a clustering-time or refinement-time fix does,
   and this investigation found no evidence the refinement approach as
   built handles it. Track the fraction of "ambiguous" families
   (species-duplicated *and* near-miss) in future studies — it was ~17%
   here; if it climbs substantially elsewhere, revisit.
5. **Known open items for anyone extending this**: (a) the BUSCO negative
   -control identifier-namespace mismatch (Task 5's Finding 3) leaves
   `fp_rate` unmeasured for pezizo_set1 throughout — every `fp_rate` cell in
   Analysis 1's table above is blank, not a verified 0.0; (b) Tier P's true
   search cost remains unmeasured (Analysis 3); (c) Tier C's own marginal
   cost independent of Tier C+H's HMM step is not isolated in the available
   trace data.
