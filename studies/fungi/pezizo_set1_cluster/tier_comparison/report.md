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


## Analysis 3: real compute cost (wall-hours)

| tier | wall_hours | note |
|---|---|---|
| P | 0.116 | PARTIAL/lower-bound only -- true all-vs-all DIAMOND_SEARCH cost is unmeasurable: storeDir cache hit means Nextflow never logs a trace row for it (unlike ordinary -resume CACHED rows). This total covers only DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, TBLASTN(+MAKEDB), BUILD_PRESENCE_MATRIX -- do NOT read this as Tier P's full cost. |
| C+H | 42.707 | PROFILE_SEARCH:* (gain-side family-profile pathway), measured directly. |
| C+H_loss | 35.394 | PROFILE_LOSS_SEARCH:* (loss-side equivalent), measured directly. |
| C | 0.0 | Free byproduct of Tier C+H's own clustering step. |
| R | 0.032 | Measured directly with `time` around refine_ambiguous_families.py's diamond step -- not read from any trace file (plain subprocess call, no trace row). |


## Analysis 4: decision framework

**Bottom line (revised after review — see note below):** on the curated
16-control set, dropping the HMM step (Tier C) does recover recall relative
to Tier C+H (0.8 vs 0.4 headline; +1, not +2, once HEX1's coincidental hit
is discounted — see Analysis 1 below). But that result does not generalize:
measured genome-wide against Tier P's own candidate lists, Tier C's recall
is only 0.41 (gains) / 0.23 (losses) — it misses the *majority* of what
pairwise search actually finds, while flagging roughly 70% of the entire
ingroup protein universe as candidates (21,258 of 30,461 gain-side proteins
in `presence_matrix.tsv`). A
filter that drops most true positives while barely shrinking the search
space is not a safe triage step. Tier C is not recommended as a pre-filter
ahead of Tier P on this evidence. Tier R adds negligible cost (0.032 wall-h)
on top of Tier C but produced no measurable improvement over it, either on
controls or genome-wide.

*(Revision note: an earlier draft of this section recommended Tier C as a
recall-preserving triage filter ahead of Tier P, reasoning from the
controls-only 0.8 recall figure. An independent review caught that this
contradicts the genome-wide recall-against-P figure already in Analysis 2
above (0.41/0.23) — the 16 controls are a small, hand-picked, favorable
sample, not representative of Tier C's real-world miss rate. The review also
caught an arithmetic error in the controls trade-off sentence and an
unsupported "essentially free" framing for a deployment scenario whose cost
was never isolated. All three are corrected below; the corrected reasoning
is a real course-change from the first draft, not a wording fix.)*

### What the controls tell us (Analysis 1)

Cluster+HMM's headline weakness (recall 0.4) is not one uniform quality
problem — it is several distinct effects that happen to combine into a net
recall gain when the HMM step is dropped, some of them more reliable than
others. Tier C flips HEX1, ADA1, and HAM5 from miss to hit but flips SPA1
from hit to miss, relative to Tier C+H — net +2 (3 gained, 1 lost), giving
Tier C's 4/5 vs. Tier C+H's 2/5 (HAM8, SPA1).

Of those three gains, only ADA1 and HAM5 are genuine fixes: both are
correctly ingroup-only clustered families that Tier C+H's permissive HMM
threshold falsely called present in Mcir/Spom — dropping the HMM step
removes that false signal. **HEX1 is not a genuine fix and should not be
counted as one when weighing Tier C's reliability**: its family also
contains Ncra's IF5A_NEUCR paralog (the mmseqs clustering-level defect Tier
R was built to address), so Tier C's "hit" here is coincidental
outgroup-blindness, not correct detection — Tier C simply cannot see that
the family, as clustered, actually is conserved in every outgroup (Tier
C+H's HMM search correctly detects this and rejects it; see Task 5's Finding
2b). SPA1, the one loss, is a real cost: mmseqs only partially clustered it
(2 of 5 ingroup orthologs merged), and Tier C+H's HMM step compensated by
finding the 3 ingroup members clustering missed — an HMM presence step can
add true signal by rescuing under-clustered orthologs, not only introduce
false signal by over-calling conserved ones. Net effect, discounting HEX1's
coincidental hit: Tier C's *reliable* gain over Tier C+H on this control set
is ADA1 + HAM5 − SPA1 = +1, not +2 — a real but modest improvement, achieved
by a tier that (per Analysis 2 below) is unreliable at genome scale.

Tier R's HEX1 "hit" is inherited unmodified from Tier C, not a genuine
paralog-split rescue — its family didn't actually split (the HEX1/eIF-5A
pair remains connected in the component graph through shared edges to other
family members beyond the one forced-cut edge; see Task 5's investigated
finding). Tier R getting credit for HEX1 here is coincidental twice over —
neither Tier C's nor Tier R's HEX1 "hit" demonstrates what either tier was
supposed to fix.

A structural blind spot neither Tier C nor Tier R can address: **POS_LAH is
`unresolved` in C+H, C, and R alike** (a real Tier P hit that never joined
any profiled mmseqs family — spec failure mode 3, clustering-sensitivity
exclusion). Any presence rule built on top of cluster membership, refined or
not, can never surface a protein clustering excluded entirely — exactly the
fast-evolving, lineage-specific class of gene this pipeline exists to find.
This is a genuine, evidenced limitation of the whole cluster-membership
family of tiers (C and R), independent of the HMM-vs-no-HMM question.

The second clade (Agaricales) confirms C+H/C/R all correctly resolve both
positive controls — but with only 2 positive/0 negative controls it cannot
differentiate the tiers from each other. Its real finding is methodological:
Tier P's own controls-scoring path structurally cannot score `fasta`
-anchored controls at all, so its recall there is permanently `n/a`, not a
biological 0. (Reproducing this section also requires `hmmer` on PATH — see
open items below.)

### What genome-wide concordance tells us (Analysis 2)

Tier C+H is conservative relative to Tier P (1,044 gain candidates,
precision 0.47, recall 0.15). Tier C is the opposite: 21,258 gain candidates
(20x more; roughly 70% of the 30,461-protein ingroup universe the presence
matrix covers), recall rises to 0.41 (2.8x, not quite triple), but precision
collapses to 0.066 — roughly 14 of every 15 Tier-C "novel" calls are not
confirmed by Tier P. Losses show the same pattern, more extreme (precision
0.42 → 0.0029, a ~144x drop). This is the genome-scale expression of the
same HMM-detection-overreach mechanism found in Analysis 1: without any step
that actually screens outgroups, raw clustering cannot distinguish a truly
lineage-restricted family from one that merely failed to cluster with its
real, conserved outgroup orthologs.

**This recall figure is the one that matters for any triage-filter
proposal, and it rules one out.** Tier C's 0.41/0.23 recall against Tier P
means it *misses 59% of Tier P's gain candidates and 77% of its loss
candidates outright* — a candidate a triage filter fails to flag never
reaches pairwise confirmation at all. Combined with Tier C flagging ~70% of
the protein universe, escalating every Tier C candidate to Tier P would cut
Tier P's own search space by at most ~30% (and Tier P's true search cost is
itself unmeasured — see Analysis 3) while still silently dropping most of
what full pairwise search would have found. That is not what a viable triage
filter looks like: a triage step must preserve recall while shrinking the
search space, and Tier C does neither strongly enough to justify the risk.
The spec's own required per-miss classification (clustering-level vs.
HMM-level, for each genome-wide miss) was not carried out in this
investigation — only the aggregate precision/recall/Jaccard numbers above
were computed — so it isn't yet known how much of Tier C's genome-wide
recall gap is LAH-style clustering exclusion vs. a different mechanism; that
gap should be closed before any triage-filter proposal is reconsidered.

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

Tier C+H's fully-measured cost — 42.7 wall-h gain-side + 35.4 wall-h loss-side,
~78 wall-h total — is the reliable baseline here. Tier C is a free byproduct
of that same run (0 wall-h marginal cost) and Tier R adds essentially nothing
on top (0.032 wall-h, under two minutes). Tier P's true all-vs-all diamond
cost could not be measured at all (its results were `storeDir`-cached before
any trace history this investigation has access to begins) — the 0.116
wall-h reported is a small fragment (self-search, hit-parsing, tblastn) that
excludes the O(N²) search this whole investigation exists to avoid, so the
central "how much compute do we save" question only has an indirect answer.
Tier C+H's own 42.7/35.4 wall-h totals also bundle mmseqs clustering together
with the HMM search step (`PROFILE_SEARCH:*` was measured as one prefix,
not broken into its `MMSEQS_FAMILY_CLUSTER` vs. `FAMILY_HMMSEARCH:*`
sub-costs) — Tier C's *own* marginal cost, if run without ever building the
family HMMs at all, is very likely much smaller than "free once you've
already paid for C+H," but this investigation did not isolate that number.

### Recommendation

1. **Do not deploy Tier C (or R) as a replacement for Tier P or Tier C+H's
   final novelty/loss calls.** Genome-wide precision is far too low to trust
   directly (0.066 gains, 0.0029 losses) — the overwhelming majority of Tier
   C's extra candidates are not confirmed by the pairwise reference.
2. **Do not deploy Tier C as a recall-preserving pre-filter ahead of Tier P
   either, on the evidence measured here.** The controls-only 0.8 recall
   figure does not generalize: genome-wide, Tier C's recall against Tier
   P's own candidates is only 0.41 (gains) / 0.23 (losses) — it would
   silently drop 59-77% of what full pairwise search finds — while flagging
   ~70% of the protein universe as candidates, so escalating everything it
   flags saves at most ~30% of Tier P's (itself unmeasured) search cost.
   That combination — high miss rate, small space reduction — does not meet
   the bar for a viable triage step. If a genuine prioritization use is
   wanted, it should be scoped narrower than "escalate every Tier C
   candidate": e.g. candidates near the `other_max_frac` boundary in Tier
   C+H's own output (small in number, plausibly enriched for true misses),
   not Tier C's full, low-precision candidate set. That narrower use was not
   itself measured in this investigation and would need its own validation
   before being recommended.
3. **Tier R is not yet earning its keep on this dataset** — it neither
   changes genome-wide precision/recall relative to Tier C, nor
   demonstrably fixes the one control case (HEX1) it was designed around
   (Tier R's HEX1 "hit" is inherited unmodified from Tier C's own
   coincidental outgroup-blindness, not a genuine split — see Analysis 1).
   Hold it back from the promotion path (per the spec's own gate) until a
   dataset with more, and more clearly ambiguous, paralog-merge cases can
   test its graph-splitting mechanism under real load — or until the
   MCL-based fallback (noted but not built here) is tried against the same
   HEX1 case, since the simple forced-single-edge-cut approach demonstrably
   failed to separate a family with more than one qualifying inter-member
   edge.
4. **What would change recommendation 2** (the one actually being made —
   not to use Tier C as a pre-filter): this recommendation should be
   revisited if a future dataset, or a differently-tuned Tier C (e.g. a
   looser `--ingroup-min-frac`), shows genome-wide recall against Tier P
   climbing well above the 0.41/0.23 measured here — a level around 0.7-0.8
   is the rough bar at which "misses most of what pairwise finds" would stop
   being true, making a pre-filter role defensible. Separately, Tier R's
   status (recommendation 3) should be revisited if the fraction of
   "ambiguous" families (species-duplicated *and* near-miss) roughly doubles
   from the ~17% measured here, to somewhere near a third of all families —
   a heuristic bar, not a derived one (this investigation has only one
   clade's data point), chosen because Tier R's effect was already too
   small to measure (Δ under 0.3%) at 17%, so "still too small to matter"
   stops being a safe assumption to carry forward unexamined once the
   ambiguous population roughly doubles in size.
5. **Known open items for anyone extending this**: (a) the BUSCO negative
   -control identifier-namespace mismatch (Task 5's Finding 3) leaves
   `fp_rate` unmeasured for pezizo_set1 throughout — every `fp_rate` cell in
   Analysis 1's table above is blank, not a verified 0.0; (b) Tier P's true
   search cost remains unmeasured (Analysis 3); (c) Tier C's own marginal
   cost independent of Tier C+H's HMM step is not isolated in the available
   trace data; (d) the spec's required per-miss clustering-vs-HMM
   classification for Analysis 2's genome-wide misses was not carried out —
   only aggregate precision/recall/Jaccard were computed; (e) reproducing
   the Agaricales section requires `hmmer` on `PATH` (fixed for this repo's
   own `pixi` environment in commit 297b763, but `score_controls.py` itself
   still swallows a missing-`hmmsearch` `FileNotFoundError` into a generic
   "no hit" note rather than a distinguishable error — a latent trap for any
   environment that doesn't have this repo's `pixi.toml` fix).
