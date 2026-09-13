# Cluster-vs-pairwise sensitivity investigation: what does mmseqs clustering cost us, and can it prioritize candidates?

Date: 2026-09-13

## Problem

`nf_NovInvenio`'s pairwise pathway (`--cluster_tool pairwise`, the default:
diamond all-vs-all + tblastn re-search + per-genome self-hit paralog-cutoff
filtering) is the pipeline's sensitivity reference for novelty (ingroup-only
gain) and loss (ingroup-only absence) candidate discovery. It is also slow —
thousands of diamond/tblastn jobs for a modest 11-genome set. The
`--cluster_tool mmseqs` pathway (mmseqs clustering → per-family HMM profile →
`hmmsearch` against every genome for presence) is materially faster and is
already running as `pezizo_set1_cluster`, built from the same `config.csv`/
`data_dir` as the pairwise `pezizo_set1` run so the two are directly
comparable.

The open question: how much sensitivity does the cheaper pathway lose, where
exactly does it lose it, and — going one step further than what either
pathway currently offers — would raw mmseqs cluster membership *alone*, with
no HMM step at all, be fast enough and accurate enough to serve as a
first-pass prioritization filter ahead of the expensive pairwise pass, for
both directions (novel gains and lineage-specific losses)? And, since one of
the root-caused failure modes below (HEX1) is specifically mmseqs merging
true orthologs with a conserved paralog — could a *targeted* pairwise
refinement step, scoped only to the small subset of families where that risk
actually shows up, recover that failure mode without reintroducing the
`O(N²)` genome-wide cost the whole investigation is trying to avoid, and
without over-splitting families into spurious false-novelty fragments?

## Ground truth already on disk (not assumptions)

Both study runs already exist and controls-scoring machinery
(`nf_NovInvenio/bin/score_controls.py`) already ran on `pezizo_set1`'s 16
curated controls (`configs/controls/pezizo_set1.controls.csv`: 6 positive
lineage-specific Pezizomycotina genes — hex-1, lah, ada-1, ham-5, ham-8, spa-1
— and 10 BUSCO single-copy negatives):

| | pairwise (`pezizo_set1`) | cluster+HMM (`pezizo_set1_cluster`) |
|---|---|---|
| recall (positives) | 1.0 (6/6) | 0.4 (2/5 resolved; 1 unresolved) |
| fp_rate (negatives) | 0.0 | 0.0 |

**Run completeness was verified, not assumed**: the published cluster-run
outputs (`family_profiles.hmm`, `presence_matrix.tsv`, `loss_presence_matrix.tsv`,
mtimes 2026-09-10) come from the 2026-09-10 10:50 nextflow run, in which the
4 initially-FAILED `PROFILE_SEARCH:BUILD_FAMILY_PROFILES:BUILD_CHUNK` tasks
were each retried and completed before those outputs were written — confirmed
by tracing every occurrence of those 4 chunk names in that run's trace file.
The recall=0.4 number reflects a fully-built family set, not a partial one.
(Separately: the 2026-09-12 20:18 resume trace shows 101/123
`PROFILE_LOSS_SEARCH` chunks failing outright even after retry — a live
pipeline-reliability bug, but it never overwrote the published results and is
out of scope for this investigation.)

### Root-causing the 4 misses (done, not proposed)

Traced each control's anchor protein through `families/families_cluster.tsv`
(raw mmseqs membership) and `presence_matrix.tsv` (post-HMM presence calls).
The 4 misses split into **three distinct failure modes**, which matters
because they imply different fixes and different answers to "would a
cluster-only tier do better or worse":

1. **Clustering-level paralog inflation (HEX1).** mmseqs merged the hex-1
   anchor into a 6-member family that also contains `IF5A_NEUCR` (eIF-5A, a
   universally-conserved translation factor) plus outgroup orthologs of
   *that* — the family is present in all 11 genomes at the raw
   cluster-membership level, before any HMM exists. This is exactly the
   paralog pairing the pairwise pipeline's per-genome self-hit
   paralog-cutoff step exists to exclude (confirmed by the controls CSV's
   own note on the fix that made HEX1 usable in pairwise mode). **A
   cluster-only tier would not fix this** — the damage happens at family
   construction, not at presence-calling.

2. **HMM-detection overreach (ADA1, HAM5).** Both anchors' mmseqs clusters
   are *correct*: exactly the 5 ingroup species (Amega, Ncra, Afum, Ztri,
   Cimm), zero outgroup members. But `presence_matrix.tsv` — built by
   `hmmsearch`-ing each family's HMM profile against every genome at
   `--hmm_presence_cov 0.3 --hmm_presence_min_residues 100` — calls "present"
   in Mcir (both) and Spom (ADA1), which alone breaks the default
   `other_max_frac <= 0.0` novelty predicate. The controls CSV's own pairwise
   notes already flagged these exact two genomes as *marginal,
   non-significant* diamond hits (e=1.8e-1, 5.7e-2) — the permissive family
   HMM converts that same weak signal into a confident "present." **A
   cluster-only tier (raw membership, no HMM) would call both of these
   correctly.**

3. **Clustering-sensitivity exclusion (LAH).** The anchor never joined any
   profiled family — a singleton below the minimum family-size floor, most
   likely because a fast-evolving/short lineage-specific sequence is
   precisely the kind of sequence mmseqs's identity/coverage clustering
   threshold struggles to place. This is a distinct failure mode from (1)
   and (2): no presence-calling method fixes it; only looser clustering
   sensitivity (or a floor rescue mechanism) would.

**Implication worth testing at scale**: on this 4-miss sample, a genuinely
cheaper, HMM-free tier would have recovered 2 of 4 misses (ADA1, HAM5) that
the current "faster" cluster+HMM pathway gets wrong — i.e. the fastest tier
may out-perform the current intermediate tier on recall, not just on speed.
That is the central hypothesis this investigation tests at scale rather than
on a 4-control anecdote.

### A fourth failure mode this investigation must also guard against

Failure mode 1 (HEX1) suggests an obvious next move: run a pairwise
diamond/phmmer search *within* a suspect family's small member set to split
true orthologs from a merged-in ancient paralog — essentially OrthoFinder's
own strategy (all-vs-all → normalized score graph → graph clustering, with
tree-based splitting on top) applied only to the family, not the genome.
That is worth testing (see Tier R below), but it is not free of risk:

4. **Over-splitting / fragmentation-induced false novelty (new, not yet
   observed — a risk to test for, not a confirmed failure).** A stricter
   within-family cutoff can cut too deep: it can fragment a real, divergent
   ortholog into its own singleton (worsening failure mode 3's problem for
   genes like LAH), or split a genuinely broadly-conserved family into an
   ingroup-only shard and an outgroup-only shard that then each look "novel"
   on their own, inflating false-novelty calls that the 16-control set
   cannot detect (the 10 BUSCO negatives are all tight single-copy conserved
   genes, exactly the kind least likely to fragment — a clean `fp_rate=0.0`
   on refinement would not prove refinement is safe). This is why Tier R
   (below) is scored against genome-wide candidate counts (Analysis 2), not
   just the curated controls.

## Four tiers under comparison

- **Tier P (pairwise)** — current default pathway (diamond all-vs-all +
  tblastn + per-genome paralog-cutoff filter). Slow; the sensitivity
  reference.
- **Tier C+H (cluster+HMM)** — current `--cluster_tool mmseqs` pathway
  (mmseqs clustering → `hmmbuild` per family → `hmmsearch` per genome).
  Faster than P, still runs one hmmsearch per family per genome.
  Already-run result: `pezizo_set1_cluster`.
- **Tier C (cluster-only, new)** — raw mmseqs family membership taken
  directly as the presence call (a genome is "present" for a family iff at
  least one of its own proteins is a member of that mmseqs cluster). No HMM
  step. For the two existing runs this requires **zero new compute** — it is
  entirely derivable from `families/families_cluster.tsv` and
  `loss_families/families_cluster.tsv`, which already exist on disk.
- **Tier R (refined clustering, new)** — Tier C, plus a targeted pairwise
  refinement pass applied only to **ambiguous families**.

  **Correction (verified against real data while writing the implementation
  plan, replacing an earlier wrong assumption in this section):** ADR-0002's
  family construction clusters **only the ingroup** — mmseqs never sees
  outgroup proteins, so no family's raw membership can "span ingroup and
  outgroup" as originally stated here; every member of every family in
  `families/families_cluster.tsv` carries an ingroup species suffix (verified
  directly: the gain-side file's ~52k member rows resolve to exactly the 5
  ingroup species suffixes, zero outgroup ones). The outgroup "presence" that
  breaks a control (Mcir/Spom for ADA1/HAM5, all 6 outgroups for HEX1) is
  entirely an HMM-search-stage phenomenon (Tier C+H only), or in HEX1's case,
  a consequence of a *within-ingroup* paralog joining the family (see below)
  whose own conservation then drives the HMM hit in every genome — not of the
  paralog itself being an outgroup member.

  The real, available-before-any-HMM signal for HEX1-style contamination is
  **within-ingroup species duplication**: a family whose member count exceeds
  the number of ingroup species means at least one species contributed ≥2
  members — confirmed exactly on HEX1 (6 members over 5 ingroup species:
  Ncra contributes both `HEX1_NEUCR` and `IF5A_NEUCR`) and correctly absent
  from ADA1/HAM5 (5 members, one per ingroup species, no duplication — Tier R
  should leave these two families untouched, matching the earlier
  prediction, now for the right reason). Species-duplication alone is too
  broad on its own to call "ambiguous" (1803 of 7281 gain-side families have
  it — recent real lineage-specific gene duplication is common and mostly not
  a HEX1-style artifact); restricting to families that are otherwise
  novelty-candidate-shaped in Tier C+H's own `presence_matrix.tsv` — ingroup
  presence fraction `>= --ingroup-min-frac` **and** outgroup presence fraction
  `> --other-max-frac` (i.e. currently rejected *only* because of outgroup
  presence) — narrows this to 1227 of 7281 families (~17%) on the gain side,
  measured directly against the published `pezizo_set1_cluster` run. That is
  the definition of "ambiguous family" this investigation uses: **species-
  duplicated AND currently a near-miss novelty candidate**. ~17% of all
  families is not the "small minority" earlier drafts of this section
  assumed, but each family is small (typically 5-30 members), so the total
  refinement cost (below) is still far below Tier P's genome-wide `O(G²)`.

  Within each ambiguous family (bounded further — `oversized_families.tsv`
  already flags pathologically large clusters, excluded from refinement
  rather than paying `O(k²)` on a family with hundreds of members), run a
  small **new** within-family diamond all-vs-all (not a reuse of Tier P's own
  all-vs-all search results, even though those already exist for this study —
  reusing them would not test whether Tier R is viable in a deployment where
  Tier P never ran at all, which is the actual point of this tier), then
  apply the pipeline's **existing, already-validated paralog-competition
  test** (`nf_NovInvenio/bin/parse_self_hits.py` + `build_presence_matrix.py`'s
  filter 2, `--paralog-competition-scope target`): for a family member `X`
  from genome `A`, look up `X`'s own within-genome paralog `P` from that
  genome's already-published `self_hits/<A>.paralog_cutoffs.tsv` (a cheap,
  already-computed, per-genome self-search — no new self-search needed). If
  `P` is *also* a member of the same family (the HEX1 case: `IF5A_NEUCR` is
  `HEX1_NEUCR`'s registered paralog and a co-member), treat the `X`-`P` edge
  as a forced split boundary and remove it from the within-family diamond
  score graph before taking connected components — mirroring "keeps calls
  where the paralog wins on a different gene" (the `target`-scope logic that
  already rescues HEX1 in Tier P) but applied as a graph edit instead of a
  per-hit disqualification, since Tier R's job is to split a family into
  subfamilies, not filter individual presence calls. Recompute presence/
  absence per resulting subfamily using Tier C's own rule (a subfamily is
  present in a genome iff one of its members' `source_proteome` matches).
  A family with no forced-split edge (no member's registered paralog is a
  co-member) passes through unchanged — this is expected to be common even
  within the 1227, since not every within-ingroup duplication pairs two
  *registered* paralogs inside the same family.

  A graph-clustering alternative (MCL on the within-family score graph,
  closer to OrthoFinder's own orthogroup-splitting step, more robust against
  chains of 3+ progressively-diverging paralogs than a single forced-split
  edge) is noted as a fallback if the simpler method proves insufficient
  (e.g. a chain where no single pair is a registered paralog but three or
  more members are each other's near-neighbors), but is not part of the first
  pass — it adds a new dependency and an inflation-parameter tuning knob the
  simpler method doesn't need.

## Analysis 1 — Tiered controls recall/FP (gains)

Extend `nf_NovInvenio/bin/score_controls.py` with a presence-source switch
(e.g. `--presence-mode {hmm,cluster_membership}`; default `hmm`, preserving
current behavior). In `cluster_membership` mode, build each family's presence
vector directly from `--cluster-tsv` membership (a family is present in a
genome iff any of its members' `source_proteome` matches that genome) instead
of reading `presence_matrix.tsv`, then reuse the existing `family_call()` /
`summarize()` unchanged — same novelty predicate, same output contract, so
the three tiers are scored identically.

Run all four tiers against:
- `pezizo_set1.controls.csv` (16 controls, the set analyzed above).
- `Agaricales.controls.csv` (spc14/spc33 fasta-anchor positive controls),
  against the existing `agaricomycetes_pairwise` and `agaricomycetes_mmseqs`
  result pairs — a second, independent ingroup/outgroup split, so findings
  aren't a one-clade artifact.

For Tier R specifically: first build the ambiguous-family detector (species-
duplicated **and** currently a near-miss novelty candidate in Tier C+H's own
`presence_matrix.tsv` — see the corrected definition above), confirm HEX1's
family is in that set (it should be — Ncra contributes 2 members) and that
ADA1/HAM5's families are *not* in it (both are 5-member, one-per-species,
non-duplicated families — Tier R should leave them exactly as Tier C already
scores them — Tier R is not expected to change ADA1/HAM5's outcome, only
HEX1's), then run the within-family forced-split refinement on the ambiguous set
and rescore.

Also fix the `--busco-map` gap: all 5 currently-`unresolved` BUSCO negatives
in *both* runs are unresolved because the map wasn't passed to this
particular scoring invocation, not because of a real method difference — the
negative side of the comparison is currently untested on 5 of 10 controls.
Generate the map (`busco_id → protein_id` from an existing BUSCO run) and
rerun so FP-rate is measured on the full negative set.

Deliverable: a per-control table (control × tier × outcome × failure-mode
classification, generalizing the manual table above) plus one recall/FP
summary row per tier per control-set.

## Analysis 2 — Genome-wide concordance (gains and losses)

16 controls (6 positive) cannot estimate genome-wide precision/recall — it
only tells us about 6 named genes. Use Tier P's full `candidates.txt` (gains)
and `loss_candidates.txt` (losses) as a provisional gold-standard proxy —
explicitly caveated: Tier P scores 1.0 on the curated controls, but that is
not the same claim as "Tier P is ground truth for every candidate it calls."

For Tier C, Tier C+H, and Tier R, compute against that proxy:
- precision / recall / Jaccard overlap of the novelty (and separately, loss)
  candidate sets, mapping between pairwise per-protein candidates and
  cluster per-family candidates via shared cluster membership.
- for every miss, classify it as clustering-level or HMM-level using the same
  anchor-tracing method as Analysis 1 (family membership vs. presence_matrix
  disagreement), so the genome-wide numbers are explained by the same two
  failure modes, not a new unexplained gap.
- **for Tier R specifically, a candidate-count-inflation check**: compare the
  total novelty (and loss) candidate count before vs. after refinement. A
  large increase — more subfamilies calling "novel" than the ambiguous-family
  count itself would predict — is the direct signature of failure mode 4
  (over-splitting), and is the metric this tier is added to watch for, since
  the curated controls cannot surface it (see the failure-mode-4 discussion
  above). Also spot-check whether refinement changed the outcome for any
  *already-resolved-correctly* family (it shouldn't, for non-ambiguous ones by
  construction, but a bug in the ambiguous-family detector would show up
  exactly as an unexpected change here).

**Losses run the identical construction**, against `loss_candidates.txt` /
`loss_presence_matrix.tsv` / `loss_families/families_cluster.tsv`. There are
currently no curated *positive* loss controls (a gene known specifically lost
in one lineage) in this repo, for either pezizo_set1 or the agaricomycetes
set. Per your steer, this investigation proceeds without them: the loss-side
comparison is concordance-only (Tier P loss candidates as the baseline), and
its sensitivity claims are explicitly weaker than the gain-side ones (no true
-positive anchors to confirm recall against) until loss controls are curated
as a follow-on.

## Analysis 3 — Real compute-cost accounting

Pull actual per-process `realtime` from the nextflow trace files already
present in both runs' `nextflow_log/` (diamond/tblastn realtime for Tier P vs.
mmseqs clustering + hmmsearch realtime for Tier C+H; Tier C's cost is the
mmseqs clustering step alone, already paid for by Tier C+H, so its marginal
cost given a Tier C+H run is zero). Aggregate to CPU-hours per tier so the
"how much compute do we save" question has a number, not just "faster."

For Tier R, cost is not free like Tier C but should stay small: measure the
actual wall-clock of the within-family diamond self-search step against the
ambiguous-family count and total member count it covers, and report it as
its own line rather than folding it into Tier C+H's cost — the point of Tier
R is that this additional cost is a small fraction of Tier P's genome-wide
`O(G²)` search, and that claim should be measured, not assumed.

## Analysis 4 — Decision framework

Given each tier's recall/FP/cost profile from Analyses 1-3, define concrete
candidate-prioritization logic, e.g.: run Tier R genome-wide as a cheap first
pass (Tier C's speed plus the small, targeted refinement cost); escalate to
full Tier P confirmation anything Tier R flags, *plus* anything it flags
whose outgroup-presence fraction sits near the `other_max_frac` boundary
(the exact pattern that broke ADA1/HAM5) rather than comfortably below it.
Whether Tier R actually earns a place ahead of plain Tier C in that pipeline
depends on Analyses 1-2: it only replaces Tier C if it recovers HEX1-style
misses (Analysis 1) without a material candidate-count inflation (Analysis
2's over-splitting check) — if Tier R's over-splitting risk turns out to be
worse than its recall gain in practice, the recommendation should say so
plainly rather than defaulting to "more refinement is always better."

State explicitly what would invalidate this framework: if genome-wide
concordance (Analysis 2) shows HEX1-style clustering-inherent paralog
inflation is common rather than a one-off *and* Tier R's targeted refinement
doesn't scale to that volume (e.g. too many families turn out ambiguous, or
the simple pairwise cutoff fails on chains of 3+ paralogs and would need the
heavier MCL fallback), no cheap tier recovers those misses — the framework
would then need a clustering-quality fix upstream of all of this (e.g.
building the paralog-competition filter into mmseqs clustering itself), not
just a downstream refinement or presence rule.

## Scope and where code lives

- `score_controls.py`'s new `--presence-mode` flag is pipeline code that any
  future study's controls-scoring can reuse → `nf_NovInvenio/bin/score_controls.py`
  (extend in place, default preserves current `hmm` behavior).
- The Analysis 2/3 comparison runner (genome-wide concordance + trace-based
  cost accounting) is reusable across any future pairwise/cluster study pair,
  not specific to pezizo_set1 → `NovInvenio_Investigations/bin/`, per this
  repo's placement rule for enrichment/comparison scripts that read a study's
  presence matrices.
- Tier R's ambiguous-family detector and within-family paralog-cutoff split
  are pipeline-adjacent logic (they read `families_cluster.tsv` and reuse the
  existing diamond self-hit + e-value-cutoff paralog-separation method), so
  they belong next to `score_controls.py` → `nf_NovInvenio/bin/`, as a new
  script (e.g. `refine_ambiguous_families.py`) rather than folded into
  `score_controls.py` itself, since its output (a refined `families_cluster`-
  style membership table) is a general input any consumer of family
  membership could use, not just controls scoring.
- Output: a single markdown/TSV report (per-control table, per-tier summary,
  genome-wide concordance numbers, cost table, and the decision-framework
  recommendation) written to this study's directory
  (`studies/fungi/pezizo_set1_cluster/`), since it compares this study's two
  runs specifically; the underlying comparison script itself stays reusable.

### Promotion path if Tier R validates

Everything above is deliberately built as standalone scripts run against
**already-published** results (`families_cluster.tsv`, `presence_matrix.tsv`,
trace files already on disk) — nothing here is wired into a Nextflow process,
and this investigation does not touch `nf_NovInvenio`'s DSL2 workflows. That
is intentional: it lets Tier R get evaluated (Analyses 1-2) before any
pipeline surface area is committed to it.

If Analyses 1-2 show Tier R recovers HEX1-style misses without material
candidate-count inflation, `refine_ambiguous_families.py`'s logic is the
starting point for a **new Nextflow process** inside the
`--cluster_tool mmseqs` pathway (`workflows/cluster.nf` /
`workflows/profile_search.nf`), following the same per-chunk pattern already
used by `BUILD_FAMILY_PROFILES:BUILD_CHUNK`, so refinement runs as a normal
scheduled step rather than a manual script a user has to remember to invoke.
It should land behind a new flag (e.g. `--refine_ambiguous_families`,
defaulting off) so existing `--cluster_tool mmseqs` runs are unaffected until
this is validated on more than one study. That Nextflow-integration work is
its own follow-on plan, scoped after this investigation's results are in
hand — not part of this spec's deliverable.

## Explicitly out of scope

- **A general, production-grade fix to mmseqs clustering itself** (e.g.
  building a paralog-competition filter into the clustering step for every
  family, or a full OrthoFinder-style all-orthogroups MCL+gene-tree
  pipeline). Tier R is a scoped experiment — targeted refinement of only the
  small ambiguous-family subset, using the pipeline's existing validated
  paralog-cutoff method — not a redesign of clustering. If Tier R's results
  show the simple cutoff is insufficient (chains of 3+ paralogs, or the
  ambiguous set turns out too large to stay cheap), that redesign becomes its
  own future spec, informed by what this investigation measures.
- The 2026-09-12 resume run's loss-side task failures (101/123 chunks) — a
  separate pipeline-reliability bug, unrelated to the published results this
  investigation is built on.
- Curating new loss-direction positive controls — noted as a follow-on, not
  part of this investigation's deliverable.
- **Porting any validated tier into an `nf_NovInvenio` Nextflow process** —
  see "Promotion path if Tier R validates" above. This spec's deliverable is
  the standalone-script comparison and its recommendation; turning a
  validated result into a first-class pipeline step (with its own flag,
  process definition, and resource profile) is separate follow-on work.
