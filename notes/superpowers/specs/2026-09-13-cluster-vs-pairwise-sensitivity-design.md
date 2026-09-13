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
both directions (novel gains and lineage-specific losses)?

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

## Three tiers under comparison

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

## Analysis 1 — Tiered controls recall/FP (gains)

Extend `nf_NovInvenio/bin/score_controls.py` with a presence-source switch
(e.g. `--presence-mode {hmm,cluster_membership}`; default `hmm`, preserving
current behavior). In `cluster_membership` mode, build each family's presence
vector directly from `--cluster-tsv` membership (a family is present in a
genome iff any of its members' `source_proteome` matches that genome) instead
of reading `presence_matrix.tsv`, then reuse the existing `family_call()` /
`summarize()` unchanged — same novelty predicate, same output contract, so
the three tiers are scored identically.

Run all three tiers against:
- `pezizo_set1.controls.csv` (16 controls, the set analyzed above).
- `Agaricales.controls.csv` (spc14/spc33 fasta-anchor positive controls),
  against the existing `agaricomycetes_pairwise` and `agaricomycetes_mmseqs`
  result pairs — a second, independent ingroup/outgroup split, so findings
  aren't a one-clade artifact.

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

For Tier C and Tier C+H, compute against that proxy:
- precision / recall / Jaccard overlap of the novelty (and separately, loss)
  candidate sets, mapping between pairwise per-protein candidates and
  cluster per-family candidates via shared cluster membership.
- for every miss, classify it as clustering-level or HMM-level using the same
  anchor-tracing method as Analysis 1 (family membership vs. presence_matrix
  disagreement), so the genome-wide numbers are explained by the same two
  failure modes, not a new unexplained gap.

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

## Analysis 4 — Decision framework

Given each tier's recall/FP/cost profile from Analyses 1-3, define concrete
candidate-prioritization logic, e.g.: run Tier C genome-wide as a cheap first
pass; escalate to full Tier P confirmation anything Tier C flags, *plus*
anything Tier C+H flags whose outgroup-presence fraction sits near the
`other_max_frac` boundary (the exact pattern that broke ADA1/HAM5) rather than
comfortably below it.

State explicitly what would invalidate this framework: if genome-wide
concordance (Analysis 2) shows HEX1-style clustering-inherent paralog
inflation is common rather than a one-off, no cheap tier recovers those
misses — the framework would need a clustering-quality fix (e.g. a
cluster-time paralog-competition filter analogous to the pairwise pipeline's
per-genome self-hit cutoff), not just a different presence rule downstream of
clustering.

## Scope and where code lives

- `score_controls.py`'s new `--presence-mode` flag is pipeline code that any
  future study's controls-scoring can reuse → `nf_NovInvenio/bin/score_controls.py`
  (extend in place, default preserves current `hmm` behavior).
- The Analysis 2/3 comparison runner (genome-wide concordance + trace-based
  cost accounting) is reusable across any future pairwise/cluster study pair,
  not specific to pezizo_set1 → `NovInvenio_Investigations/bin/`, per this
  repo's placement rule for enrichment/comparison scripts that read a study's
  presence matrices.
- Output: a single markdown/TSV report (per-control table, per-tier summary,
  genome-wide concordance numbers, cost table, and the decision-framework
  recommendation) written to this study's directory
  (`studies/fungi/pezizo_set1_cluster/`), since it compares this study's two
  runs specifically; the underlying comparison script itself stays reusable.

## Explicitly out of scope

- Fixing the mmseqs clustering paralog-inflation problem itself (HEX1-style)
  — this investigation measures and classifies that failure mode, it does
  not design its fix.
- The 2026-09-12 resume run's loss-side task failures (101/123 chunks) — a
  separate pipeline-reliability bug, unrelated to the published results this
  investigation is built on.
- Curating new loss-direction positive controls — noted as a follow-on, not
  part of this investigation's deliverable.
