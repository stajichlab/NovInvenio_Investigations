# Tier-2 clustering wiring + bin/ generalization recommendations (2026-09-16)

Recommendations only — not implemented here. The user has also asked the
coccidioides-testing session ("genome-data-migration-spec") to look at these same
two open items on the `NovInvenio` `pangenome-profiling-module` branch, so this file
exists to hand off findings, not to claim the work.

## 1. Tier-2 clustering wiring

**Current state** (`NovInvenio/bin/pangenome_cluster_backend.py`,
`NovInvenio/modules/pangenome/prefix_and_cluster.nf`): tier-1 clustering is fully
wired into the subworkflow (`CLUSTER_TIER1` → `PRESENCE_MATRIX`,
`FAMILY_POSITIONS`, `SELECT_BACKGROUND_REPS`, etc., see
`workflows/pangenome_profile.nf`). Tier-2 exists only as two CLI subcommands
(`mmseqs-tier2`/`diamond-tier2`) plus a library function `two_tier_families()`
(`pangenome_cluster_backend.py:85-96`) that **no script or CLI subcommand ever
calls** — dead code from the subworkflow's point of view, reachable only by
importing the module directly. Never run in the original study either (no
`results/*tier2*` output anywhere) and not being worked on by the coccidioides
track as of this check (zero tier2/superfamily references there).

**What's needed to wire it in, three concrete steps:**

1. **New `CLUSTER_TIER2` process** in `modules/pangenome/prefix_and_cluster.nf`,
   modeled directly on `CLUSTER_TIER1` (same mmseqs/diamond backend dispatch
   pattern, lines 68-119): input `CLUSTER_TIER1.out.rep_fasta`, calls
   `pangenome_cluster_backend.py mmseqs-tier2`/`diamond-tier2`, outputs
   `tier2_cluster.tsv`. Nearly a copy-paste of the existing process with
   different defaults (0.4 identity vs 0.9) — low risk.

2. **A CLI entry point for `two_tier_families()` — the actual gap, not just
   missing wiring.** The combiner function exists but has never been exposed as
   a runnable subcommand or written to a file anywhere, in either repo. Add a
   `combine-tiers` subcommand to `pangenome_cluster_backend.py` (or a small new
   script) that calls `two_tier_families(tier1_tsv, tier2_tsv)` and writes a
   `family_superfamily.tsv` (columns: tier1_family, superfamily, member_count).
   Needs its own unit test before trusting it on real data — nothing currently
   exercises this function end-to-end.

3. **Wire the output in as annotation only, not into the presence/frequency
   path.** Per the docstring, tier-2 is deliberately a superfamily *label*,
   never used for presence/frequency calls — so the safe integration point is
   alongside `PRESENCE_MATRIX`/`workflow.emit` (line 258-259 pattern), adding
   `superfamily_tsv = CLUSTER_TIER2.out.cluster_tsv` as an extra emitted
   output, plus a new `params.pangenome_run_tier2` boolean gate (default
   `false`, matching how other optional stages are gated in this same file) so
   it costs nothing for studies that don't want it.

**Sequencing recommendation**: do step 2 (the CLI wrapper + its test) before
step 1's Nextflow wiring — there's no way to validate the process's output shape
without a script that actually produces the combined TSV first. Small,
self-contained task: one new function, one new subcommand, one process, one
param — not a redesign.

## 2. `bin/` generalization / namespace question

Checked what's already been ported: **19 of this study's ~34 `bin/` scripts
already exist in `NovInvenio/bin/`** under a flat `pangenome_<name>.py` prefix
convention (not a subdirectory namespace — matches that repo's existing flat
`bin/` layout: `pangenome_cluster_backend.py`, `pangenome_domain_enrichment.py`,
`pangenome_dereplicate_strains.py`, etc.). The generalization is already
underway and the naming convention is already established — no new naming
decision needed, just continuation of the existing pattern.

Already ported (rough correspondence, study name → NovInvenio name):
`assign_clades.py`, `build_family_positions.py`, `build_gene_positions.py`,
`build_presence_matrix.py`, `cluster_backend.py`, `cooccurrence.py`,
`dereplicate_strains.py`, `extract_absent_family_queries.py`,
`extract_rescue_positions.py`, `frequency_bins.py`, `pair_classification.py`,
`rescue_pass.py`, `find_accessory_islands.py`+`synteny_windows.py` (combined into
`pangenome_build_islands.py`), `summarize_island_functions.py` (verbatim port as
`pangenome_domain_enrichment.py`, confirmed by its own docstring).

**Not yet ported** (study-only today, split by whether the logic is
dataset-generic or *A. fumigatus*-specific):

- **Generic, worth porting**: `detect_trans_modules.py` (Leiden module
  detection — currently the biggest capability gap, since the Nextflow port has
  no module-detection step at all), the SwissProt/GO/Pfam-URL annotation chain
  (`annotate_families_with_swissprot.py`, `add_swissprot_to_islands.py`,
  `map_pfam_to_go.py`, `add_pfam_urls.py` — pure joins against standard
  databases, no *A. fumigatus* assumptions), `add_island_locus.py`, and
  `cluster_homology_test.py` (already written generic — takes any
  query/reference FASTA pair).
- **Study-specific, keep local, do not port**: `benchmark_scorecard.py`,
  `fetch_paper_reference_proteins.py`, `parse_starship_supplement.py`,
  `hac_screen.py`/`hac_reference_screen.py`, `build_ground_truth_tables.py` —
  these encode this one paper's ground truth (Starship supplement tables,
  HAC/hacA screen, the AF293 benchmark control), not reusable pangenome
  methodology.
- Also not ported, not yet triaged: `dereplication_stability.py`,
  `leiden_stability.py` (built 2026-09-16 for parameter-validation, arguably
  generic methodology but tied to this study's specific ground-truth files —
  worth a second look before porting as-is).

## Status

Recommendations only, per explicit request — no code changes made against
either repo for either item in this pass.
