# Generalized accessory-island + Pfam functional-enrichment pipeline step

Date: 2026-09-16
Repo: `nf_NovInvenio` (local checkout: `/bigdata/stajichlab/jstajich/projects/NovInvenio`,
branch `pangenome-profiling-module`)
Reviewed with: two rounds of Fable review (module granularity/correctness), both
findings folded in below.

## Problem

The Afumigatus pangenome study (`NovInvenio_Investigations/studies/fungi/
Afumigatus_pangenome/`) manually built a real, working accessory-island detection +
Pfam functional-enrichment analysis as ad hoc `bin/` scripts and raw SLURM shell
scripts — never wired into the Nextflow pipeline. This spec generalizes that proven
methodology into a permanent, flag-gated step in `nf_NovInvenio`'s
`pangenome_profile.nf` subworkflow, so any future pangenome study (starting with
Coccidioides) gets it as a pipeline capability, not a copy-pasted script chain.

**Explicitly not new statistical work**: the co-occurrence permutation test
(`pangenome_cooccurrence.py`'s `exact_stratified_pvalue()`) and the pairwise FDR
(`fdr_q < fdr_alpha` gate in `pangenome_cooccurrence.py`) are already implemented,
exact (not Monte Carlo), and unchanged by this spec. This step *consumes* their
output (`pair_classification.tsv`'s `classification` column) — it does not
recompute or replace them. The two new statistical elements here are (a) island
construction (a deterministic algorithm, not a hypothesis test) and (b) a *separate*
Fisher-exact-per-Pfam-domain enrichment test with its own BH correction — a
different question (is this domain over-represented in flagged islands vs.
background?) from the co-occurrence test.

## Fact base (verified against real code/data, 2026-09-16)

### Island construction — `synteny_windows.py`'s `accessory_islands()`

Real, tested algorithm (`studies/fungi/Afumigatus_pangenome/bin/synteny_windows.py`,
lines 62-113): given one strain's gene order (a list of `(gene_id, contig, ...)`
tuples in contig-then-start-sorted order) and an `is_core: dict[gene_id, bool]` map,
walks the list and merges **maximal consecutive runs of non-core genes on the same
contig** into one island. Any core gene, or a contig boundary, ends the current run.
A gene absent from `is_core` defaults to **core** (conservative — an unrecognized
gene more likely reflects an ID mismatch than a genuine accessory gene), and the
function **warns if >10% of genes are missing from `is_core`** — that's a real
safety mechanism (catches a GFF3-ID-vs-presence-matrix-ID format mismatch that would
otherwise silently collapse every island to nothing) and must be ported as-is, not
dropped.

`find_accessory_islands.py` (the real-data driver, 236 lines) reuses
`family_positions.rescued.tsv`'s already-computed per-strain `rank` column directly
as the gene-order input (no GFF3 re-parsing — `accessory_islands()` only needs
`(gene_id, contig)` pairs in traversal order, and family_positions.tsv's rank column
already is that order). `is_core` comes from `frequency_table.tsv`'s `bin` column:
`bin in ("core", "soft_core")`. Islands are then filtered to keep only those
containing **at least one pair whose `pair_classification.tsv` `classification` is
in `{starship_explained, unexplained_physical, ambiguous_linkage}`** (the
`PHYSICAL_CLASSIFICATIONS` frozenset, `find_accessory_islands.py` line ~47) — this
*is* the significance gate; verified directly, it never touches raw
`permutation_p`/`fdr_q`. That's a legitimate simplification, not a corner cut: those
three labels come from `classify_pair()`'s `linkage_fraction` threshold (`>=0.5`
physical, `0.05-0.5` ambiguous — `pangenome_pair_classification.py` lines 149-180,
`DEFAULT_PHYSICAL_THRESHOLD=0.5`), and every row in `pair_classification.tsv`
already passed the upstream `fdr_q < fdr_alpha` filter in `pangenome_cooccurrence.py`
(line ~294) before classification ever runs — so "significant, physically-linked
pair" is already fully encoded in `classification`, nothing to re-derive.

Real Afumigatus output (`results/accessory_islands/significant_islands.tsv`,
confirmed on disk, 8.4 MB): 12,861 distinct significant islands, size 2-837
(median 7). Columns: `n_strains, example_strain, island_size, member_families,
n_supporting_pairs, classifications, has_captain_gene, has_sm_backbone_gene`.

### Two genuinely separate marker searches — NOT one column

Verified directly (`find_accessory_islands.py` CLI: `--captain_tblout`,
`--sm_backbone_tblout`, both optional, independent files): `has_captain_gene` comes
from the existing `CAPTAIN_HMMSEARCH` pipeline module (DUF3435, Starship-specific).
`has_sm_backbone_gene` comes from a **second, unrelated `hmmsearch`** —
`results/sm_backbone/SM_backbone.hmm` (confirmed on disk) is PF00109
(ketoacyl-synt) + PF00668 (Condensation) concatenated, run against the full
all-strains proteome, entirely independent of the captain search. The 2×2
breakdown (11,741 neither / 904 captain-only / 151 SM-only / 65 both) was a
headline Afumigatus result. **A single `has_marker_gene` column would silently
collapse two real, independent signals** — this was the design's original mistake,
caught by Fable review, and is fixed below (§ Marker-gene handling).

### Pfam domain enrichment — `summarize_island_functions.py`

Real algorithm: background = **all shell+cloud-bin families** (`load_eligible_families()`,
line ~70 — matches `pangenome_cooccurrence.py`'s own co-occurrence-eligibility
selection, `row["bin"] in ("shell", "cloud")` — NOT core/soft_core/singleton).
Singleton island members are explicitly excluded from enrichment testing (a
singleton can never appear in the shell+cloud background set, since background must
be a superset of anything being tested against it — `domain_enrichment()` line ~100
comment). For each Pfam domain seen in the domtblout: one-sided Fisher exact test
(island-member families with that domain vs. background families with that
domain), then `scipy.stats.false_discovery_control(..., method="bh")` across every
domain tested (901 domains in the real run, 71 significant at q<0.05).

`run_island_pfam_scan.sh`/`run_background_pfam_scan.sh` (raw SLURM scripts, not
pipeline code) ran `hmmscan --domtblout ... -E 1e-3 --cpu N Pfam-A.hmm <reps.fa>`
**twice** — once for 9,224 island-member family reps, once for 2,620 remaining
shell+cloud "background" reps (together, all 11,052 eligible families). This
two-scan split has no principled reason (background must already be a superset of
island members for the Fisher test to be valid) and is corrected in this design
(§ Architecture) to one scan over the full eligible set.

### Report format

`results/REPORT.md` (21 KB, real file) — Markdown with `![...](figures/*.png)`
embeds, parallel `figures_pdf/*.pdf` vector copies, hand-assembled narrative
(correction history, synthesis, open items — not purely templated). Figures via
**matplotlib** (`Agg` backend, `bin/plot_pangenome_summary.py`): frequency
histogram, core/shell/cloud pie, presence/absence raster (`imshow`, not per-cell —
kept fast at tens-of-thousands-of-families scale), rarefaction curve (20
random-order permutations), pair-classification bar chart. `bin/build_summary_report.py`
is pure aggregation (reads existing TSVs → Markdown tables, no recomputation) →
`SUMMARY.md` (a smaller, fully machine-generated companion, distinct from the
hand-assembled `REPORT.md`). `lib/report_template.py`/`core_report_template.py`/
`losses_report_template.py` (in the pipeline repo) are **not** used for any of this
— they generate the interactive HTML novelty/loss reports, an unrelated report
family. This design's `PANGENOME_REPORT` output must be **fully templated/generated**
(no hand-written narrative), since it has to work for any future study
automatically, not just Afumigatus.

## Architecture

All new scripts (`bin/pangenome_build_islands.py`,
`bin/pangenome_select_background_reps.py`, `bin/pangenome_domain_enrichment.py`,
`bin/pangenome_report_tables.py`, `bin/pangenome_report_render.py`) live in
`nf_NovInvenio`'s own shared `bin/` — this is pipeline code, not study-specific
(unlike the Coccidioides onboarding scripts, which live in that study's own
`bin/` per this repo's study-specific-vs-shared convention). Same for the new
Nextflow modules and workflow wiring.

Six new modules in `modules/pangenome/`, wired into `workflows/pangenome_profile.nf`
after the existing `GENE_POSITIONS`/`FAMILY_POSITIONS`/`PAIR_CLASSIFICATION` steps,
all flag-gated behind a new required param `--pangenome_pfam_hmm` (reusing that
param name from the existing novelty/loss workflow for consistency — off by
default, matching the existing `--pangenome_captain_hmm` optional-branch pattern).

```
FAMILY_POSITIONS.out ─┐
PAIR_CLASSIFICATION.out ─┤
FREQUENCY_BINS.out ─────┼──▶ BUILD_ISLANDS ──▶ significant_islands.tsv ─┐
CLUSTER_TIER1.out (cluster_tsv) ┘                                       │
MARKER_HMMSEARCH (0+, named) ──────────────────────────────────────────┘
                                                                          │
CLUSTER_TIER1.out (tier1_rep_seq.fasta) ─┐                               │
FREQUENCY_BINS.out ───────────────────────┼─▶ SELECT_BACKGROUND_REPS      │
                                          │    ──▶ background_reps.fa    │
                                          │         │                    │
                          --pangenome_pfam_hmm ──▶ FAMILY_PFAM_SCAN      │
                                                     ──▶ pfam.domtblout   │
                                                              │          │
                                                              ▼          ▼
                                                        DOMAIN_ENRICHMENT
                                                   (partitions domtblout by
                                                    island membership itself)
                                                              │
                                                              ▼
                                              island_pfam_enrichment.tsv
                                                              │
              significant_islands.tsv ──────────────────────┤
              pair_classification.tsv, frequency_table.tsv,  │
              presence_matrix.tsv ────────────────────────────┤
                                                              ▼
                                                       REPORT_TABLES
                                                   ──▶ tidy TSVs (islands
                                                       joined to enriched
                                                       domains, size dist,
                                                       classification counts)
                                                              │
                                                              ▼
                                                       REPORT_RENDER
                                            ──▶ report.md, figures/*.png,
                                                figures_pdf/*.pdf
```

Key correction from the first draft (both from Fable review, verified above):
**one** Pfam scan over the full shell+cloud background set, not two — island vs.
background partitioning happens in `DOMAIN_ENRICHMENT`, which decouples the scan
from `BUILD_ISLANDS` entirely (they run in parallel) and matches
`summarize_island_functions.py`'s actual logic exactly.

### Module specs

**`BUILD_ISLANDS`**
- Script: `bin/pangenome_build_islands.py` (new, ports `synteny_windows.py`'s
  `accessory_islands()` verbatim including the 10%-missing warning, and
  `find_accessory_islands.py`'s driver logic, generalized per marker handling below)
- Inputs: `family_positions.tsv`, `frequency_table.tsv`, `pair_classification.tsv`,
  `cluster_tsv` (needed to map marker-hit proteins back to families — present in
  the real script, missing from the first draft), 0+ named marker hit tables
  (`--marker_tblout name=path`, repeatable)
- Params consumed: `--pangenome_island_min_size` (int, default 2)
- Output: `significant_islands.tsv` — columns `n_strains, example_strain,
  island_size, member_families, n_supporting_pairs, classifications,
  has_<marker_name>` (one column per named marker passed in, zero if none)

**`MARKER_HMMSEARCH`** (0+ instances, one per named marker set)
- Reuses the *existing* `HMMFETCH_CAPTAIN`/`CAPTAIN_HMMSEARCH` modules
  (`modules/pangenome/captain.nf`) as-is — already generic (fetch-by-Pfam-name or
  direct HMM path, `hmmsearch` vs. the full concatenated proteome). No new module
  code; this is a wiring/aliasing change (one call per entry in a new params list),
  not new pipeline logic.
- Params: `--pangenome_marker_names` (comma list, e.g. `captain,sm_backbone`),
  `--pangenome_marker_hmm_paths` (parallel comma list of HMM file paths, matched by
  position — same comma-list convention already used elsewhere in this pipeline,
  e.g. `--groups 'IN,OUT'`). Empty by default (no marker searches run). The
  workflow must hard-error at samplesheet-parsing time (matching this pipeline's
  existing "validate now, not deep inside a later process failure" convention —
  see `pangenome.nf`'s GFF3-existence check) if the two lists have different
  lengths, rather than silently zipping to the shorter one.
- This directly replaces the single `has_marker_gene` column mistake: each named
  marker gets its own `has_<name>` column in `significant_islands.tsv`, preserving
  the captain-vs-SM-backbone (or any future named marker) distinction.

**`SELECT_BACKGROUND_REPS`**
- Script: `bin/pangenome_select_background_reps.py` (new, trivial — filters
  `tier1_rep_seq.fasta` to headers whose family is `shell` or `cloud` per
  `frequency_table.tsv`)
- Output: `background_reps.fa` (the full eligible set — island members are a subset
  of this, not a separate file)

**`FAMILY_PFAM_SCAN`**
- One `hmmscan --domtblout ... -E {--pangenome_pfam_domain_evalue, default 1e-3}`
  call against `background_reps.fa`, using `--pangenome_pfam_hmm`
- Output: `pfam.domtblout` (one file, not two)

**`DOMAIN_ENRICHMENT`**
- Script: `bin/pangenome_domain_enrichment.py` (new, ports
  `summarize_island_functions.py`'s `domain_enrichment()` verbatim: one-sided
  Fisher exact per domain, BH correction via `scipy.stats.false_discovery_control`
  across all domains tested)
- Inputs: `pfam.domtblout`, `significant_islands.tsv` (to know which families are
  island members), `frequency_table.tsv` (background set definition)
- Output: `island_pfam_enrichment.tsv` — columns `pfam_domain, n_island,
  n_background, odds_ratio, fisher_p, fdr_q`

**`REPORT_TABLES`**
- Script: `bin/pangenome_report_tables.py` (new, ports
  `annotate_islands_with_enrichment.py`'s join logic)
- Inputs: `significant_islands.tsv`, `island_pfam_enrichment.tsv`,
  `frequency_table.tsv`, `pair_classification.tsv`
- Output: tidy, testable-without-matplotlib TSVs (islands joined to their enriched
  domains, island-size distribution, classification counts) — also independently
  useful to this repo's `docs/` publishing pipeline, not just the report step

**`REPORT_RENDER`**
- Script: `bin/pangenome_report_render.py` (new, ports `plot_pangenome_summary.py`'s
  matplotlib figure generation — frequency histogram, core/shell/cloud pie,
  presence/absence raster, rarefaction curve, classification bar chart, **plus**
  new figures for island-size distribution and top enriched domains — and a
  templated Markdown generator, fully automated, no hand-written narrative)
- Inputs: `REPORT_TABLES`'s output, `presence_matrix.tsv`
- Output: `report.md`, `figures/*.png`, `figures_pdf/*.pdf` (both formats per figure,
  matching Afumigatus's convention)

## Params summary (new)

| Param | Type | Default | Purpose |
|---|---|---|---|
| `--pangenome_pfam_hmm` | path | none (branch off) | Enables the whole island+Pfam step |
| `--pangenome_island_min_size` | int | 2 | Minimum island size to report |
| `--pangenome_pfam_domain_evalue` | float | 1e-3 | hmmscan domain-level E-value cutoff |
| `--pangenome_marker_names` | comma list | empty | Named marker searches (e.g. `captain,sm_backbone`) |
| `--pangenome_marker_hmm_paths` | comma list | empty | Parallel HMM paths for each named marker |

## Explicitly out of scope for this spec

- Publishing a Claude Artifact — that's an agent-side action performed per-study
  after the pipeline runs, not a pipeline capability (a Nextflow module cannot
  invoke the Artifact tool). `REPORT_RENDER`'s Markdown+PNG+PDF output is what a
  study's own agent session would read to build an Artifact from, separately.
- The MULE/DDE-transposase marker idea (deferred — see
  `nf_NovInvenio/todo/mule-dde-transposase-island-marker.md`) — this spec's
  `--pangenome_marker_names` mechanism is exactly what that idea would plug into
  later, once specific Pfam accessions are identified, but no accessions are
  chosen here.
- A new island-level significance test (e.g. "is this island's size/clustering
  more significant than a null model") — not requested; the existing two-pass
  approach (pairwise FDR gates island flagging, separate domain-enrichment FDR)
  is what's being generalized, not replaced or extended with a third test.
- Applying this to Coccidioides' actual data (a separate execution step, once this
  spec is approved and turned into an implementation plan).
