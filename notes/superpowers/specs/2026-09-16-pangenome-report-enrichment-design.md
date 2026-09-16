# Pangenome Report Enrichment Bundle — Design

**Sub-project A** of a 4-part porting effort from `studies/fungi/Afumigatus_pangenome/`'s
newest analysis scripts into `nf_NovInvenio`'s generalized pangenome pipeline step
(island detection + Pfam enrichment, landed 2026-09-16 on `pangenome-profiling-module`,
now merged to `main`). Sub-projects B (SwissProt family annotation), C
(`cluster_homology_test.py` standalone port), and D (dereplication/Leiden stability
tools) are separate, later specs.

Repo: `/bigdata/stajichlab/jstajich/projects/NovInvenio` (nf_NovInvenio), branch TBD at
plan time (not main/master without explicit consent).

## Motivation

The Afumigatus_pangenome study's newest work (dereplication/Leiden stability,
`Asfu_H1106` outlier resolution, SwissProt/Pfam2GO cross-validation, InterPro
hotlinking) produced four capabilities that are genuinely generalizable — no
Afumigatus-specific hardcoding — and cheap to port, with no new runtime
infrastructure beyond one operator-supplied reference file:

1. InterPro hotlinks for Pfam domains in the rendered report
2. Pfam2GO GO-term annotation of enriched domains
3. Island genomic-locus assignment (citable `contig:start-end` coordinates)
4. Per-strain outlier flag (z-score on family-bin counts, catching assembly-
   fragmentation artifacts like the real `Asfu_H1106` case)

None of these were in scope for the original island/Pfam-enrichment plan; they were
identified during a post-hoc review of what else the Afumigatus study had since done
manually. See `research/` notes (agent reports) from 2026-09-16 for the full survey
of PORTABLE vs. STUDY-SPECIFIC Afumigatus code.

## Global constraints (carried from the prior island/Pfam-enrichment plan)

- A Nextflow process can only be invoked ONCE per workflow — multiplicity via channel
  fan-out, never a Groovy loop.
- Every new gating param is distinct from any pre-existing param of a similar name —
  no silent reuse/collision (this bit the prior plan once: `pangenome_pfam_hmm` vs.
  `pangenome_island_pfam_hmm`).
- No hand-written narrative/prose in the generated report — fully templated from
  tabular data, matching `render_report_markdown()`'s existing convention.
- Real 10-strain execution smoke test required before any task touching
  `workflows/pangenome_profile.nf` is considered done — `--help` alone is not
  sufficient (this is what caught F2/F3 in the prior plan; both were report-content
  bugs invisible to a `--help`-only check).
- Prefer `bin/pangenome_*.py` naming (repo convention, 100% consistent today per the
  2026-09-16 bin/lib reorg survey — do not deviate).
- `lib/compressed_io.py` stays unprefixed and untouched by this work (explicit user
  decision 2026-09-16: it's intended to become genuinely cross-pipeline shared, not
  pangenome-only, revisit separately).

## Component 1 — InterPro hotlinks (always on, zero new params)

**Files:**
- Modify: `bin/pangenome_domain_enrichment.py` — add a `pfam_url` column to
  `island_pfam_enrichment.tsv`'s existing output, derived from the already-parsed
  `pfam_accession` column (`parse_domtblout_accessions()`, existing function).
  URL form: `https://www.ebi.ac.uk/interpro/entry/pfam/<PFxxxxx>/` — strip any
  version suffix (`.9` in `PF13577.9`) before building the URL. When
  `pfam_accession` is `-` (accession lookup failed), `pfam_url` is also `-`.
- Modify: `bin/pangenome_report_render.py` — `render_report_markdown()`'s "Pfam
  domain enrichment" table (currently plain `| {domain} | {fisher_p} | {fdr_q} |`
  rows, `render_report_markdown()` around line 331-338) renders
  `[{domain}]({pfam_url})` instead of a bare domain name when `pfam_url` is present
  and not `-`; falls back to the bare name otherwise. No change to
  `plot_domain_enrichment()` (figures don't carry hyperlinks).

**No new Nextflow process, no new param.** This is an in-place augmentation of an
existing script's output plus a render-time formatting change.

**Testing:** unit test for the URL-derivation helper (version-suffix stripping,
`-` passthrough), unit test for the render function choosing link vs. bare-name
based on column presence/value.

## Component 2 — Pfam2GO annotation (gated on `--pangenome_pfam2go`)

**Files:**
- Create: `bin/pangenome_pfam2go.py` — ported from
  `studies/fungi/Afumigatus_pangenome/bin/map_pfam_to_go.py`. Parses a standard
  `pfam2go` file (format: `Pfam:PFxxxxx <name> > GO:<term> ; GO:nnnnnnn` lines) into
  `{pfam_accession: [go_id, ...]}`, reads `island_pfam_enrichment.tsv` (Component 1's
  output, so it has `pfam_accession` already), and writes the same table with two
  added columns: `n_go_terms`, `go_terms` (semicolon-joined GO IDs, or `-` if none
  mapped). CLI: `--island_pfam_enrichment --pfam2go --output`.
- Create: `modules/pangenome/pfam2go.nf` — new `PFAM2GO` process, `input: path
  island_pfam_enrichment, path pfam2go`, `output: path("*.go.tsv"), emit: annotated`.
  Only invoked when `params.pangenome_pfam2go` is truthy (mirrors the existing
  `params.pangenome_island_pfam_hmm`-gated block's conditional-inclusion pattern in
  `workflows/pangenome_profile.nf`).
- Modify: `workflows/pangenome_profile.nf` — when `params.pangenome_pfam2go` is set,
  route `DOMAIN_ENRICHMENT.out.enrichment` through `PFAM2GO` before it reaches
  `REPORT_TABLES`; otherwise pass `DOMAIN_ENRICHMENT.out.enrichment` straight through
  unchanged. Same conditional-channel pattern already used for the marker-search
  branch.
- Modify: `nextflow.config` — add `pangenome_pfam2go = null` (default off).
- Modify: `pangenome.nf` help text + `README.md` — document `--pangenome_pfam2go`
  with explicit operator instructions: download
  `http://current.geneontology.org/ontology/external2go/pfam2go` and pass its local
  path. State plainly that this is a manually-obtained reference file, not
  auto-fetched, and that the step is skipped entirely (no error) when the param is
  unset.
- Modify: `bin/pangenome_report_render.py` — when `go_terms`/`n_go_terms` columns
  are present in the domain-enrichment input, add them to the existing "Pfam domain
  enrichment" markdown table (extra columns, not a new section) — omit entirely when
  absent (matches Component 1's presence-conditional pattern).

**Testing:** unit tests for `bin/pangenome_pfam2go.py`'s parser (real pfam2go-format
fixture, multiple GO terms per accession, an accession with none) and its
join/output logic; workflow-level test via the real-execution smoke test with and
without `--pangenome_pfam2go` set.

## Component 3 — Island genomic locus (always on, no new param)

**Files:**
- Create: `bin/pangenome_island_locus.py` — ported from
  `studies/fungi/Afumigatus_pangenome/bin/add_island_locus.py`. Reads
  `significant_islands.tsv` (columns: `n_strains, example_strain, island_size,
  member_families, n_supporting_pairs, classifications[, has_<marker>...]` — exact
  header confirmed from `bin/pangenome_build_islands.py`'s `main()`), inverts
  `cluster_tsv` (tier1_cluster.tsv, rep→members via the existing
  `lib/pangenome_matrix.read_cluster_tsv`) to resolve each island's
  `member_families` down to the specific proteins belonging to that island's own
  `example_strain`, then joins against `gene_positions.tsv` (columns: `Short,
  protein_id, contig, start, end`, confirmed from
  `bin/pangenome_build_gene_positions.py`) to compute a min-start/max-end span.
  Adds columns: `locus_id` (`{example_strain}:{contig}:{start}-{end}`),
  `locus_contig`, `locus_start`, `locus_end`, `n_members_with_coordinates`,
  `n_contigs_in_locus`. Members with no resolvable coordinate (e.g. a rescue-pass
  `genome_only` call with a family-level tblastn position but no annotated
  `protein_id`) are excluded from the span and counted, not treated as an error.
  `n_contigs_in_locus > 1` is reported, not silently collapsed (probable
  paralog-copy pull-in — a known, documented limitation carried over from the
  original script). CLI: `--significant_islands --cluster_tsv --gene_positions
  --output`.
- Create: `modules/pangenome/island_locus.nf` — new `ISLAND_LOCUS` process,
  `input: path significant_islands, path cluster_tsv, path gene_positions`,
  `output: path("*.locus.tsv"), emit: annotated`. Always invoked (no gating param)
  when the whole island/Pfam step is enabled (`params.pangenome_island_pfam_hmm`
  set) — same top-level gate as the rest of this feature family, no separate flag.
- Modify: `workflows/pangenome_profile.nf` — insert `ISLAND_LOCUS` between
  `BUILD_ISLANDS` and everything downstream that currently consumes
  `BUILD_ISLANDS.out.significant_islands` (`SELECT_BACKGROUND_REPS` does NOT consume
  it — check at plan time which downstream consumers actually need the locus-
  augmented file vs. the original; `REPORT_TABLES` does). Add `gene_positions.tsv`
  as an input to this new process (it already exists in the workflow as
  `POSITIONS.out.gene_positions` or equivalent — confirm exact channel name at plan
  time).

**Testing:** unit tests matching the original `test_add_island_locus.py`'s coverage
(multi-contig flagging, unresolvable-member exclusion+counting, single vs. multi-
strain islands), ported and adapted to this repo's fixtures/column names.

## Component 4 — Per-strain outlier flag (always on, no new param)

**Files:**
- Modify: `bin/pangenome_report_tables.py` — `per_strain_summary()` (currently
  emits `Short, n_families, core, soft_core, shell, cloud, singleton` per strain)
  additionally computes, per bin column (`core, soft_core, shell, cloud,
  singleton`), the cross-strain mean/stdev and each strain's z-score, and adds one
  new column: `singleton_z` (z-score of the `singleton` count specifically — this
  is the bin the real `Asfu_H1106` case flagged; not computing z for every bin
  keeps the table narrow and matches what the real investigation actually needed).
  Add a second column `is_outlier` (`Y`/`N`), threshold `abs(singleton_z) > 3`
  (matches the "extreme outlier" convention; the real case was z=15.58, comfortably
  clearing this bar — a threshold this conservative won't false-positive on normal
  cross-strain variance). With fewer than ~3 strains, z-scores are undefined
  (stdev could be 0 or the sample too small to be meaningful) — in that case emit
  `-` for both new columns rather than a fabricated number or a crash (matches
  the small-cohort-fragility pattern already fixed in F8 of the prior plan).

**Testing:** unit test with a synthetic multi-strain fixture containing one clear
outlier (verify `is_outlier=Y` only for that strain), a fixture with no outliers
(verify all `N`), and a <3-strain fixture (verify `-` fallback, no crash).

## Wiring summary (workflows/pangenome_profile.nf)

```
BUILD_ISLANDS ──> ISLAND_LOCUS ──┐
                                  ├──> REPORT_TABLES ──> REPORT_RENDER
SELECT_BACKGROUND_REPS ──> FAMILY_PFAM_SCAN ──> DOMAIN_ENRICHMENT ──> [PFAM2GO if set] ──┘
```

`REPORT_TABLES`'s `per_strain_summary()` change (Component 4) is internal to that
process — no new channel wiring needed.

## New params (nextflow.config)

- `pangenome_pfam2go = null` — path to a local pfam2go mapping file; when unset,
  Component 2 is skipped entirely (no GO columns emitted, no error).

No other new params — Components 1, 3, 4 are unconditional once the existing
`pangenome_island_pfam_hmm`-gated step is enabled.

## Out of scope for this spec

- Sub-projects B (SwissProt), C (`cluster_homology_test.py`), D (stability tools) —
  separate specs.
- Any auto-fetch mechanism for the pfam2go file (explicitly rejected — operator-
  supplied only, per 2026-09-16 decision).
- Renaming/moving `lib/compressed_io.py` (explicitly deferred — user wants it to
  stay a candidate for genuine cross-pipeline sharing beyond just pangenome).
- Running this against real Coccidioides data (same boundary as the prior plan —
  the real-execution smoke test uses the existing 10-strain subset only).
