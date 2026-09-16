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
  bugs invisible to a `--help`-only check). Note: the prior plan's smoke test was
  built ad hoc under `/scratch/.../island_pfam_smoketest*` during its own
  execution, not committed to the repo as a reusable fixture — a repo-wide grep
  found no 10-strain fixture or launcher checked in. The plan for this spec must
  name the exact samplesheet/data-dir path it will use (rebuild it the same way, or
  locate wherever it was left) rather than assume one already exists.
- Prefer `bin/pangenome_*.py` naming (repo convention, 100% consistent today per the
  2026-09-16 bin/lib reorg survey — do not deviate).
- `lib/compressed_io.py` stays unprefixed and untouched by this work (explicit user
  decision 2026-09-16: it's intended to become genuinely cross-pipeline shared, not
  pangenome-only, revisit separately).

## Component 1 — InterPro hotlinks (always on, zero new params)

**Files:**
- Modify: `bin/pangenome_domain_enrichment.py` — add a shared helper
  `bare_pfam_accession(accession: str) -> str` (strips the version suffix, e.g.
  `PF13577.9` → `PF13577`; passes through `-` unchanged). `pfam_accession` itself
  stays versioned in `island_pfam_enrichment.tsv` (unchanged — do not break the
  existing test at `tests/test_pangenome_domain_enrichment.py:96`, which asserts
  `data_row[1] == "PF13577.9"`). Add a new `pfam_url` column immediately after
  `pfam_accession`, built from `bare_pfam_accession()`'s output:
  `https://www.ebi.ac.uk/interpro/entry/pfam/<bare_accession>/`, or `-` when the
  accession itself is `-`. Update the existing header-assertion test deliberately
  to include `pfam_url` in the expected column list.
- Modify: `bin/pangenome_report_render.py` — `render_report_markdown()`'s "Pfam
  domain enrichment" table (`render_report_markdown()`, exact lines 331-339) renders
  `[{domain}]({row['pfam_url']})` instead of a bare domain name when
  `row.get("pfam_url")` is present and not `-`; falls back to the bare name
  otherwise (use `.get()`, not direct indexing — existing tests build `top_domains`
  rows by hand without a `pfam_url` key, e.g. `tests/test_pangenome_report_render.py:15`,
  and must keep passing unchanged). No change to `plot_domain_enrichment()` (figures
  don't carry hyperlinks).

**No new Nextflow process, no new param.** This is an in-place augmentation of an
existing script's output plus a render-time formatting change.

**Testing:** unit test for `bare_pfam_accession()` (version-suffix stripping,
`-` passthrough), unit test for the render function choosing link vs. bare-name
based on column presence/value (including the no-`pfam_url`-key case, to guard the
`.get()` fallback).

## Component 2 — Pfam2GO annotation (gated on `--pangenome_pfam2go`)

**Files:**
- Create: `bin/pangenome_pfam2go.py` — ported from
  `studies/fungi/Afumigatus_pangenome/bin/map_pfam_to_go.py`. Parses a standard
  `pfam2go` file (format: `Pfam:PFxxxxx <name> > GO:<term> ; GO:nnnnnnn` lines) into
  `{bare_pfam_accession: [go_id, ...]}` (bare, unversioned keys — the pfam2go file
  itself has no version suffixes). Imports `bare_pfam_accession()` from
  `pangenome_domain_enrichment.py` (same helper Component 1 uses — do not
  reimplement the stripping logic) to convert each input row's versioned
  `pfam_accession` before the lookup. Reads `island_pfam_enrichment.tsv`
  (Component 1's output — has both `pfam_accession` and `pfam_url` already), and
  writes the same table with two added columns: `n_go_terms`, `go_terms`
  (`;`-joined bare GO IDs, no spaces, e.g. `GO:0004930;GO:0005886`; drops the GO
  term names the original Afumigatus script kept — accession-only, matching this
  pipeline's existing "accession is the stable key" convention — or `-` if none
  mapped). CLI: `--island_pfam_enrichment --pfam2go --output`.
- Create: `modules/pangenome/pfam2go.nf` — new `PFAM2GO` process, `input: path
  island_pfam_enrichment, path pfam2go`, `output: path("island_pfam_enrichment.go.tsv"),
  emit: annotated` (fixed filename, matching every sibling process's convention —
  not a glob). Only invoked when `params.pangenome_pfam2go` is truthy.
- Modify: `workflows/pangenome_profile.nf` — when `params.pangenome_pfam2go` is
  set, define `enrichment_for_report = PFAM2GO(DOMAIN_ENRICHMENT.out.enrichment,
  file(params.pangenome_pfam2go)).out.annotated`; otherwise
  `enrichment_for_report = DOMAIN_ENRICHMENT.out.enrichment`. Feed
  `enrichment_for_report` to **both** `REPORT_TABLES` (its existing
  `--island_pfam_enrichment` input, currently a dead/unused argparse arg per
  `bin/pangenome_report_tables.py:122` — note this as a pre-existing minor
  cleanup, not something this spec fixes beyond keeping it fed with the right
  data for whenever it is used) **and** `REPORT_RENDER`, which currently receives
  `DOMAIN_ENRICHMENT.out.enrichment` directly at `workflows/pangenome_profile.nf:263`
  — that call site must be changed to `enrichment_for_report`, or the GO columns
  never reach the rendered report at all. This is the same if/else
  channel-reassignment pattern already used for `rescued_matrix`
  (`workflows/pangenome_profile.nf:104-111`), not the marker-search branch's
  `Channel.value([])` pattern.
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
  absent (matches Component 1's presence-conditional pattern, and the same
  `row.get(...)` defensive-read requirement).

**Testing:** unit tests for `bin/pangenome_pfam2go.py`'s parser (real pfam2go-format
fixture, multiple GO terms per accession, an accession with none, and one whose
input `pfam_accession` is versioned to confirm `bare_pfam_accession()` is actually
used for the lookup) and its join/output logic; workflow-level test via the
real-execution smoke test with and without `--pangenome_pfam2go` set, explicitly
checking `report.md` for the GO column (not just the intermediate TSV) in the
with-pfam2go run.

## Component 3 — Island genomic locus (always on, no new param)

**Design note (revised from an earlier draft of this spec after review):** rather
than a new standalone process, this folds into `REPORT_TABLES`, which already does
the equivalent kind of per-island join (`annotate_islands_with_domains()`) against
the same `islands_with_domains.tsv` output. `BUILD_ISLANDS.out.islands` (the real
emit name — not `significant_islands`) has exactly two consumers today:
`DOMAIN_ENRICHMENT` (`workflows/pangenome_profile.nf:247`, reads only
`member_families`, does not need locus columns) and `REPORT_TABLES`
(`:250`). Leaving `DOMAIN_ENRICHMENT` on the unmodified `BUILD_ISLANDS.out.islands`
(no change there) and adding the locus computation only inside `REPORT_TABLES`
avoids a third near-duplicate published TSV and one more process on the critical
path, without touching `DOMAIN_ENRICHMENT` at all.

**Files:**
- Modify: `bin/pangenome_report_tables.py` — add a new function
  `add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")`
  (logic ported from `studies/fungi/Afumigatus_pangenome/bin/add_island_locus.py`),
  called from `main()` after `annotate_islands_with_domains()` and before writing
  `islands_with_domains.tsv`. Inverts `cluster_tsv` (already loaded via the
  existing `lib/pangenome_matrix.read_cluster_tsv`, same as
  `bin/pangenome_build_islands.py:213`'s precedent) to resolve each island's
  `member_families` down to the specific proteins belonging to that island's own
  `example_strain`, splitting on `id_sep` (default `"|"`, **not hardcoded** — see
  the `--id_sep` CLI addition below, matching the F6 fix already applied to
  `pangenome_build_islands.py`/`lib/pangenome_matrix.py`), then joins against
  `gene_positions.tsv` (columns: `Short, protein_id, contig, start, end`) to
  compute a min-start/max-end span. Adds columns: `locus_id`
  (`{example_strain}:{contig}:{start}-{end}`, or `-` if unresolvable — matching
  this repo's existing sentinel convention, not the original script's `""`),
  `locus_contig`, `locus_start`, `locus_end` (each `-` when unresolvable),
  `n_members_with_coordinates`, `n_contigs_in_locus`. Members with no resolvable
  coordinate (e.g. a rescue-pass `genome_only` call with a family-level tblastn
  position but no annotated `protein_id`) are excluded from the span and counted,
  not treated as an error. `n_contigs_in_locus > 1` is reported, not silently
  collapsed (probable paralog-copy pull-in — a known, documented limitation
  carried over from the original script).
- Modify: `bin/pangenome_report_tables.py`'s `main()` — add two new required CLI
  args, `--cluster_tsv` and `--gene_positions`, plus `--id_sep` (default `"|"`).
  The zero-island fallback header (currently a hardcoded list at
  `bin/pangenome_report_tables.py:143-146`) must be extended to include the 6 new
  locus columns, so a zero-island run's `islands_with_domains.tsv` has the same
  header shape as a non-zero run.
- Modify: `modules/pangenome/report.nf` — `REPORT_TABLES`'s `input:` block gains
  `path cluster_tsv, path gene_positions`; its script block passes them plus
  `--id_sep '${params.pangenome_id_sep}'` (the pre-existing param, default `'|'`,
  already defined at `nextflow.config:192` — no new param needed here).
- Modify: `workflows/pangenome_profile.nf` — `REPORT_TABLES`'s call gains two new
  args: `CLUSTER_TIER1.out.cluster_tsv` (existing emit) and
  `GENE_POSITIONS.out.positions` (existing emit — the real process/emit name;
  not `POSITIONS.out.gene_positions`).
- Modify: `bin/pangenome_report_render.py` — in the existing "## Accessory
  islands" section (`render_report_markdown()`, after the count/figure lines,
  currently ~line 307-312), add a "top islands" markdown table (columns:
  `locus_id`, `island_size`, `n_strains`, `pfam_domains`; capped at the top 20 by
  `island_size` descending, matching `plot_domain_enrichment()`'s existing
  `top_n=20` convention) — rendered only when at least one row has a `locus_id`
  that is not `-`; omitted entirely otherwise (this is the report surface for the
  locus data — without it, `locus_id` would only ever appear in the TSV, which
  defeats the "citable coordinates" motivation).

**Testing:** unit tests matching the original `test_add_island_locus.py`'s coverage
(multi-contig flagging, unresolvable-member exclusion+counting, single vs. multi-
strain islands, `-` sentinel on total failure, custom `--id_sep` correctly
resolving members), ported and adapted to this repo's fixtures/column names; a
render test for the new top-islands table's presence/absence conditional.

## Component 4 — Per-strain outlier flag (always on, no new param)

**Design note (revised from an earlier draft after review):** a plain
mean/stdev z-score with threshold `abs(z) > 3` is mathematically unreachable for
small cohorts — for a sample of size n, the maximum possible `|z|` for any single
point is bounded (`sqrt(n-1)` with population stdev, `(n-1)/sqrt(n)` with sample
stdev), which equals exactly 3.00/2.85 at n=10. This pipeline's mandated
real-execution smoke test uses a 10-strain subset, so a fixed `>3` cutoff on
mean/stdev could **never** fire there or on any similarly-sized real study,
regardless of how extreme a strain actually is. The real `Asfu_H1106` case (the
motivating precedent) was also computed against the panel **median**, not the
mean, over 295 strains — the original investigation's own method, not a mean/stdev
z-score. Use a robust statistic instead: the modified z-score
(`0.6745 * (x - median) / MAD`, Iglewicz-Hoaglin), which is not bounded by n the
same way and matches the median-based reasoning the original investigation
actually used.

**Files:**
- Modify: `bin/pangenome_report_tables.py` — `per_strain_summary()` (currently
  emits `Short, n_families, core, soft_core, shell, cloud, singleton` per strain)
  additionally computes, across all strains' `singleton` counts only (not every
  bin — this is the bin the real `Asfu_H1106` case flagged, and keeping the table
  narrow matches what the real investigation actually needed): the median, the
  median absolute deviation (MAD), and each strain's modified z-score
  `0.6745 * (singleton_count - median) / MAD`. Adds column `singleton_z` (rounded
  to 2 decimals) and `is_outlier` (`Y`/`N`, threshold `abs(singleton_z) > 3.5` —
  the standard Iglewicz-Hoaglin cutoff). Exact fallback rules, both required
  (state explicitly, do not leave to implementer judgment):
  - Fewer than 3 strains: emit `-` for both columns (a modified z-score is not
    meaningful below this; matches the F8 small-cohort-fragility precedent).
  - MAD == 0 (most strains share the same singleton count, so the statistic is
    undefined/would divide by zero): emit `-` for both columns rather than `inf`
    or a crash.
- Modify: `bin/pangenome_report_render.py` — in the existing "## Accessory
  islands" (or a location the plan should confirm against the current section
  layout) area, add one templated line when at least one strain has
  `is_outlier == "Y"`: `**Outlier strains (singleton-count modified z-score
  beyond threshold):** <comma-joined Short names>` — this is templated data (a
  literal listing driven by the TSV), not hand-written narrative, matching the
  existing convention used for e.g. the classification-breakdown bullet list.
  Omit the line entirely when no strain is flagged (no "none found" filler).

**Testing:** unit test with a synthetic multi-strain fixture (>=3 strains)
containing one clear outlier (verify `is_outlier=Y` only for that strain, and that
the fixture size is actually large enough for the modified-z statistic to fire —
do not repeat the earlier mean/stdev draft's mistake of an untested-at-scale
threshold), a fixture with no outliers (verify all `N`), a <3-strain fixture
(verify `-` fallback, no crash), and a fixture where all strains share the same
singleton count (verify the MAD==0 fallback, no crash/`inf`).

## Wiring summary (workflows/pangenome_profile.nf)

```
BUILD_ISLANDS.out.islands ──────────────────────┬──> DOMAIN_ENRICHMENT (unchanged)
                                                 │         │
CLUSTER_TIER1.out.cluster_tsv ──────────────────┤         v
GENE_POSITIONS.out.positions ───────────────────┤   DOMAIN_ENRICHMENT.out.enrichment
                                                 │         │
                                                 │   [PFAM2GO if params.pangenome_pfam2go]
                                                 │         │
                                                 v         v
                                          REPORT_TABLES <──┴── enrichment_for_report
                                                 │
                                                 v
                                          REPORT_RENDER <── enrichment_for_report (also, directly)
```

`enrichment_for_report` (new channel variable) is `PFAM2GO.out.annotated` when
`params.pangenome_pfam2go` is set, else `DOMAIN_ENRICHMENT.out.enrichment`
unchanged — fed to **both** `REPORT_TABLES` and `REPORT_RENDER` (the latter's call
site currently takes `DOMAIN_ENRICHMENT.out.enrichment` directly and must be
updated, or Component 2's GO columns never reach the rendered report). Component
3 (island locus) and Component 4 (outlier flag) are both internal to
`REPORT_TABLES`/`REPORT_RENDER` — Component 3 needs two new inputs
(`cluster_tsv`, `gene_positions`) threaded into `REPORT_TABLES`; Component 4 needs
no new channel wiring at all.

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
