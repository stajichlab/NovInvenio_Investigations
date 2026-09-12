# NovInvenio_Investigations (NII) — Design Plan

Status: draft, agreed in conversation on 2026-09-06. Nothing has been executed yet —
no `nf_NovInvenio` rename, no history rewrite, no UniProt pull, no repo creation on
GitHub. This document is the record of what was decided, so implementation can
proceed from here without re-litigating settled questions.

## 1. Purpose

Split NovInvenio into:
- **`nf_NovInvenio`** — the pipeline (code only), invoked as a standard Nextflow
  pipeline: `nextflow run stajichlab/nf_NovInvenio -r <tag>`.
- **`NovInvenio_Investigations` (NII)** — the science: trait data, curated species/config
  setup, per-study results, and the analysis extensions (UniProt/GO-based inputs,
  functional-enrichment reporting) that build on top of the pipeline.

Rationale: `NovInvenio`'s `.git` history had accumulated ~1.5 GB, of which ~230 MB
was real unique content — 206 MB in `view/` (self-contained HTML reports, re-committed
at full size on every regeneration) and 10.5 MB in `results/`. Config/data and pipeline
code were never cleanly separated, so every study's config, trait tables, and generated
reports lived in the same repo and same history as the pipeline logic. NII exists so
that doesn't happen again: the pipeline repo stays code-only and small; investigation
repos hold data, configs, and results and can grow/shrink/be pruned independently of the
pipeline's own history.

## 2. `NovInvenio` → `nf_NovInvenio` migration

**Decision: rename in place + purge history (git filter-repo), not a fresh repo.**
Rationale: repo is public but has 0 forks and 1 star — negligible external-dependency
risk — so a disruptive history rewrite is safe, and doing it properly once (rename +
purge) beats carrying 1.3 GB of dead weight forever or fragmenting into two repos with
no rename-redirect benefit.

Mechanics:
1. GitHub repo rename `NovInvenio` → `nf_NovInvenio` (Settings → rename). Auto-redirects,
   preserves all issues/PRs/stars.
2. `git filter-repo` to strip, from **all history**:
   - `view/` (206 MB, 136 blobs) — self-contained reports, fully regenerable.
   - `results/` (10.5 MB legacy tracked files) — same.
   - `pipeline/` (2 dead pre-Nextflow prototype scripts) — confirmed dead, drop entirely.
   - `stories/` — untrack from git (stop version-controlling), but leave the file
     (`stories/hex/Pyronema_query.fa`) on local disk; revisit later for validation use,
     not yet decided whether it's needed.
   - `db/modelorgs/*.csv` (7.4 MB, FungiDB gene-name lookups) — moves to NII (see §3),
     not investigation-agnostic pipeline code, drop from `nf_NovInvenio` entirely.
   - `configs/`, `config_support/` — moves to NII as current-state seed content (no
     history import needed there; NII starts its own history).
3. Repack (`git gc --aggressive`) after the filter-repo pass. Expected result: `.git`
   drops from 1.5 GB to low tens of MB.
4. Force-push the rewritten history. Every existing local clone (both `/bigdata` and
   `/rhome` copies) needs a fresh `git clone` afterward — old SHAs are gone.
5. Update `.github/workflows/docker-build.yml` and `nextflow.config`'s
   `container_version`/image-name references from `ghcr.io/stajichlab/novinvenio` to
   match the new repo name (exact new image name TBD at implementation time — not
   blocking this design).
6. What stays in `nf_NovInvenio`: `bin/`, `lib/`, `modules/`, `workflows/`, `main.nf`,
   `nextflow.config`, `conf/`, `pixi.toml`, `.github/`, `docs/`, `tests/` (including its
   existing tiny 4-species synthetic demo under `tests/data/` — already exactly the
   "small demo" a code-only pipeline repo needs), `.living/`, `todo/`, `skillpacks/`
   (pipeline-engineering memory, not scientific results — stays with the code repo).

## 3. NII repository structure

**Decision: one repo, domain subdirectories** (not one repo per domain). Rationale: the
UniProt/GOA-pull tooling, taxonomy→ingroup/outgroup logic, and GO/Pfam/InterPro ORA
scripts are identical regardless of domain — splitting into 4+ repos means duplicating
that tooling or building another shared-code dependency, re-creating the exact
entanglement problem just solved for `nf_NovInvenio`. One repo also gives one coherent
GitHub Pages gallery across all domains.

```
NovInvenio_Investigations/
├── DESIGN.md                        ← this file
├── CLAUDE.md / AGENTS.md             ← provenance rules (§4), mycelium conventions
├── bin/                              ← NII-specific tooling (not pipeline code):
│   │                                    UniProt pull scripts, .dat.gz → GO/Pfam/InterPro
│   │                                    extraction, ORA (goatools + custom hypergeometric)
│   └── ...
├── lib/                              ← shared Python for bin/ scripts
├── data/                             ← GITIGNORED. Ephemeral, recipe-driven pulls:
│   │                                    raw UniProt reference-proteome .fasta.gz/.dat.gz,
│   │                                    go-basic.obo, any GOA files if ever added later.
│   │                                    Never archived — regenerated on demand by bin/
│   │                                    scripts. Not backed up; reproducibility relies on
│   │                                    the provenance record (§4), not on keeping the
│   │                                    raw bytes around.
├── config_support/                   ← TRACKED. Moved from NovInvenio:
│   ├── traits/                          (funguild_reference.csv, trait_definitions.yaml,
│   │                                     traits.csv, README.md)
│   ├── animal_pool.csv, source_db.csv, 1KFG_*.csv, Chaetothyriales_samples.csv
│   └── modelorgs/                       (moved from NovInvenio's db/modelorgs/*.csv)
├── studies/
│   └── fungi/
│       ├── pezizo_set1/
│       │   ├── config.csv               ← NII-generated run config (UniProt-sourced)
│       │   ├── provenance.yaml          ← per-study provenance record (§4)
│       │   └── (pipeline outputs land here transiently; only small derived
│       │        artifacts get committed — see §7 archive/ )
│       └── (pezizo_set2, ... as they're defined)
│   ├── animal/       ← empty until scoped
│   ├── plant/
│   └── bacteria/
├── docs/                              ← GitHub Pages source (§8)
└── LICENSE                            ← CC-BY-4.0 (§8)
```

Configs previously in `NovInvenio/configs/` and `NovInvenio/config_support/` move here
as **current-state seed content** (copied, not history-imported — NII starts its own
git history). They'll largely be superseded by fresh UniProt-based configs per study,
but remain useful reference/comparison material.

## 4. Data provenance & tracking rules

Two classes of data, strictly enforced:

- **Ephemeral / recipe-driven** — raw UniProt reference-proteome files, `go-basic.obo`,
  any future GOA pull. Never committed (`data/` is gitignored). The *recipe* (pull
  script + its parameters) is the checked-in artifact; its output is disposable and
  re-fetchable on demand. No versioning needed because nothing is tracked.
- **Curated / tracked** — trait tables, curated sample/species lists, taxonomy↔proteome
  mapping tables, ingroup/outgroup definitions, and small derived summary tables meant
  for publishing (presence matrices, enrichment results). These live in git and **must**
  carry a provenance record.

Mechanism: **Mycelium** (`mycelium:ingest` skill + `DATA_MANIFEST` convention), the same
framework `NovInvenio` already uses. Confirmed lightweight enough to adopt: the plugin
itself is a 42 MB *shared, one-time* install (not per-project weight); per-project
footprint is tracked markdown/YAML (`.living/decisions.md`, `learnings.md`,
`conventions.md`, `todo/TODO_REGISTRY.md`) plus a gitignored `.mycelium/` run-state dir
and a handful of hook scripts already wired for exactly this (`mycelium-data-tracker.sh`
etc).

Every tracked data file's provenance record captures, at minimum:
- Exact source URL (the actual file path pulled, not a landing page)
- Source version/release (e.g. UniProt reference-proteome release `2026_02`, a specific
  `go-basic.obo` release date)
- Access/download date
- Checksum (sha256) of the upstream file it derives from
- License/redistribution terms (UniProt + GOA are CC BY 4.0)
- If derived rather than raw: which script + parameters + NII commit produced it

**Refresh ≠ silent overwrite.** Re-pulling from a newer upstream release bumps the
manifest's recorded version/date as a normal, visible commit — never a silent replace.

**Agent-facing hard rule** (goes in NII's `CLAUDE.md`/`AGENTS.md`): no agent commits a
new or changed tracked data file without first recording/updating its provenance entry.
If source URL/version/date is unknown, stop and ask — never guess or omit.

## 5. UniProt data acquisition

**Source: UniProt reference proteomes**
(`https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/`),
current release **2026_02** (10-Jun-2026), 3,528 Eukaryota reference proteomes at time
of writing.

Per-species pull is **two files**, both under
`reference_proteomes/<Domain>/<Proteome_ID>/`:
- `{Proteome_ID}_{taxid}.fasta.gz` — canonical protein sequences.
- `{Proteome_ID}_{taxid}.dat.gz` — SwissProt/TrEMBL flat-file format, and (confirmed by
  direct inspection of the *N. crassa* reference proteome, UP000001805) carries
  **everything else needed in this one file**:
  - `OC` lines — full taxonomic lineage (e.g. `Eukaryota; Fungi; Dikarya; Ascomycota;
    Pezizomycotina; Sordariomycetes; Sordariomycetidae; Sordariales; Sordariaceae;
    Neurospora.`) — **no separate NCBI taxonomy dump needed** for ingroup/outgroup
    construction.
  - `OX` line — NCBI taxon ID.
  - `GN` lines — gene names / ORF IDs.
  - `DR GO;` lines — GO terms with evidence codes (28,780 for *N. crassa* alone).
  - `DR Pfam;` / `DR InterPro;` lines — Pfam and InterPro cross-references (10,036 /
    26,095 for *N. crassa*) — **every protein already carries precomputed
    InterProScan-derived annotation**, including ingroup proteins that become novelty/
    loss candidates.

**Decision: GOA is not needed for v1.** The separate GOA GAF files (all of UniProtKB,
tens of GB) would only add more-current or broader annotation than what's baked into a
reference-proteome release — not a blocker for a first working pipeline. Revisit only
if a real annotation gap is found later.

**Decision: add a third pull for genome/DNA (resolves the Fable-review gap above).**
UniProt's per-proteome JSON endpoint (`https://rest.uniprot.org/proteomes/{Proteome_ID}.json`)
carries a `genomeAssembly.assemblyId` field — a standard GenBank accession (e.g.
`GCA_000182925.2` for *N. crassa*/UP000001805) — confirmed present for all 11
`pezizo_set1` species (§6). Genome + GFF3 pulled via `ncbi-datasets-cli` (bioconda,
confirmed available, `datasets download genome accession <GCA_...> --include genome,gff3`),
keyed off that same accession so genome and protein annotation are guaranteed
version-matched (this was the whole point — an assembly/annotation mismatch would
manufacture false TBLASTN "genome hit, no protein" signals). Same provenance treatment
as the other pulls (§4): recipe-driven, not archived, accession + access date recorded.

**Decision: no separate InterProScan run for v1.** Confirmed empirically that reference-
proteome `.dat.gz` files carry `DR Pfam`/`DR InterPro`/`DR GO` for essentially every
protein, including the ones that become candidates (they're UniProt entries from the
ingroup's own proteome). InterProScan-the-tool becomes a fallback capability for later,
needed only for sequences that aren't in UniProt at all (a genome not yet indexed there,
or genuinely novel gene models) or if per-match coordinate detail beyond the compact `DR`
lines is ever needed.

**Open / deferred: parser performance.** Biopython's `Bio.SwissProt` parser can read
`.dat.gz` directly, but its throughput on files of this size (multi-MB to tens-of-MB
compressed per species) is untested. Decision explicitly deferred until the pull +
extract framework exists to benchmark against; if Biopython is too slow, options are a
hand-rolled line-oriented parser for just the `AC`/`OC`/`OX`/`GN`/`DR GO`/`DR Pfam`/
`DR InterPro` lines actually needed (skipping the rest of the SwissProt-format record),
dumping the extracted per-protein summary into DuckDB for fast per-study querying rather
than re-parsing per report.

**New dependency not previously discussed:** `go-basic.obo`
(`https://purl.obolibrary.org/obo/go/go-basic.obo`, from geneontology.org, not UniProt) —
needed for GO-DAG structure/true-path propagation if `goatools` is to do full GO
enrichment properly. This needs the same provenance treatment as UniProt (§4): recipe-
driven pull, version/date recorded, not archived.

**2026-09-11 update — `species.csv`/`bin/build_study_config.py` generalized beyond
UniProt+NCBI-only.** This section (and §4) describe the UniProt/NCBI pulls this
repo started with; those remain the default path, but `bin/build_study_config.py`
now dispatches each species' protein/genome/GFF3 *independently* via its own
`Protein_Source`/`Genome_Source`/`GFF3_Source` `species.csv` columns — `uniprot`
(as above), `ncbi` (protein pulled from the same NCBI Datasets genome package,
for species with no UniProt reference proteome), `local_faa`/`local_genome`/
`local_gff3` (an already-downloaded file, copied in with a provenance record per
§4's rule, no fetch). This exists because bacteria studies (`UHM_Koxytoca`,
`UHM_lachnoNovelclade`) routinely combine sources this way — one species'
protein already sitting on disk, its genome still needing an NCBI pull. Full
rationale, schema, and migration record:
`notes/superpowers/specs/2026-09-11-study-onboarding-design.md`.

## 6. First study: `fungi/pezizo_set1`

Originally scoped as `pezizo5` minus Nirr/Mcir/Amega (8 species), but the independent
Fable-model review (§10 addendum) flagged that dropping those three leaves an all-Dikarya
outgroup with no polarizing power against Dikarya-ancestral gene loss — **reversed**:
`pezizo_set1` now replicates the full original `pezizo5` config, all 11 species, same
ingroup/outgroup split:

| Short | Species (strain) | Group | UniProt Proteome ID | Taxon ID | Assembly (GCA) |
|---|---|---|---|---|---|
| Nirr | *Neolecta irregularis* (DAH-3) | OUT | UP000186594 | 1198029 | GCA_001929475.1 |
| CneoH99 | *Cryptococcus neoformans* (H99) | OUT | UP000010091 | 235443 | GCA_000149245.3 |
| Ccin | *Coprinopsis cinerea* (Okayama-7/130) | OUT | UP000001861 | 240176 | GCA_000182895.1 |
| Mcir | *Mucor circinelloides* f. circ. (1006PhL) | OUT | UP000014254 | 1220926 | GCA_000401635.1 |
| Spom | *Schizosaccharomyces pombe* (972h) | OUT | UP000002485 | 284812 | GCA_000002945.2 |
| Scer | *Saccharomyces cerevisiae* (S288C) | OUT | UP000002311 | 559292 | GCA_000146045.2 |
| Amega | *Arthrobotrys megalospora* (TWF281) | IN | UP001658139 | 2528406 | GCA_036971655.1 |
| Ncra | *Neurospora crassa* (OR74A) | IN | UP000001805 | 367110 | GCA_000182925.2 |
| Afum | *Aspergillus fumigatus* (Af293) | IN | UP000002530 | 330879 | GCA_000002655.1 |
| Ztri | *Zymoseptoria tritici* (IPO323) | IN | UP000008062 | 336722 | GCA_000219625.1 |
| Cimm | *Coccidioides immitis* (RS) | IN | UP000001261 | 246410 | GCA_000149335.2 |

All 11 confirmed present in the UniProt Eukaryota reference-proteome set (release
`2026_02`), and all 11 resolve to a GenBank assembly accession via the proteomes JSON
endpoint (`genomeAssembly.assemblyId`) — no fallback-to-non-reference-proteome or
missing-assembly cases in this set.

Still open from the Fable review, not yet resolved (tracked in §9): Afum/Cimm
Eurotiomycete non-independence in the ingroup at `ingroup_min_frac 0.75`, and per-species
annotation-era heterogeneity (Ncra/Cimm ~2007 Broad, Ztri 2011 JGI, Ccin 2010 Broad vs.
continuously-curated Scer/Spom).

Further sets (`pezizo_set2`, etc.) to be defined the same way, guided per-set as this one
was.

### Additional `fungi/` sets (2026-09-07)

Three more pairwise IN/OUT contrasts, same UniProt-reference-proteome sourcing, same
`--run_tool diamond`/no-Pfam-hmmscan pattern as `pezizo_set1`:

- **`zoosporic_dikarya`** — IN: 4 zoosporic (chytrid) fungi (Batrdend, Spizpunc, Syncendo,
  Cateangu) vs OUT: 5 Dikarya (Scer, Ncra, Anid, Spom, Ccin). Originally scoped for 6
  zoosporic species (matching the old BFD-era `zoosporic_opisthokont_dikarya` config), but
  **Chytriomyces hyalinus and Paraphysoderma sedebokerense have no UniProt reference
  proteome** (confirmed against the Eukaryota reference-proteome listing — only
  non-reference entries exist for either) — dropped rather than breaking the
  reference-proteome-only sourcing principle.
- **`yeast_filamentous`** — IN: 4 yeasts, Saccharomycotina (Scer, Ylip) + Taphrinomycotina
  (Spom, Nirr), vs OUT: 4 filamentous Pezizomycotina (Ncra, Afum, Ztri, Cimm — reused
  directly from `pezizo_set1`'s ingroup). Direction chosen as "genes specific to yeast
  growth form, absent from filamentous relatives" — flip the GROUP column if the reverse
  question (filamentous-specific genes) turns out to be the one actually wanted.
- **`mushrooms_tremella`** — IN: 2 Agaricomycotina (mushroom-forming) species (Ccin,
  Agbis — *Agaricus bisporus* var. *burnettii*, the only UniProt-reference-proteome
  strain), vs OUT: 2 Tremellomycetes (CneoH99, Tremes — *Tremella mesenterica*). Thinnest
  of the three (2v2) — could be expanded with more Agaricomycotina/Tremellomycetes species
  later if the initial result looks background-noise-dominated.

All species/UniProt-proteome/GCA-accession mappings verified live against
`rest.uniprot.org` and the Eukaryota reference-proteome FTP listing before use, same
verification standard as `pezizo_set1`.

## 7. Functional enrichment

**Decision: ORA (over-representation analysis) only, hypergeometric test — no GSEA.**
Classic GSEA needs a continuous per-gene ranking (e.g. fold-change) and tests whether a
gene set skews toward one end of it via permutation; novelty/loss candidates are a
discrete set-membership call (candidate vs. not), which is what ORA is for. Applied to:

- **GO terms** — via `goatools` (handles GO-DAG propagation given `go-basic.obo`, and has
  a built-in hypergeometric/Fisher's-exact ORA implementation, so this is "use the
  library" rather than "build the stats").
- **Pfam domains** — hypergeometric test, custom (goatools is GO-specific), same
  candidate-vs-background framing.
- **InterPro domains** — same, custom hypergeometric.

## 8. Publishing

**2026-09-11 update — release-asset publishing now covers the whole *data-bearing* part
of a study's bundle, not just alignment shards; the small, stable parts stay committed.**
The rest of this section (structure diagram, archive format) still describes the site
layout accurately, but its original text below treated
`novelties.html`/`core.html`/`losses.html`/`summary.pdf` as fine to commit directly and
reserved the release-asset mechanism for `alignments/`/`loss_alignments/` alone. That
undersold the actual bloat risk: every study *rerun* recommits a fresh copy of these
(sequences dominate `novelties.html`'s size — one real case reached 105MB, over GitHub's
100MB file limit, 2026-09-10), not just once but on every iteration of a study while
it's being tuned — the exact repeated-bloat pattern this repo already hit once (the
~230MB/~200-commit history this repo's own retrospective cites). `.gitignore` already
reflects the corrected line (implemented ahead of this doc catching up):

- **Release-asset only, never committed:** `docs/*/*/novelties.html`, `core.html`,
  `losses.html`, `summary.pdf` (candidate-sequence-bearing, grows with the data),
  `alignments/`, `loss_alignments/` (real sequence text, multi-MB-gzipped per genome).
  Published via `bin/publish_report_release.sh` (report bundle) and
  `bin/publish_alignment_release.sh` (alignment shards) — both GitHub Release assets,
  `--clobber`-overwritten on every rerun (no run-history accumulation). A CI workflow
  (`.github/workflows/static.yml`) downloads both and assembles them into `docs/` only
  at Pages-deploy time — the committed repo never grows from either, no matter how many
  times a study reruns.
- **Still committed normally:** `report.html` (the small, stable per-study landing
  page) and `alignment.html` (the small static viewer shell) — neither grows with
  candidate data — plus `archive/*.tsv.gz` (small derived tables, thousands of rows,
  not multi-MB) and the two gallery `index.html` levels. `species.csv` / `config.csv` /
  `DATA_MANIFEST.yaml` (Sec 4's "curated/tracked" class) are unaffected either way —
  small, hand-curated, and don't regenerate on every pipeline rerun.

The dividing line is "does this file's size scale with candidate/sequence count," not
"is this file under `docs/`" — a landing page and a viewer shell don't, a report table
with embedded protein sequences does.

**Required one-time repo setting, not captured anywhere in code — GitHub Pages
"Build and deployment" source must be `GitHub Actions`, not `Deploy from a branch`.**
`.github/workflows/static.yml`'s `actions/deploy-pages` step only actually serves
traffic when Pages is configured this way; with the branch-based ("legacy") source,
GitHub builds and serves `docs/` directly from what's committed to `main`, ignoring
the workflow's runtime merge of `alignments-*`/`reports-*` release assets entirely —
the workflow reports success (a real deployment is created) while the live site
silently never shows anything beyond the small, always-committed `report.html`/
`alignment.html`. This exact failure happened 2026-09-09 through 2026-09-11: every
study's `novelties.html`/`core.html`/`losses.html`/`summary.pdf`/`alignments/` 404'd
on the live site the whole time, undetected because the workflow's own run history
showed nothing but green. Verify with:
```
gh api repos/stajichlab/NovInvenio_Investigations/pages --jq '.build_type'
# must print "workflow", not "legacy"
```
Fix (if it ever reverts — e.g. after a Pages setting reset or repo transfer):
```
gh api -X PUT repos/stajichlab/NovInvenio_Investigations/pages -f build_type=workflow
gh workflow run static.yml   # re-trigger so the fix takes effect immediately
```

**Companion rule, `nf_NovInvenio` side:** the pipeline source repo never holds analysis
run output at all, not even transiently pre-publish — see its own `CLAUDE.md`. A study
run's `--outdir` must resolve outside that checkout (`bin/run_study.sh` already does
this correctly: `--outdir "$NII_ROOT/results"`); an ad hoc run launched directly against
a `nf_NovInvenio` checkout with an in-repo `--outdir` also pollutes *that* repo's own
git-tracked `docs/` (its real, hand-written ADRs/Sphinx site — `Helpers.docsDir()`
resolves as a sibling of `--outdir`'s parent, so an in-repo `--outdir` lands there too,
not just the gitignored `results/`). Real incident, 2026-09-10/11: an agent-run
comparison study launched directly against the pipeline checkout did exactly this;
caught and cleaned up before anything was committed, but it's exactly the class of
mistake this two-repo split exists to make structurally hard, not just documented
against.

**Site structure** — two-tier gallery, extending `NovInvenio`'s `docs/<project>/`
publishing convention (`Helpers.docsDir()` — renamed from `view/` in 2026-09; the old
`view/index.html` + `view/generate_index.py` gallery tool this section originally
described has since been removed from `NovInvenio`) by one level, domain grouping above
the existing per-study grouping:

```
docs/
├── index.html                 ← top-level: one card per domain (Fungal, Animal, Plant,
│                                 Bacteria, Other — greyed out until populated)
├── fungi/
│   ├── index.html             ← domain landing: one card per set (pezizo_set1, ...);
│   │                             card = set name, one-line hypothesis, species count/
│   │                             ingroup-outgroup summary, date generated, link in
│   └── pezizo_set1/
│       ├── report.html        ← set landing (mirrors nf_NovInvenio's own
│       │                         docs/<project>/report.html)
│       ├── novelties.html
│       ├── core.html
│       ├── losses.html
│       ├── alignment.html             ← standalone TBLASTN alignment viewer
│       │                                 (nf_NovInvenio issue #86)
│       ├── alignments/                ← TBLASTN alignment shards (novelty direction) --
│       │                                 gitignored, NOT committed (see below)
│       ├── loss_alignments/           ← same, loss direction -- also gitignored
│       ├── go_enrichment.html         ← new
│       ├── domain_enrichment.html     ← new (Pfam + InterPro ORA)
│       └── archive/                   ← *.tsv.gz downloads for this set
├── animal/index.html          ← domain landing, empty-state card grid until sets exist
├── plant/index.html
└── bacteria/index.html
```

`studies/fungi/pezizo_set1/` (working dir) and `docs/fungi/pezizo_set1/` (published
site) share the same two-level key, so "where do I find X" stays consistent between repo
and live site.

**Archive format:** plain `.tsv.gz`, not zstd. Every tool (pandas, `zcat`, R) reads
gzip natively with zero extra dependency; the size difference on tables this small
(presence matrices, enrichment results — thousands of rows) doesn't justify requiring
`zstd` downstream. (This describes the *content* format only — per the 2026-09-11
update above, `archive/*.tsv.gz` is no longer committed either; it ships inside the same
per-study release asset as everything else.)

**`alignments/`/`loss_alignments/` were the first directories moved to release-asset
publishing, not the only ones (extended 2026-09-11 above).** They were the forcing
case: TBLASTN alignment shards carry real sequence text and can run into the
multi-MB-gzipped range per genome, so they compound the fastest per study rerun — the
`view/`-style repo-bloat mistake `NovInvenio` itself already made once (see the
retrospective note below), which this repo exists partly to avoid repeating. Once that
mechanism existed, extending it to `novelties.html`/`core.html`/`losses.html`/
`summary.pdf` (same recompound-per-rerun problem, just a slower fuse — one real instance
still crossed GitHub's 100MB file limit) was the smaller step, not a separate decision —
see the 2026-09-11 note above for the exact dividing line (`report.html`/
`alignment.html`/`archive/*.tsv.gz` stay committed; they don't grow with candidate
data). `bin/publish_alignment_release.sh` (shards) and `bin/publish_report_release.sh`
(report bundle) both push via `--clobber` and are downloaded and merged back into
`docs/` only at Pages-deploy time (`.github/workflows/static.yml`) — the committed repo
never grows from either no matter how many times a study reruns.

**License:** CC-BY-4.0 for data/pages (matches UniProt/GOA upstream terms, requires
attribution, permits redistribution/adaptation). Every published page carries a standing
attribution line ("Includes data from UniProt (release 2026_02),
https://www.uniprot.org"). Code (any scripts in `bin/`/`lib/`) may warrant a separate,
more permissive license (MIT/BSD) — not yet decided, low priority.

## 9a. Branding assets

Logo (`NI_icon.png`, 320×320, RGBA) originated in `NovInvenio/assets/` (untracked).
Canonical copy now lives in NII at `assets/logo/NI_logo_source.png`, since NII's `docs/`
site is what actually renders card grids/favicons and gets the most use out of it.
Generated sizes (ImageMagick, all in `assets/logo/`):

| File | Size | Use |
|---|---|---|
| `NI_logo_source.png` | 320×320 | canonical original |
| `NI_logo_512.png` | 512×512 | header logo, top-level/domain landing pages |
| `NI_logo_apple-touch-180.png` | 180×180 | `apple-touch-icon` |
| `NI_logo_card-96.png` | 96×96 | per-card icon on domain/set gallery cards |
| `NI_logo_favicon-32.png` / `-16.png` | 32×32 / 16×16 | favicon (PNG) |
| `NI_logo_favicon.ico` | 16/32/48 multi-res | favicon (`.ico`, older browsers) |

Duplicated (not URL-referenced) into `nf_NovInvenio`'s own `assets/logo/` too — same
files, same names, in both repos — since `nf_NovInvenio`'s README/docs should render the
logo without depending on NII's repo being reachable/public.

## 9. Explicitly deferred / open items

- ~~Parser performance~~ **Resolved.** `bin/extract_dat_annotations.py` is a hand-rolled,
  single-pass parser (only AC/GN/OX/DR lines, resets on `//`) rather than Biopython's
  full SwissProt parser. Benchmarked against the real *N. crassa* `.dat.gz` (9,759
  proteins, 9.7 MB compressed): **1.75s wall time.** Fast enough that DuckDB storage
  isn't needed for a single proteome; revisit only if a study's combined per-species
  TSVs get large enough that repeated re-parsing across enrichment runs becomes the
  bottleneck (not the case yet).
- GOA integration — only if a real annotation-completeness gap is found.
- Full InterProScan execution — fallback for non-UniProt inputs only, not built for v1.
- `stories/` folder — untracked from git per §2, kept on local disk, revisit later for
  possible validation use.
- Timing/scope for Animal/Plant/Bacteria domains — Fungal (`pezizo_set1`, and further
  `pezizo_setN`) goes first, end-to-end, before replicating the pattern elsewhere.
- **TBLASTN cost at Animal/Plant genome scale.** Fungal genomes are ~30-40 Mb; the
  `nf_NovInvenio` TBLASTN validation step runs a translated search per candidate-cluster
  representative against *each outgroup genome*. Animal/plant genomes are 1-30+ Gb
  (human ~3.1 Gb, many plants 5-15+ Gb) — 2-3 orders of magnitude larger per genome, and
  the domain will also have more outgroup genomes than a typical fungal study. Before
  scaling this pattern to Animal/Plant, investigate whether full-genome TBLASTN
  validation is still worth its cost there, versus alternatives (validate against a
  representative subset of outgroup genomes rather than all of them, a faster translated
  search such as DIAMOND `blastx`/`mmseqs`, or dropping genomic validation and relying on
  the proteome-level presence/absence call alone for these domains). Not resolved — flag
  and revisit when Animal/Plant scoping actually starts.
- Exact new container image name for `nf_NovInvenio`'s Docker build (currently
  `ghcr.io/stajichlab/novinvenio`) — implementation detail, not blocking.
- `ncbi-datasets-cli` (bioconda, confirmed available, v18.36.0) needs adding to NII's own
  `pixi.toml` for the genome/GFF3 pull in §5 — not a `nf_NovInvenio` dependency, since the
  pipeline itself stays agnostic to where its input FASTAs came from.
- UniProt's per-proteome JSON completeness field (`proteomeCompletenessReport`) exists
  but the exact BUSCO-score JSON path wasn't nailed down during this session — worth
  finding before implementing the "record completeness in provenance.yaml" recommendation
  from the Fable review.
- **Phase 2: open-ended clade/taxon discovery — deferred.** Tasks 1-4 implement Phase 1
  only: given a species already known to belong to a study, auto-discover its UniProt
  proteome ID and NCBI genome accessions via existing taxonomy/nomenclature. Phase 2
  would be open-ended "find me good Ascomycota outgroups I haven't heard of" automation
  via UniProt clade/phylogeny queries or other discovery heuristics — valuable long-term,
  but a separate scope. See `notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md`'s
  own "Scope: Phase 1 only" section for the reasoning.

## 10. Expert evaluation

### As a bioinformatics reviewer

1. **The ORA background set is the single biggest determinant of what "enriched" means,
   and it isn't decided yet.** Candidate-vs-what? Whole ingroup proteome, whole presence
   matrix (ingroup+outgroup), or just the other candidates? These give materially
   different, sometimes contradictory, results. Recommend: background = every protein in
   the presence matrix for that direction (matches what `make_novelties.py` already
   considers the population) — but this needs to be pinned per report, not left implicit,
   or two people will get two different enrichment tables from the "same" candidate list.
2. **Family/paralog inflation.** `mmseqs`-clustered candidate families already exist in
   this pipeline (`clusters_cluster.tsv`); if ORA counts every individual candidate
   *protein* rather than one representative per family, a single expanded gene family
   with one shared Pfam domain can look like "significant enrichment" when it's actually
   one evolutionary event. Count at the family level, not the protein level.
3. **GO-term multiplicity.** Testing GO + Pfam + InterPro ORA in parallel triples the
   multiple-testing burden across three highly *correlated* annotation vocabularies
   (Pfam domains often map directly to specific GO terms via Pfam2GO) — expect redundant
   "hits" across all three; apply BH/FDR correction (goatools does this for GO; the
   custom Pfam/InterPro tests need the same) and report all three together so a reader
   can see the overlap, not as three independent claims.
4. **Reference-proteome churn.** UniProt reference proteomes update every ~8 weeks;
   re-running `pezizo_set1` against a later release can shift accessions/isoform calls
   even for the "same" genome. The provenance pinning in §4 covers this, but worth
   restating: a study's *config* should record the exact release, and re-running against
   a newer release is a new study run, not an in-place refresh.
5. **Strain/proteome-selection rule for future sets.** `pezizo_set1`'s 11 species all had
   an unambiguous single reference proteome to pick, but note *Mucor circinelloides*
   itself (now included, §6) has 4-5 competing UniProt proteome entries for closely
   related strains/subspecies — the pick was made by matching the exact strain
   (1006PhL) the prior BFD-based `pezizo5` config used. Needs this same explicit,
   repeatable selection rule (prefer the UniProt "reference" flag, then match the known
   prior strain) applied consistently as sets with messier taxonomy get added.

### As a data-management reviewer

1. **"Ephemeral, not archived" is in tension with reproducibility**, not fatally, but
   worth naming: UniProt doesn't guarantee indefinite retention of old releases, so a
   from-scratch reproduction of `pezizo_set1` five years from now may not be able to
   re-fetch byte-identical `2026_02` files. The provenance record (checksums, release ID,
   accessions actually used) is what makes the *derived*, tracked artifacts the durable
   source of truth even if raw re-fetch ever becomes impossible — that's already the
   design, just worth being explicit that "not archived" means "not bit-reproducible
   forever," and that's an accepted tradeoff, not an oversight.
2. **Don't let NII repeat `NovInvenio`'s exact mistake.** The whole reason this split
   exists is that `view/`-style self-contained reports quietly bloated a repo to 1.5 GB
   over ~200 commits with nobody deciding that should happen. NII's `docs/archive/`
   pattern is lower-risk (small compressed TSVs, not multi-MB HTML with embedded
   sequences) but will still compound across studies × domains × reruns. Worth a
   lightweight periodic size check (even just `du -sh .git` in a session-start hook)
   rather than rediscovering this the same way, later, at larger scale.
3. **Mycelium hook scope.** The per-tool-call hooks (`mycelium-read-tracker.sh` etc.)
   should explicitly exclude `data/` (the gitignored ephemeral UniProt/GOA cache) —
   parsing multi-MB `.dat.gz` files repeatedly during development shouldn't generate
   tracking noise for content that's deliberately untracked and disposable.
4. **New external dependency needs the same discipline as UniProt.** `go-basic.obo`
   (§5) comes from a different source (geneontology.org, not UniProt) and wasn't part of
   the original UniProt/GOA framing — make sure it gets a provenance record too, not just
   UniProt-derived files, or the provenance rule in §4 has a silent gap on day one.
5. **License diligence is a one-line check, not an assumption.** CC-BY-4.0 for UniProt/
   GOA is confirmed; Pfam-A/InterPro data specifically (surfaced via UniProt's own `DR`
   cross-references, not a separate pull) should get the same one-line verification
   before the first public GitHub Pages publish, rather than assuming it inherits
   UniProt's license by association.

### Independent second-opinion review (Fable model, 2026-09-06)

Requested as a critical, read-only pass over this document, explicitly told not to
repeat §10 above. Flagged several things §10 missed:

**Input-source biases (UniProt reference proteomes vs. BFD assemblies)**
- **Missing genome FASTA — a real gap, not yet resolved.** §5's per-species pull
  (`.fasta.gz` + `.dat.gz`) has no DNA. TBLASTN validation (and GFF3 chrom/start) need
  the genome assembly, which isn't in either file. Needs a third pull (ENA/NCBI
  assembly, version-matched to the annotation the UniProt proteome was built from) with
  its own provenance record — an assembly/annotation version mismatch would manufacture
  exactly the "genome hit, no protein call" signal TBLASTN exists to catch. **See open
  question below — this needs a decision, not just a note.**
- **Annotation heterogeneity now varies by species**, confounded with ingroup/outgroup
  membership (Ncra/Cimm ~2007-era Broad gene models, Ztri 2011 JGI, Ccin 2010 Broad, vs.
  continuously curated Scer/Spom). A missed gene model in one outgroup's annotation
  reads as a "novelty." Recommend recording each proteome's UniProt-reported BUSCO/CPD
  completeness in `provenance.yaml`, and treating TBLASTN as a hard filter (not
  reporting-only) for this study specifically, given the heterogeneity.
- Reference proteomes include `Fragment`/`NON_TER` entries — filter or flag these, they
  inflate both spurious novelties and losses.
- Canonical-isoform-only fixes the isoform-inflation problem, but UniProt's canonical
  pick isn't always the longest isoform — results won't be directly length-comparable to
  the old `pezizo5` run.
- Once IDs are UniProt accessions, the SwissProt-diamond annotation step becomes a
  self-hit, and `GN` lines supersede `config_support/modelorgs/*.csv` for these species —
  need one clear `function_source` priority order, not multiple overlapping sources.

**Species choices for `pezizo_set1` — outgroup polarization concern**
- Dropping Nirr, Mcir, and Amega leaves an all-Dikarya outgroup (two genome-reduced
  yeasts + two Basidiomycota). A Dikarya-ancestral gene lost in both basidiomycetes and
  both yeasts is indistinguishable from a true Pezizomycotina gain. Neolecta
  (non-reduced, filamentous Taphrinomycotina) was the most polarizing outgroup in the
  original 11-species set; its removal should be reconsidered, or a non-Dikarya outgroup
  restored.
- Afum + Cimm are both Eurotiomycetes; at `ingroup_min_frac 0.75` (3/4), a gene present
  in both Eurotiomycetes plus one other already clears the "Pezizomycotina-wide" bar —
  ingroup isn't taxonomically independent as currently chosen.

**Enrichment statistics — underspecified in ways that would bias results**
- Novelty candidates are systematically short/fast-evolving/low-complexity
  (homology-detection failure, not real absence — Weisman et al. 2020); ORA on such a
  set mostly reports depletion of everything conserved unless this is controlled for.
- Candidates are disproportionately unannotated by construction — the ORA universe must
  be restricted to proteins/families with ≥1 term in the tested vocabulary, and
  per-species annotation depth varies a lot (Scer/Ncra >> Ztri/Cimm/Ccin), which will
  dominate a pooled cross-species result.
- §10's "count at the family level" fix needs the **background** clustered the same way
  as the candidates, not just the candidate set — otherwise orthologs are 4× in the
  universe and 1× in the candidate set.
- TrEMBL `DR GO` is almost entirely `IEA` (electronic, InterPro2GO/UniRule-derived) —
  i.e. GO is largely *derived from* Pfam/InterPro, which is the actual mechanism behind
  §10's GO/Pfam/InterPro redundancy point. Decide whether IEA evidence is included
  (excluding it leaves Ztri/Cimm near-zero GO) and state the choice per report.
- Needs pinning: one-sided (over-representation) test, `propagate_counts=True`, BP/MF/CC
  tested separately, a minimum term size, and how InterPro parent/child hierarchy + Pfam
  clans get collapsed (otherwise nested InterPro entries register as multiple
  independent hits). Per-species vs. pooled-candidate-set ORA is also still undecided.
- ORA-only (no GSEA) reconfirmed as the right call — no natural continuous ranking
  exists here for v1.

**Precomputed DR annotations**
- `.dat.gz` doesn't state which InterPro/Pfam *release* it embeds — record that
  alongside the UniProt release in provenance, so DR-derived Pfam IDs can be reconciled
  against the pipeline's own bundled `Pfam-A.hmm` version if both are ever used together.
