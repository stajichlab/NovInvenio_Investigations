# CLAUDE.md

This file guides Claude Code (and any other agent) working in this repository. See
`DESIGN.md` for the full design record (why this repo exists, what was decided, and the
still-open items). This file is the operational ruleset that follows from it.

## What this repo is

`NovInvenio_Investigations` (NII) holds the *science*: trait data, curated species/config
setup, per-study results, and analysis extensions (UniProt/GO-based inputs, functional
enrichment) built on top of the `nf_NovInvenio` pipeline (code lives there, not here —
see `DESIGN.md` Sec 2). NII depends on `nf_NovInvenio` as an external Nextflow pipeline,
invoked by reference (`nextflow run stajichlab/nf_NovInvenio`), never vendored/copied in.

Organized as one repo, domain subdirectories: `studies/<domain>/<set_name>/` (e.g.
`studies/fungi/pezizo_set1/`). See `DESIGN.md` Sec 3 for the full layout and Sec 8 for
the published-site structure this mirrors.

## Data provenance & tracking — hard rules

Three classes of data. Getting this wrong is the one thing this repo exists to avoid
repeating (`NovInvenio`'s own history had ~230 MB of regenerable report/data blobs
committed across ~200 commits before anyone decided that shouldn't happen).

1. **Ephemeral / recipe-driven** — anything under `data/` (gitignored). Raw UniProt
   reference-proteome pulls, NCBI genome/GFF3 pulls, `go-basic.obo`, any future GOA pull.
   **Never committed.** The pull script (`bin/fetch_uniprot_proteome.py`,
   `bin/fetch_genome_assembly.py`, etc.) is the checked-in artifact; its output is
   disposable and re-fetchable on demand. No versioning needed because nothing here is
   tracked.

2. **Curated / tracked** — `config_support/` (trait tables, curated sample/species
   lists), `studies/*/*/species.csv` + `config.csv` + `DATA_MANIFEST.yaml`. These live in
   git and **must** carry a provenance record. Also in this class, without a provenance
   record (they're generated, not sourced, but don't regenerate-and-grow the way class 3
   does): each study's `docs/<domain>/<set>/report.html` (landing page), `alignment.html`
   (static viewer shell), `archive/*.tsv.gz` (small derived tables), and the two gallery
   `index.html` levels. What makes something class 2 vs. class 3 is whether its size
   scales with candidate/sequence count — none of these do.

3. **Published / release-asset-only (extended 2026-09-11)** — the candidate/sequence
   -bearing parts of a study's `docs/<domain>/<set>/`: `novelties.html`, `core.html`,
   `losses.html`, `summary.pdf` (one real case reached 105MB, over GitHub's 100MB file
   limit), plus `alignments/`/`loss_alignments/` (real sequence text,
   multi-MB-gzipped per genome). **Never committed, regardless of size** — these
   regenerate on every study rerun, which is exactly the repeated-recommit pattern that
   caused this repo's original bloat incident. Publish via `bin/publish_report_release.sh`
   (report bundle) / `bin/publish_alignment_release.sh` (alignment shards) — both GitHub
   Release assets, `--clobber`-overwritten per rerun — instead of committing. `docs/`'s
   `report.html` (per-study landing page) and `alignment.html` (static viewer shell)
   stay in class 2 and are committed normally: neither grows with candidate data.
   `archive/*.tsv.gz` (small derived tables) and the two gallery `index.html` levels
   are also class 2. See `DESIGN.md` Sec 8's 2026-09-11 update for the full mechanism
   and the exact dividing line (grows with candidate/sequence count, or doesn't).

**No agent (Claude Code or otherwise) commits a new or changed tracked data file
without a provenance record for it.** If source URL/version/date is unknown, stop and
ask — never guess or omit. Mechanism: this repo uses Mycelium (`mycelium:ingest` +
`DATA_MANIFEST`); `lib/provenance.py` is the same idea applied directly inside the
pull scripts (`build_record()`/`write_record()`/`append_manifest()`), since those need
to run non-interactively as part of a study build, not just as an agent-driven ingest.

Every tracked data file's provenance record captures, at minimum:
- Exact source URL (the actual file pulled, not a landing page)
- Source version/release (e.g. `UniProt reference proteomes release 2026_02`)
- Access/download date
- Checksum (sha256) of the file it derives from
- License/redistribution terms
- If derived rather than raw: which script + parameters produced it (`derived_by`)

**Refresh ≠ silent overwrite.** Re-pulling from a newer upstream release bumps the
manifest's recorded version/date as a normal, visible commit — never a silent replace.
Re-running a study against a newer UniProt release is a new study run, not an in-place
refresh of the old one (`DESIGN.md` Sec 10).

## Running a study

```bash
bin/run_study.sh fungi/pezizo_set1 --run_tool diamond --pfam_hmm /path/to/Pfam-A.hmm
```

Builds `studies/<domain>/<set>/config.csv` + `data_dir/` from `species.csv` (via
`bin/build_study_config.py`, which drives the two fetch scripts per species) if not
already built, then invokes `nf_NovInvenio` against it. See `bin/run_study.sh`'s header
comment for the `NII_PIPELINE` override needed until the `NovInvenio` → `nf_NovInvenio`
rename actually happens.

## Where new code goes

| What | Where |
|---|---|
| A new external-data pull recipe (new database, new domain) | `bin/fetch_<source>.py`, following `fetch_uniprot_proteome.py`'s provenance pattern |
| Shared Python logic | `lib/` |
| A new study | `studies/<domain>/<set_name>/species.csv`, then `bin/build_study_config.py --study-dir studies/<domain>/<set_name>` |
| A script specific to one study (custom `species.csv` schema, hardcoded accession lists, a one-off model-organism pull, etc.) | `studies/<domain>/<set_name>/bin/`, not `bin/` — see "Study-specific vs. shared scripts" below |
| Enrichment/analysis scripts (GO/Pfam/InterPro ORA) | `bin/`, reading a study's annotated presence matrix — not yet built, see `DESIGN.md` Sec 7 |
| Site generation | `docs/` — release-asset published, never committed; see `DESIGN.md` Sec 8 |

### Study-specific vs. shared scripts

`bin/` is for scripts any study could call, parameterized by CLI args (`fetch_uniprot_proteome.py --proteome-id ...`, `build_study_config.py --study-dir ...`). A script belongs in `studies/<domain>/<set_name>/bin/` instead the moment it stops being reusable in that sense — the giveaway is usually right there in the name or the body:

- The filename itself names a study, species, or one-off dataset (`build_koxytoca_config.py`, `fetch_kpn78578_modelorg.py`) rather than a general source/action (`fetch_uniprot_proteome.py`).
- It hardcodes a species list, an accession list, or a specific study's directory layout instead of taking them as arguments.
- It exists to solve one study's particular data-shape problem (e.g. proteomes that are already local files instead of a UniProt/NCBI pull-per-species) rather than a recipe any study would reuse.

Study-specific scripts still follow the shared-`lib/` and provenance conventions (import `lib/provenance.py`, write `.provenance.yaml` sidecars / fold into that study's `DATA_MANIFEST.yaml`) — only their own location, and any hardcoded reference back to `NII_ROOT`/`NII_BIN`, changes. Since they don't live next to the scripts they invoke via `subprocess` (e.g. `fetch_uniprot_proteome.py`) or the `lib/` they import, hardcode absolute paths to those (`NII_ROOT = Path("/bigdata/.../NovInvenio_Investigations")`) rather than deriving them from `__file__` or cwd — same reasoning as `~/.claude/CLAUDE.md`'s BASH_SOURCE/SLURM warning, just triggered by directory distance instead of a SLURM work dir.

## Mycelium

This repo adopts the Mycelium living-repo framework (`.living/`, `todo/`,
`mycelium:ingest`), same as `nf_NovInvenio`. Confirmed lightweight: the plugin itself is
a shared, one-time install; per-project footprint is a handful of tracked markdown/YAML
files plus a gitignored `.mycelium/` run-state dir.
