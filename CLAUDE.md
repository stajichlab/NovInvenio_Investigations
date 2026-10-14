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
   **Since 2026-09-24 every study publishes per run, at `docs/<domain>/<study>/<run>/`**
   (spec `notes/superpowers/specs/2026-09-24-study-run-site-layout-design.md`);
   a folder is published only if it has `studies/<domain>/<folder>/publish.yaml`.
   The class 2/3 split above applies one level down. Set-level
   `novelties/core/losses.html` may only be committed as old-URL redirect stubs.
   pangenome.nf runs (`docs/<domain>/<set>/<run>/`, staged by
   `bin/sync_pangenome_report.py`): only the run's `report.html` + `run.json` are
   committed; `figures/`, `figures_pdf/`, `archive/`, `clinker/`, `island_synteny.html`,
   `assembly_quality.html` are release-asset only (DESIGN.md Sec 8, 2026-09-24;
   `clinker/` added 2026-09-27, NovInvenio#202's clinker synteny panel --
   grows with drawn-locus count like the other release-asset-only entries).
   `clinker/` itself only stages the top 25 loci by default
   (`--clinker_publish_top`, N4, NovInvenio#202); the rest stay in the
   pipeline's own `results/<run>/.../pangenome/clinker/` dir (never copied
   into `docs/`), and `island_synteny.html` gets a small injected script so
   an unpublished locus's panel points there instead of a broken iframe.

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
`bin/build_study_config.py`, which dispatches each species' protein/genome/GFF3
independently to a UniProt fetch, an NCBI Datasets fetch, or a local-file copy per
its `Protein_Source`/`Genome_Source`/`GFF3_Source` columns — see
`notes/superpowers/specs/2026-09-11-study-onboarding-design.md`) if not
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
| Enrichment/analysis scripts (GO/Pfam/InterPro ORA) | `bin/go_enrichment.py` (GO, via goatools) + `bin/domain_enrichment.py` (Pfam/InterPro, custom hypergeometric) — both built, hypergeometric + BH-FDR, per the design in `DESIGN.md` Sec 7. Standalone CLI tools (candidates/background accession lists + `extract_dat_annotations.py` output as input) — not yet wired into `bin/run_study.sh`/`bin/sync_reports.sh` for automatic per-study invocation, and no dedicated `tests/test_go_enrichment.py`/`tests/test_domain_enrichment.py` yet |
| Site generation | `docs/` — release-asset published, never committed; see `DESIGN.md` Sec 8 |
| Pfam annotation for proteomes without a UniProt record (NCBI/JGI sources) | `bin/annotate_proteomes_pfam.py` (SLURM hmmsearch jobs + finalize) and `bin/hmm_domtbl_to_annotations.py` (HMMER domtblout -> `annotations/<Protein_Accession>.tsv.gz`); see "Annotating non-UniProt proteomes" below |

### Study-specific vs. shared scripts

`bin/` is for scripts any study could call, parameterized by CLI args (`fetch_uniprot_proteome.py --proteome-id ...`, `build_study_config.py --study-dir ...`). A script belongs in `studies/<domain>/<set_name>/bin/` instead the moment it stops being reusable in that sense — the giveaway is usually right there in the name or the body:

- The filename itself names a study, species, or one-off dataset (`build_akkermansia_config.py`, `fetch_kpn78578_modelorg.py`) rather than a general source/action (`fetch_uniprot_proteome.py`).
- It hardcodes a species list, an accession list, or a specific study's directory layout instead of taking them as arguments.
- It exists to solve one study's particular data-shape problem (e.g. proteomes that are already local files instead of a UniProt/NCBI pull-per-species) rather than a recipe any study would reuse.

Study-specific scripts still follow the shared-`lib/` and provenance conventions (import `lib/provenance.py`, write `.provenance.yaml` sidecars / fold into that study's `DATA_MANIFEST.yaml`) — only their own location, and any hardcoded reference back to `NII_ROOT`/`NII_BIN`, changes. Since they don't live next to the scripts they invoke via `subprocess` (e.g. `fetch_uniprot_proteome.py`) or the `lib/` they import, hardcode absolute paths to those (`NII_ROOT = Path("/bigdata/.../NovInvenio_Investigations")`) rather than deriving them from `__file__` or cwd — same reasoning as `~/.claude/CLAUDE.md`'s BASH_SOURCE/SLURM warning, just triggered by directory distance instead of a SLURM work dir.

### Annotating non-UniProt proteomes

Why: `bin/build_study_config.py` writes `annotations/<Protein_Accession>.tsv.gz` only for
`Protein_Source == uniprot` (from the UniProt `.dat.gz`). NCBI/JGI proteomes get no file,
so the enrichment analysis has no terms for them.

What: a Pfam-A hmmsearch (HMMER 3.4, `--cut_ga`) per proteome, converted into the same
file format. hmmsearch, not hmmscan: on 500 PiroE2sp proteins (4 cores, job 29943954)
hmmsearch took 81 s, hmmscan 362 s, with the same 525 (protein, Pfam) pairs. Measured
full-proteome time with `--cpu 8` on milan nodes: about 0.1 s per protein (20,214
proteins 1,911 s; 23,663 proteins 2,448 s; 27,675 proteins 2,765 s).

Commands:
```bash
bin/annotate_proteomes_pfam.py --study fungi/<study> --dry-run    # which species, how many jobs
bin/annotate_proteomes_pfam.py --study fungi/<study> --submit     # chunk jobs + dependent finalize job
bin/annotate_proteomes_pfam.py --study fungi/<study> --finalize   # by hand, if the finalize job did not run
# existing domain table (hmmsearch or hmmscan layout, plain/.gz/.zst):
bin/hmm_domtbl_to_annotations.py --domtbl X.domtbl.zst --fasta data_dir/pep/<Short>.pep.fa \
    --taxon-id <Taxon_ID> --out annotations/<Protein_Accession>.tsv.gz
```
`bin/run_study.sh` runs the `--submit` step after the config build only when
`NII_ANNOTATE_PFAM=1` is set (default: off). The jobs run on their own; the
enrichment run must wait until the annotation files exist.

Rules and outputs:
- Selection: species with `Protein_Source != uniprot` and no annotation file. UniProt
  files are never written or overwritten. `--force` redoes existing non-UniProt files.
- Job size: at most 35,000 proteins per hmmsearch run (larger proteomes are split by sequence count and the chunk tables are merged). Runs are packed into batch jobs of about 75 min estimated time (`--batch-minutes`, 0.1 s per protein). One finalize job (`afterany`) converts all species. Jobs use `--constraint=milan` and
  `$SCRATCH`; absolute paths only.
- Output name: `annotations/<name>.tsv.gz`. `<name>` comes from `lib/annotation_names.py`: the Protein_Accession for UniProt/NCBI accessions; for a local FASTA path, its basename without the proteome and compression suffix; Short as the fallback. Every reader and writer uses this helper. The scripts refuse to write outside `annotations/`, the work dir or `annotations_orig/`.
- Output: `annotations/<name>.tsv.gz` (UniProt columns; Pfam ids/names
  filled; two extra columns with best domain i-Evalue and score) + a
  `.provenance.yaml` sidecar (domain table path + sha256, HMMER version and command,
  Pfam release/path, FASTA sha256, script command). Work files:
  `analysis/pfam_hmmsearch/`. `annotations/` is gitignored (regenerable).
- `--fill-empty-uniprot` (opt-in): fills `pfam_ids`/`pfam_names` of UniProt files that have no terms. The original goes to `annotations_orig/`. Rerunning `build_study_config.py` reverts the file.
- Limitation: GO and InterPro are empty for these species. Pfam and Pfam-clan tests
  cover them; GO and InterPro tests do not. `enrich_novelties_v2.py` records this in
  PROVENANCE.yaml (`vocabulary_coverage`), and `bin/generate_enrichment_page.py` adds a
  caveat line ("GO and InterPro cover only N of M species: ...") to the page.
- Report sync: `bin/sync_reports.sh` merges these files into the `uniprot_*` columns
  together with the UniProt files. `bin/merge_uniprot_annotations.py` labels each row
  in a last column `annotation_source` (`uniprot` or `hmmer-pfam`; detection in
  `lib/annotation_sources.py`). `bin/annotation_source_note.py` then writes a note under
  the header of novelties/core/losses.html that names the species with `hmmer-pfam` rows
  in that page's matrix ("Pfam from a HMMER search (no GO/InterPro) for: ...").

## Mycelium

This repo adopts the Mycelium living-repo framework (`.living/`, `todo/`,
`mycelium:ingest`), same as `nf_NovInvenio`. Confirmed lightweight: the plugin itself is
a shared, one-time install; per-project footprint is a handful of tracked markdown/YAML
files plus a gitignored `.mycelium/` run-state dir.
