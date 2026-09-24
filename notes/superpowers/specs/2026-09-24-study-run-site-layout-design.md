# Study/run site layout — design (2026-09-24)

Status: **approved 2026-09-24**; code and local migration implemented on branch
`study-run-layout` (see "Implementation notes" at the end). Release re-tagging
and old-tag deletion (Migration steps 3-5) not done. Nothing below is implemented yet except the
pangenome.nf half of step 1 (`bin/sync_pangenome_report.py`, merged `a4a5222`).

## Goal

Every published study uses one layout: `docs/<domain>/<study>/<run>/`. A run is
one pipeline result for the study's species set. Old runs stay published as
archived versions, so results from different parameters, dates, or pipeline
versions can be compared.

## Decisions (made 2026-09-24)

| Question | Decision |
|---|---|
| Layout | One layout for all studies: `docs/<domain>/<study>/<run>/` |
| Scope | **Site only.** `studies/` folders do not move. Each folder declares which study/run it publishes as. |
| Run names | Free slug, `[A-Za-z0-9._-]+`. Date, parameters, pipeline commit are recorded in `run.json`, not required in the name. |
| Current run | Each study marks one run as current. The gallery card links to it; other runs are listed as archived. |

## Facts this design rests on (checked 2026-09-24)

- Variant folders of one study have identical `species.csv` (same md5) and
  symlinked `config.csv`/`data_dir`. They differ only in `run_params.txt`.
  Exception: `agaricomycetes_novelty_discovery` has a different `species.csv`
  (9 lines vs 8).
- `main.nf` pages link to each other by sibling-relative paths
  (`report.html` → `novelties.html`/`core.html`/`losses.html`;
  `alignment.html` fetches `alignments/...`). A set can move down one level
  without breaking its own links.
- 11 flat sets are tracked in git (3 files each: `report.html`,
  `alignment.html`, `index.html`). 25 `alignments-*`/`reports-*` releases exist.

## Layout

```
docs/<domain>/<study>/
  report.html        run list; current run highlighted          (committed)
  index.html         redirect -> current run's report.html      (committed)
  study.json         {"current": "<run>"}                       (committed)
  <run>/
    report.html      run landing page                            (committed)
    run.json         run metadata, see below                    (committed)
    alignment.html   main.nf viewer shell                        (committed)
    index.html       redirect -> report.html                     (committed)
    novelties.html core.html losses.html summary.pdf alignments/ loss_alignments/
    figures/ figures_pdf/ archive/ island_synteny.html assembly_quality.html
                     release asset only (gitignored)
```

`run.json` fields: `run`, `published` (UTC date), `pipeline` (`main.nf` or
`pangenome.nf`), `source_dir`, `params` (non-comment lines of the folder's
`run_params.txt`, main.nf only), `report_sha256`, counts (`n_families`,
`n_strains` for pangenome; `n_species` for main.nf), `pipeline_commit` when it is
known (optional — not recorded anywhere reliable today for main.nf runs).

## Declaring study and run: `studies/<domain>/<folder>/publish.yaml`

```yaml
study: pezizo_set1      # docs/<domain>/<study>/
run: mmseqs-cov0.3      # docs/<domain>/<study>/<run>/  (main.nf studies)
```

- **Absent file: the folder is not published.** The sync step prints a
  message saying `publish.yaml` is missing and skips the docs/ step (exit 0,
  so `bin/run_study.sh`'s post-run chain does not fail). Publishing is opt-in
  per folder.
- pangenome.nf folders also need `publish.yaml`, with `study` only; the run
  comes from `--run <results/ subdir>`, because one folder holds several runs.
- CLI flags (`--study`, `--run`) override the file.

## Code changes

1. `bin/sync_reports.sh`: read `publish.yaml`; write to
   `docs/<domain>/<study>/<run>/` instead of `docs/<domain>/<folder>/`; write
   `run.json`; rebuild the study run list. `results/<folder>/` and the
   pipeline's flat `docs/<folder>/` are unchanged inputs.
2. Shared study-page code: move the run-list/`run.json` handling from
   `bin/sync_pangenome_report.py` into `lib/` so both sync scripts use it. Add
   `study.json` + `--current` (and `bin/set_current_run.py <domain>/<study> <run>`).
   A study with one run makes it current automatically.
3. `bin/publish_alignment_release.sh`: accept `<domain>/<study>/<run>` (tag
   `alignments-<domain>-<study>--<run>`), same as the report script already does.
4. `bin/generate_docs.py` + `lib/site_pages.py`: study cards come from
   `docs/<domain>/<study>/study.json` + `*/run.json`, plus pending/staged
   folders from `studies/` grouped by `publish.yaml`. Card shows the run count
   and links to the current run.
5. `.gitignore`: move the `novelties/core/losses.html`, `summary.pdf`,
   `alignments/`, `loss_alignments/` rules down to `docs/*/*/*/`. The set level
   holds only small generated files (run list, redirects, `study.json`).
6. `bin/publish_all_studies.sh`: iterate runs, not sets.

`static.yml` already places a release by its manifest's `run` field (merged in
`a4a5222`); it needs no further change.

## Old URLs

GitHub Pages has no server-side redirects. For each migrated flat set, commit
small meta-refresh stubs at the old paths (`index.html`, `report.html`,
`novelties.html`, `core.html`, `losses.html`, `alignment.html`) pointing to the
same file in the new run folder. Where the old set name becomes the new study
name (e.g. `pezizo_set1`), `report.html` becomes the run list and the other stubs
point to the current run.

Conflict to resolve in implementation: stubs named `novelties.html` etc. at set
level would match today's `.gitignore` rules. Step 5 moves those rules to run
level, so the stubs can be committed. A test must check that every set-level
`novelties.html`/`core.html`/`losses.html` is a redirect stub (e.g. < 2 KB and
contains `http-equiv="refresh"`), so a real report can never be committed there
by mistake.

Old `alignments/` shard URLs are not redirected (the viewer fetches them
relative to its own page, which moves with them).

## Migration

Mapping for the 11 tracked flat sets (edited and approved 2026-09-24). The
default run name is `pairwise-diamond`. `sordariales_shallow` has no
`run_params.txt` lines; it is named `pairwise-diamond` per this decision.

| Old folder | Study | Run | Current? | Params (from `run_params.txt`) |
|---|---|---|---|---|
| pezizo_set1 | pezizo_set1 | pairwise-diamond | yes | `--run_tool diamond` |
| pezizo_set1_cluster | pezizo_set1 | mmseqs-cov0.3 | | `--run_tool diamond --cluster_tool mmseqs --hmm_presence_cov 0.3 --hmm_presence_min_residues 100` |
| sordariales_shallow | sordariales_shallow | pairwise-diamond | yes | (none) |
| sordariales_shallow_cluster | sordariales_shallow | mmseqs-cov0.3 | | same as pezizo_set1_cluster |
| agaricomycetes_pairwise | agaricomycetes | pairwise-diamond | yes | `--cluster_tool pairwise --run_tool diamond` |
| agaricomycetes_mmseqs | agaricomycetes | mmseqs-cov0.3 | | `... --cluster_tool mmseqs --hmm_presence_cov 0.3 --hmm_presence_min_residues 100 --hmm_presence_domain_evalue 0.01` |
| agaricomycetes_novelty_discovery | agaricomycetes | novelty-discovery | | `--run_tool diamond --cluster_tool novelty_discovery --hmm_presence_domain_evalue 0.01` — different species.csv (9 vs 8 lines). Kept as a run of `agaricomycetes`; probably an old run, likely removed later. |
| mushrooms_tremella | mushrooms_tremella | pairwise-diamond | yes | `--run_tool diamond` |
| yeast_filamentous | yeast_filamentous | pairwise-diamond | yes | `--run_tool diamond` |
| zoosporic_dikarya | zoosporic_dikarya | pairwise-diamond | yes | `--run_tool diamond` |
| UHM_Akkermansia / UHM_Koxytoca / UHM_lachnoNovelclade | same | pairwise-diamond | yes | `--run_tool diamond --cluster_tool pairwise` + Pfam/SwissProt/modelorgs |

Not migrated: `*_dmnd_{default,sensitive,very_sensitive}` (temporary benchmark
folders, to be removed later; never committed, no releases). Pangenome studies
have nothing published yet; they start in the new layout.

Migration order (each step is visible; steps 3 and 5 change the live site):

1. Write `publish.yaml` into each old folder per the table (committed).
2. `git mv` each flat set's committed files into `<study>/<run>/`; write
   `run.json`, `study.json`, run lists, and the old-URL stubs. Move the
   gitignored local files the same way (plain `mv`).
3. For each old `reports-<d>-<folder>` / `alignments-<d>-<folder>` release:
   download it, rewrite `manifest.json` with the new study and run, and upload
   it as `reports-<d>-<study>--<run>` / `alignments-<d>-<study>--<run>`. No
   pipeline rerun needed.
4. Commit, deploy, and check every old URL and every new run page on the live
   site.
5. **Only after step 4 passes, and with explicit approval:** delete the 25 old
   release tags, then redeploy. Until then the old tags keep merging old
   content into the old set folders, over the stubs.

## Answered questions (2026-09-24)

1. Run names and current runs: as in the migration table. Default run name is
   `pairwise-diamond`.
2. `agaricomycetes_novelty_discovery`: a run of `agaricomycetes` (probably an
   old run; likely removed later).
3. Folder without `publish.yaml`: not published. The sync step asks for one
   and skips.

## Implementation notes (2026-09-24)

Changes from the design above, found while implementing:

- `set-current` is `bin/study_pages.py set-current <domain>/<study> <run>`, not a
  separate `bin/set_current_run.py`. The same CLI has `target`, `record-run`,
  and `rebuild`.
- **Set-level `.gitignore` rules are kept**, not moved. Run-level copies were
  added. Reason: old flat folders and the temporary `*_dmnd_*` folders still
  hold real `novelties.html` etc. at set level; un-ignoring those names would
  let an `git add -A` commit them. The redirect stubs are force-added
  (`git add -f`) by `bin/migrate_docs_to_runs.py`, and
  `tests/test_study_runs.py` checks every tracked set-level
  `novelties/core/losses.html` is a stub.
- **Old-URL stubs point at the migrated run**, not the current run. An old URL
  meant that specific result, and the stub then never needs rewriting when the
  current run changes.
- A complete gallery card is a `<div>` with two links (title -> current run,
  "N runs" -> run list), not an `<a>`: links cannot nest, and main.nf report
  pages have no link back to the study, so the run list needs its own link.
- **Folders without `publish.yaml` get no gallery card at all** (before, every
  `studies/` folder with `species.csv` got a pending/staged card). Cards that
  disappear with this change: `fusarium_FOXY`, `cyanobacteria` (and the
  variant/test folders, which are now runs or unpublished). Add a
  `publish.yaml` to bring a card back.
- `run.json` for a migrated run: `published` = date of the last commit that
  touched the old `report.html`; `note` says it was migrated and that `params`
  come from `run_params.txt` at migration time.
- The migration is re-runnable. In a checkout that received the committed
  migration by merge, a second `--apply` moves the leftover gitignored files
  (`novelties.html`, `alignments/`, ...) into the run folders and writes and
  force-adds their old-URL stubs.
