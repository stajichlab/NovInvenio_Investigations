# `ni pangenome`: one launcher for pangenome.nf runs — design (2026-09-26)

Status: **design approved in chat 2026-09-26 (sections 1-3). Nothing is implemented.**

## Goal

A study runs pangenome.nf through `bin/ni`, by run name, from one committed
runs file. The same command works for a new dataset, for a two-species
comparison (either direction), and for a single-species run. It replaces the
six hand-copied `studies/fungi/*/run_pangenome.sh` scripts and their
hand-written `submit_nextflow_head.sh` files.

## Facts this design rests on (checked 2026-09-25/26)

- `bin/ni run` calls `bin/run_study.sh`, which runs `main.nf` only
  (`nextflow run "$PIPELINE"`, line 67). No `ni` command runs pangenome.nf.
- Six studies have their own `run_pangenome.sh`: Afumigatus_pangenome,
  coccidioides_pangenome, FSSC_ambrosia, FSSC_solani, fusarium_FOL,
  fusarium_FOXY_vs_FSSC. They differ in how they pin the pipeline (live
  checkout, detached worktree, `git archive` snapshot), in their variant
  names, in how they set rescue (env var or params YAML), and in whether they
  pass Pfam. None passes a species tree.
- Nextflow 26.04.6 keeps one clone per commit:
  `nextflow pull stajichlab/NovInvenio -r <full sha>` created
  `~/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/<sha>/`, which contains
  `pangenome.nf`. A short SHA fails: "Remote does not have 0dddb42 available
  for fetch."
- The pipeline's non-container `beforeScript` (`nextflow.config`, about line
  455) activates pixi from `${projectDir}/pixi.toml` and links
  `${projectDir}/db` and `${projectDir}/configs` into each task dir.
- The clone has a tracked `db/` directory. Its Pfam files are tracked
  symlinks to `/srv/projects/db/pfam/2026-01-27-Pfam38.2/`. pangenome.nf and
  its modules/workflows never read a `db/` or `projectDir` path. Every study
  passes Pfam as an absolute path.
- `bin/sync_pangenome_report.py` finds `results/<run>/*/pangenome/report/report.md`
  and takes `--study`, `--run`, `--current`, `--pangenome-dir`.
- Boolean params given on the CLI arrived as strings before #191 (fixed in
  NovInvenio `3ac3919`). The Fusarium studies pass params by YAML file.
- `bin/ni` needs Python 3.10+ (`str | None` annotations). UCR HPCC's default
  `python3` is 3.9. Since 2026-09-26 `bin/ni` has `#!/usr/bin/env python3.12`
  and exits with an error under < 3.10. The generated head job calls
  `/usr/bin/python3.12` by absolute path.

## Decisions

| Question | Decision |
|---|---|
| Where runs are defined | One tracked file per study: `studies/<domain>/<set>/pangenome_runs.yaml` (class 2). |
| Which runs are listed | New runs only. Existing results folders are not back-filled. |
| Pipeline pinning | Remote commit through Nextflow: `nextflow run stajichlab/NovInvenio -r <full sha>`. |
| Tool environment | `ni` runs `pixi install` once in the commit's clone. No `db/` setup. |
| Launch | `ni` writes and submits a SLURM head job. On success the head job stages docs locally. Release upload is a separate command. |
| Reuse of a results folder | Refused when it was made at a different commit. |
| Subcommands | `list`, `check`, `run`, `stage`, `publish`. |

## 1. Runs file

`studies/<domain>/<set>/pangenome_runs.yaml`:

```yaml
defaults:
  pipeline: stajichlab/NovInvenio          # GitHub repo
  pipeline_commit: 0dddb421720e0eb0a3045231608b27b113d49be6   # full SHA, required
  pfam_hmm: /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm
  queue_config: stajichlab_queue.config    # study-relative
  data_dir: data_dir                       # study-relative, default data_dir
  head_job: {partition: stajichlab, account: null, time: 30-00:00:00, mem: 8G, cpus: 1}
  params: {pangenome_cluster_backend: mmseqs}
runs:
  immitis_in_posadasii_out_0dddb42:
    samplesheet: config_immitis_in_posadasii_out.csv
    species_tree: analysis/species_tree/coccidioides_coccidioides_only.full.nwk
    params: {pangenome_rescue_enable: true}
    publish: true
    current: true
  immitis_only_0dddb42:
    samplesheet: config_immitis.csv
    params: {pangenome_rescue_enable: false}
```

Rules:

- A run inherits every `defaults` key and can override any of them. `params`
  merges one key at a time (run value wins). `head_job` merges the same way.
- Allowed run keys: `samplesheet` (required), `species_tree`, `params`,
  `publish` (default false), `current` (default false), plus any `defaults` key.
  Any other key is an error.
- The run name must match `SLUG_RE` (`^[A-Za-z0-9._-]+$`, `lib/study_runs.py`).
  It is the folder name in `results/<run>/`, `.nf_launch/<run>/` and
  `docs/<domain>/<study>/<run>/`.
- Relative paths resolve against the study directory. `pfam_hmm` must be
  absolute.
- `publish: true` means the run is meant for the site. The head job stages it
  after success, and `ni pangenome publish` accepts it. A run with
  `publish: true` needs the study's `publish.yaml`.
- `current: true` passes `--current` to `sync_pangenome_report.py`. At most
  one run per file may set it.

## 2. CLI

New subcommand group in `bin/ni`. The logic goes in a new `lib/ni_pangenome.py`.
Every command takes `--study-dir studies/<domain>/<set>`.

| Command | Does |
|---|---|
| `ni pangenome list` | One line per run: run name, short commit, samplesheet row count, state (`not started`, `submitted <jobid>`, `failed`, `done`, `staged`). State comes from `ni_run.json` and `docs/.../<run>/run.json`. |
| `ni pangenome check --run R` | Loads the file and runs every check in section 3. Prints each failure. Exit 1 on any failure. Submits nothing. |
| `ni pangenome run --run R [--dry-run] [--foreground] [-- extra nextflow args]` | `check`, pull, install, write `params.yaml` and the head job, then `sbatch` it (section 4). |
| `ni pangenome stage --run R` | Calls `bin/sync_pangenome_report.py --study <d>/<s> --run R [--current]`. |
| `ni pangenome publish --run R` | Requires `publish: true` and a staged run (`docs/<d>/<study>/<run>/run.json` exists). Calls `bin/publish_report_release.sh <d>/<study>/<run>`. Never called by any other command. |

`<study>` is the target study name from `publish.yaml`. For studies without
`publish.yaml`, `stage` and `publish` refuse with a message.

## 3. Checks (`check`, and the first step of `run`)

1. The YAML parses. It has only known keys. Run names match `SLUG_RE`.
2. `pipeline_commit` is 40 hex characters.
3. The samplesheet exists and has a `Short` column. `Short` values are unique.
4. Every protein, genome and GFF3 file the samplesheet names exists under
   `data_dir/` (`pep/`, `dna/`, `gff3/`).
5. If `species_tree` is set: the file exists, and its tip set equals the
   samplesheet `Short` set. The error lists the extra and the missing names.
6. `pfam_hmm` exists.
7. `queue_config` exists.
8. No input path resolves inside `~/.nextflow/assets/`.
9. If `.nf_launch/<run>/ni_run.json` exists and records a different commit,
   or `results/<run>/` exists without an `ni_run.json`: refuse. The message
   says to choose a new run name.

## 4. Pinning, setup and launch

Steps of `ni pangenome run`, in order. Steps 1-3 run on the submitting node,
so errors appear before any job is queued.

1. **Pull.** `nextflow pull <pipeline> -r <sha>`. The clone path is
   `~/.nextflow/assets/.repos/<org>/<repo>/clones/<sha>/`. If the commit is not
   on GitHub, the pull fails and `ni` stops.
2. **Install.** If `<clone>/.ni_provisioned` does not hold the sha256 of
   `<clone>/pixi.lock`: run `pixi install --frozen --manifest-path <clone>/pixi.toml`,
   then write the marker. `ni` prints the install time and the size of
   `<clone>/.pixi`. (Env size is not yet measured. The run worktrees are about
   10 GB each. It is an assumption that most of that is the pixi env.)
3. **Write `.nf_launch/<run>/params.yaml`.** Content: the merged `params`, plus
   `pangenome_samplesheet`, `pangenome_data_dir`, `outdir: <study>/results/<run>`,
   and `pangenome_island_pfam_hmm` / `pangenome_species_tree` when set. All paths
   absolute. YAML booleans stay booleans. `--pangenome_project` is not set, so
   outputs land in `results/<run>/output/pangenome/` (NovInvenio #194: publishDir
   ignores the project name pangenome.nf derives). Verified by the FOL smoke run
   at 0dddb42 (2026-09-26). A 2026-09-26 correction to "<samplesheet stem>/" was
   wrong and is reverted. `sync_pangenome_report.py` finds `*/pangenome/` either way.
4. **Write `.nf_launch/<run>/submit_nextflow_head.sh`** from a template in
   `lib/ni_pangenome.py`:
   - `#SBATCH` lines from `head_job`: `-J nf-<study>-<run>`, `-p`, `-A` (only
     if set), `-c`, `--mem`, `-t`, and `-o`/`-e` to `.nf_launch/<run>/slurm-%j.{out,err}`.
   - Absolute paths only. No `BASH_SOURCE`.
   - Body: `cd` to `.nf_launch/<run>/`, then
     ```
     nextflow run <pipeline> -r <sha> -main-script pangenome.nf \
       -params-file <launch>/params.yaml \
       -profile slurm -c <clone>/conf/ucr_hpcc_slurm.config -c <study>/<queue_config> \
       -with-trace <study>/results/<run>/trace.txt -resume <extra args>
     ```
     then record the exit status in `ni_run.json`. If the exit status is 0 and
     `publish: true`, run
     `/usr/bin/python3.12 <NII>/bin/ni pangenome stage --study-dir <abs> --run <run>`.
   - The script is regenerated on every `ni pangenome run`. To resume, rerun
     `ni pangenome run` with the same run name (same commit, so check 9 passes).
5. **Write `.nf_launch/<run>/ni_run.json`**: `run`, `pipeline`, `pipeline_commit`,
   `clone_path`, `samplesheet` and its sha256, sha256 of `params.yaml`,
   `submitted_at` (UTC), `slurm_job_id`, `exit_status` (null until the job ends).
6. **Submit.** `sbatch <launch>/submit_nextflow_head.sh`. The job ID goes into
   `ni_run.json`.

Flags:

- `--dry-run`: run the checks, then print `params.yaml` and the head-job script.
  No pull, no install, no files written, no `sbatch`.
- `--foreground`: steps 1-5, then run the script body in the current shell
  instead of `sbatch`.

## 5. Stage and publish

- `stage` runs `bin/sync_pangenome_report.py`. It passes `--current` when the
  run sets `current: true`.
- Change to `bin/sync_pangenome_report.py`: if `.nf_launch/<run>/ni_run.json`
  exists, copy its `pipeline` and `pipeline_commit` into the `run.json` it
  writes. Runs made before `ni` are unchanged.
- `publish` calls `bin/publish_report_release.sh <domain>/<study>/<run>`. It is
  an explicit, outward-facing step (GitHub Release plus Pages deploy).

## 6. Adoption

- Runs files list new runs only.
- A study keeps its old `run_pangenome.sh` and submit scripts until its first
  `ni` run succeeds. They are then removed in a commit for that study.
- First adopters:
  - `fusarium_FOL` norescue (25 strains): the live smoke test.
  - `coccidioides_pangenome`: new reciprocal runs (both directions), plus
    immitis-only and posadasii-only runs, with the species tree. These run at a
    NovInvenio commit that contains the island locus view.

## 7. Tests and verification

Unit tests in `tests/test_ni_pangenome.py` (pytest, pattern of
`tests/test_ni_run.py`). Each test builds a fake study in `tmp_path`.

1. Runs-file loading: inheritance, one-key `params` and `head_job` merge,
   unknown keys rejected, bad slug rejected, two `current: true` rejected.
2. One test per check in section 3, each asserting the error text.
3. `params.yaml`: booleans stay booleans; paths are absolute.
4. Head-job script: compared with a checked-in expected file. No `BASH_SOURCE`.
   `stage` runs only on exit 0 and only for `publish: true`.
5. `run` with patched subprocess (fake `nextflow`, `pixi`, `sbatch`): order is
   pull, install, write, submit; install skipped when the marker matches;
   `--dry-run` calls nothing external and writes nothing; `--foreground` does
   not call `sbatch`.
6. `stage`/`publish`: the right arguments reach the two scripts; `publish`
   refuses an unstaged run and a run without `publish: true`.
7. `tests/test_sync_pangenome_report.py`: `pipeline_commit` from `ni_run.json`
   reaches `run.json`.

Live verification, reported with real output:

1. `ni pangenome check` and `run --dry-run` on the Cocci runs file.
2. The real pull and install for one commit. Record install time and env size.
3. `fusarium_FOL` norescue through `ni pangenome run`. It passes when the run
   completes, `stage` writes `docs/fungi/fusarium_FOL/<run>/`, and `run.json`
   holds the commit. Compare family count and matrix size with the existing
   `results/mmseqs_norescue` (pipeline `af6fd68`). Differences are reported,
   not treated as failure, because the commits differ. No publish unless asked.

## Out of scope

- Generating reciprocal or per-species samplesheets
  (`build_reciprocal_configs.py`, `filter_config_by_taxon.py` stay as they are).
- Container (`-profile singularity`) runs.
- Changes to `main.nf` studies, `bin/run_study.sh`, or pipeline code.
- Automatic cleanup of old clones in `~/.nextflow/assets`.

## Open items

- pixi env size and install time per commit: measured in live step 2.
- NovInvenio #193: `BUILD_ISLANDS` runs only when `--pangenome_island_pfam_hmm`
  is set. Until it is fixed, a run without `pfam_hmm` has no islands. `check`
  prints a warning when `pfam_hmm` is unset.
- NovInvenio #189: `BUILD_PRESENCE_MATRIX` fixed at 4 GB. Studies handle it in
  their `queue_config` until it is fixed.
