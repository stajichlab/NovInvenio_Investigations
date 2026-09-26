# `ni pangenome` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `bin/ni pangenome {list,check,run,stage,publish}` so any NII study runs pangenome.nf by run name from a committed `pangenome_runs.yaml`.

**Architecture:** All logic lives in a new `lib/ni_pangenome.py` (pure functions plus one orchestration function that takes an injectable `runner` for subprocess calls). `bin/ni` only parses arguments and calls it. `bin/sync_pangenome_report.py` gets one small change: copy the pipeline commit from the run record into `run.json`.

**Tech Stack:** Python 3.12 (`/usr/bin/python3.12`), PyYAML, pytest (run via `pixi run python -m pytest` in the NII repo root), Nextflow 26.04.6, pixi, SLURM `sbatch`.

**Spec:** `notes/superpowers/specs/2026-09-26-ni-pangenome-mode-design.md`

## Global Constraints

- Repo root: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations`. Run all commands from there.
- Tests: `pixi run python -m pytest -q <path>`. The system `python3` is 3.9 and lacks scipy; never use it.
- `bin/ni` shebang is `#!/usr/bin/env python3.12` with a `< (3, 10)` guard. Keep both.
- No `BASH_SOURCE` in any generated script. All paths in generated files are absolute.
- Runs file name: `studies/<domain>/<set>/pangenome_runs.yaml`. Top-level keys: `defaults`, `runs`.
- `pipeline_commit` must match `^[0-9a-f]{40}$`.
- Run names must match `^[A-Za-z0-9._-]+$` (reuse `study_runs.validate_slug`).
- Default `pipeline`: `stajichlab/NovInvenio`. Default `data_dir`: `data_dir`. Default `head_job`: `{partition: stajichlab, account: null, time: 30-00:00:00, mem: 8G, cpus: 1}`.
- Clone path: `~/.nextflow/assets/.repos/<org>/<repo>/clones/<sha>/`.
- Outputs: `results/<run>/output/pangenome/` (no `--pangenome_project`). Launch dir: `studies/<d>/<s>/.nf_launch/<run>/`.
- `publish` is never called by any other command. Staging never publishes.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit to `main` (repo practice); do not push (the controller pushes).

## Review Focus

1. A samplesheet with blank `DNA`/`GFF3` cells (proteome-only strains): those cells must be skipped by the file-exists check, not reported as missing `data_dir/dna/`. Test in Task 2.
2. A Newick tree with quoted labels, branch lengths, support values and internal node labels: only leaf labels count as tips. Test in Task 2.
3. Resuming the same run at the same commit after a failed head job: `check` must pass (ni_run.json exists, same commit) and `run` must not reinstall pixi. Tests in Task 2 and Task 4.
4. A study with no `publish.yaml` and a run with `publish: true`: `check` must fail with a message naming `publish.yaml`. Test in Task 2.
5. `sbatch` fails (non-zero, or output not a job ID): `run` must return non-zero and leave `slurm_job_id` null in `ni_run.json`. Test in Task 4.

---

## File Structure

| File | Responsibility |
|---|---|
| Create `lib/ni_pangenome.py` | Runs-file loading (`RunSpec`), checks, rendering (params YAML, head-job script, run record), orchestration (`run_pipeline`), stage/publish command builders, `record_exit` |
| Create `tests/test_ni_pangenome.py` | Unit tests for all of the above |
| Create `tests/data/ni_pangenome/expected_head_job.sh` | Expected head-job script (golden file) |
| Modify `bin/ni` | Add the `pangenome` subcommand group; dispatch to `lib/ni_pangenome.py` |
| Modify `bin/sync_pangenome_report.py` | `stage_run` copies `pipeline`/`pipeline_commit` from `.nf_launch/<run>/ni_run.json` |
| Modify `tests/test_sync_pangenome_report.py` | Test for the above |
| Create `studies/fungi/fusarium_FOL/pangenome_runs.yaml` | Smoke-test runs file (Task 7) |
| Create `studies/fungi/coccidioides_pangenome/pangenome_runs.yaml` | Cocci runs file, checked with `--dry-run` only (Task 7) |

---

### Task 1: Runs-file loading

**Files:**
- Create: `lib/ni_pangenome.py`
- Test: `tests/test_ni_pangenome.py`

**Interfaces:**
- Consumes: `study_runs.validate_slug(name: str, what: str) -> str` (raises `ValueError`).
- Produces:
  - `class RunsFileError(ValueError)`
  - `@dataclass(frozen=True) class RunSpec` with fields `name: str, study_dir: Path, pipeline: str, pipeline_commit: str, samplesheet: Path, data_dir: Path, pfam_hmm: Path | None, queue_config: Path | None, species_tree: Path | None, head_job: dict, params: dict, publish: bool, current: bool`
  - `RUNS_FILE = "pangenome_runs.yaml"`
  - `load_runs_file(study_dir: Path) -> dict[str, RunSpec]` (keeps file order)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_ni_pangenome.py`:

```python
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))

import ni_pangenome as nip  # noqa: E402

SHA = "0dddb421720e0eb0a3045231608b27b113d49be6"


def write_runs(study: Path, doc: dict) -> None:
    study.mkdir(parents=True, exist_ok=True)
    (study / nip.RUNS_FILE).write_text(yaml.safe_dump(doc, sort_keys=False))


def base_doc(**run_overrides) -> dict:
    run = {"samplesheet": "config.csv"}
    run.update(run_overrides)
    return {
        "defaults": {
            "pipeline_commit": SHA,
            "pfam_hmm": "/abs/Pfam-A.hmm",
            "queue_config": "queue.config",
            "head_job": {"partition": "stajichlab", "time": "3-00:00:00"},
            "params": {"pangenome_cluster_backend": "mmseqs", "pangenome_rescue_enable": True},
        },
        "runs": {"r1": run},
    }


def test_run_inherits_defaults_and_resolves_paths(tmp_path):
    study = tmp_path / "studies" / "fungi" / "s1"
    write_runs(study, base_doc())
    spec = nip.load_runs_file(study)["r1"]
    assert spec.name == "r1"
    assert spec.pipeline == "stajichlab/NovInvenio"
    assert spec.pipeline_commit == SHA
    assert spec.samplesheet == study / "config.csv"
    assert spec.data_dir == study / "data_dir"
    assert spec.queue_config == study / "queue.config"
    assert spec.pfam_hmm == Path("/abs/Pfam-A.hmm")
    assert spec.species_tree is None
    assert spec.publish is False and spec.current is False


def test_params_and_head_job_merge_one_key_at_a_time(tmp_path):
    study = tmp_path / "s1"
    write_runs(study, base_doc(params={"pangenome_rescue_enable": False},
                               head_job={"mem": "16G"}))
    spec = nip.load_runs_file(study)["r1"]
    assert spec.params == {"pangenome_cluster_backend": "mmseqs", "pangenome_rescue_enable": False}
    assert spec.head_job == {"partition": "stajichlab", "account": None,
                             "time": "3-00:00:00", "mem": "16G", "cpus": 1}


def test_unknown_run_key_rejected(tmp_path):
    study = tmp_path / "s1"
    write_runs(study, base_doc(samplesheat="typo.csv"))
    with pytest.raises(nip.RunsFileError, match="samplesheat"):
        nip.load_runs_file(study)


def test_unknown_top_level_key_rejected(tmp_path):
    study = tmp_path / "s1"
    doc = base_doc()
    doc["extra"] = 1
    write_runs(study, doc)
    with pytest.raises(nip.RunsFileError, match="extra"):
        nip.load_runs_file(study)


def test_bad_run_name_rejected(tmp_path):
    study = tmp_path / "s1"
    doc = base_doc()
    doc["runs"] = {"bad name": {"samplesheet": "c.csv"}}
    write_runs(study, doc)
    with pytest.raises(nip.RunsFileError, match="bad name"):
        nip.load_runs_file(study)


def test_two_current_runs_rejected(tmp_path):
    study = tmp_path / "s1"
    doc = base_doc(current=True)
    doc["runs"]["r2"] = {"samplesheet": "c.csv", "current": True}
    write_runs(study, doc)
    with pytest.raises(nip.RunsFileError, match="current"):
        nip.load_runs_file(study)


def test_missing_samplesheet_key_rejected(tmp_path):
    study = tmp_path / "s1"
    doc = base_doc()
    doc["runs"]["r1"] = {"publish": True}
    write_runs(study, doc)
    with pytest.raises(nip.RunsFileError, match="samplesheet"):
        nip.load_runs_file(study)


def test_relative_pfam_rejected(tmp_path):
    study = tmp_path / "s1"
    doc = base_doc()
    doc["defaults"]["pfam_hmm"] = "db/Pfam-A.hmm"
    write_runs(study, doc)
    with pytest.raises(nip.RunsFileError, match="pfam_hmm"):
        nip.load_runs_file(study)


def test_missing_runs_file(tmp_path):
    with pytest.raises(nip.RunsFileError, match=nip.RUNS_FILE):
        nip.load_runs_file(tmp_path)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'ni_pangenome'`

- [ ] **Step 3: Write the implementation**

Create `lib/ni_pangenome.py`:

```python
"""ni pangenome: run nf_NovInvenio's pangenome.nf for an NII study by run name.

Spec: notes/superpowers/specs/2026-09-26-ni-pangenome-mode-design.md.
A study declares its runs in studies/<domain>/<set>/pangenome_runs.yaml.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from study_runs import validate_slug

RUNS_FILE = "pangenome_runs.yaml"
DEFAULT_PIPELINE = "stajichlab/NovInvenio"
DEFAULT_HEAD_JOB = {"partition": "stajichlab", "account": None, "time": "30-00:00:00",
                    "mem": "8G", "cpus": 1}
DEFAULT_KEYS = {"pipeline", "pipeline_commit", "pfam_hmm", "queue_config", "data_dir",
                "head_job", "params"}
RUN_KEYS = DEFAULT_KEYS | {"samplesheet", "species_tree", "publish", "current"}


class RunsFileError(ValueError):
    """The runs file is missing or malformed."""


@dataclass(frozen=True)
class RunSpec:
    name: str
    study_dir: Path
    pipeline: str
    pipeline_commit: str
    samplesheet: Path
    data_dir: Path
    pfam_hmm: Path | None
    queue_config: Path | None
    species_tree: Path | None
    head_job: dict
    params: dict
    publish: bool
    current: bool


def _study_path(study_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    p = Path(value)
    return p if p.is_absolute() else study_dir / p


def load_runs_file(study_dir: Path) -> dict[str, RunSpec]:
    study_dir = Path(study_dir)
    path = study_dir / RUNS_FILE
    if not path.is_file():
        raise RunsFileError(f"{path} not found (a pangenome study needs {RUNS_FILE})")
    doc = yaml.safe_load(path.read_text()) or {}
    extra = set(doc) - {"defaults", "runs"}
    if extra:
        raise RunsFileError(f"{path}: unknown top-level key(s): {sorted(extra)}")
    defaults = doc.get("defaults") or {}
    bad = set(defaults) - DEFAULT_KEYS
    if bad:
        raise RunsFileError(f"{path}: unknown defaults key(s): {sorted(bad)}")
    runs = doc.get("runs") or {}
    if not runs:
        raise RunsFileError(f"{path}: no runs defined")

    specs: dict[str, RunSpec] = {}
    for name, run in runs.items():
        try:
            validate_slug(str(name), "run name")
        except ValueError as e:
            raise RunsFileError(f"{path}: {e}") from None
        run = run or {}
        bad = set(run) - RUN_KEYS
        if bad:
            raise RunsFileError(f"{path}: run {name!r}: unknown key(s): {sorted(bad)}")
        merged = {**defaults, **run}
        merged["params"] = {**(defaults.get("params") or {}), **(run.get("params") or {})}
        merged["head_job"] = {**DEFAULT_HEAD_JOB, **(defaults.get("head_job") or {}),
                              **(run.get("head_job") or {})}
        if not merged.get("samplesheet"):
            raise RunsFileError(f"{path}: run {name!r}: samplesheet is required")
        if not merged.get("pipeline_commit"):
            raise RunsFileError(f"{path}: run {name!r}: pipeline_commit is required")
        pfam = merged.get("pfam_hmm")
        if pfam is not None and not Path(pfam).is_absolute():
            raise RunsFileError(f"{path}: run {name!r}: pfam_hmm must be an absolute path")
        specs[name] = RunSpec(
            name=name,
            study_dir=study_dir,
            pipeline=merged.get("pipeline", DEFAULT_PIPELINE),
            pipeline_commit=str(merged["pipeline_commit"]),
            samplesheet=_study_path(study_dir, merged["samplesheet"]),
            data_dir=_study_path(study_dir, merged.get("data_dir", "data_dir")),
            pfam_hmm=Path(pfam) if pfam else None,
            queue_config=_study_path(study_dir, merged.get("queue_config")),
            species_tree=_study_path(study_dir, merged.get("species_tree")),
            head_job=merged["head_job"],
            params=merged["params"],
            publish=bool(merged.get("publish", False)),
            current=bool(merged.get("current", False)),
        )
    currents = [s.name for s in specs.values() if s.current]
    if len(currents) > 1:
        raise RunsFileError(f"{path}: more than one run sets current: true: {currents}")
    return specs
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: 9 passed

- [ ] **Step 5: Commit**

```bash
git add lib/ni_pangenome.py tests/test_ni_pangenome.py
git commit -m "ni pangenome: runs-file loading (pangenome_runs.yaml)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Checks

**Files:**
- Modify: `lib/ni_pangenome.py` (append)
- Test: `tests/test_ni_pangenome.py` (append)

**Interfaces:**
- Consumes: `RunSpec`, `load_runs_file` (Task 1); `study_runs.load_publish_target(folder_dir: Path) -> PublishTarget | None`.
- Produces:
  - `newick_tips(text: str) -> set[str]`
  - `launch_dir(spec: RunSpec) -> Path` → `spec.study_dir / ".nf_launch" / spec.name`
  - `results_dir(spec: RunSpec) -> Path` → `spec.study_dir / "results" / spec.name`
  - `RUN_RECORD = "ni_run.json"`
  - `check_run(spec: RunSpec, *, assets_root: Path) -> tuple[list[str], list[str]]` → `(errors, warnings)`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_pangenome.py`:

```python
import json


def make_study(tmp_path, *, tree=None, publish_yaml=True, rows=None, **run_overrides):
    """A study that passes every check unless a test breaks one thing."""
    study = tmp_path / "studies" / "fungi" / "s1"
    for sub in ("pep", "dna", "gff3"):
        (study / "data_dir" / sub).mkdir(parents=True)
    rows = rows or [("A", "A.pep.fa", "A.dna.fa", "A.gff3"), ("B", "B.pep.fa", "B.dna.fa", "B.gff3")]
    lines = ["GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup"]
    for short, pep, dna, gff in rows:
        lines.append(f"IN,Sp one,{short},{pep},{dna},{gff},{short},")
        for sub, fname in (("pep", pep), ("dna", dna), ("gff3", gff)):
            if fname:
                (study / "data_dir" / sub / fname).write_text(">x\nM\n")
    (study / "config.csv").write_text("\n".join(lines) + "\n")
    pfam = tmp_path / "Pfam-A.hmm"
    pfam.write_text("HMMER3\n")
    (study / "queue.config").write_text("process {}\n")
    if publish_yaml:
        (study / "publish.yaml").write_text("study: s1\n")
    doc = base_doc(**run_overrides)
    doc["defaults"]["pfam_hmm"] = str(pfam)
    if tree is not None:
        (study / "tree.nwk").write_text(tree)
        doc["runs"]["r1"]["species_tree"] = "tree.nwk"
    write_runs(study, doc)
    return study


def check(study, tmp_path):
    spec = nip.load_runs_file(study)["r1"]
    return nip.check_run(spec, assets_root=tmp_path / "assets")


def test_clean_study_passes(tmp_path):
    study = make_study(tmp_path)
    errors, warnings = check(study, tmp_path)
    assert errors == [] and warnings == []


def test_short_sha_rejected(tmp_path):
    study = make_study(tmp_path, pipeline_commit="0dddb42")
    errors, _ = check(study, tmp_path)
    assert any("40" in e and "0dddb42" in e for e in errors)


def test_missing_samplesheet(tmp_path):
    study = make_study(tmp_path, samplesheet="nope.csv")
    errors, _ = check(study, tmp_path)
    assert any("nope.csv" in e for e in errors)


def test_duplicate_short(tmp_path):
    study = make_study(tmp_path, rows=[("A", "A.pep.fa", "", ""), ("A", "A2.pep.fa", "", "")])
    errors, _ = check(study, tmp_path)
    assert any("duplicate Short" in e and "A" in e for e in errors)


def test_missing_data_file_reported(tmp_path):
    study = make_study(tmp_path)
    (study / "data_dir" / "gff3" / "B.gff3").unlink()
    errors, _ = check(study, tmp_path)
    assert any("gff3/B.gff3" in e for e in errors)


def test_blank_dna_and_gff3_cells_are_skipped(tmp_path):
    study = make_study(tmp_path, rows=[("A", "A.pep.fa", "", ""), ("B", "B.pep.fa", "", "")])
    errors, _ = check(study, tmp_path)
    assert errors == []


def test_tree_tips_must_match_samplesheet(tmp_path):
    study = make_study(tmp_path, tree="(A:0.1,C:0.2);")
    errors, _ = check(study, tmp_path)
    msg = [e for e in errors if "species_tree" in e]
    assert msg and "B" in msg[0] and "C" in msg[0]


def test_newick_tips_ignores_internal_labels_and_lengths():
    text = "(('A x':0.1,B:0.2)95:0.3,(C,D)[&support=1]:0.4)root;"
    assert nip.newick_tips(text) == {"A x", "B", "C", "D"}


def test_matching_tree_passes(tmp_path):
    study = make_study(tmp_path, tree="((A:0.1,B:0.2)100:0.5);")
    errors, _ = check(study, tmp_path)
    assert errors == []


def test_missing_pfam_is_error_unset_pfam_is_warning(tmp_path):
    study = make_study(tmp_path)
    (tmp_path / "Pfam-A.hmm").unlink()
    errors, _ = check(study, tmp_path)
    assert any("pfam_hmm" in e for e in errors)
    study2 = make_study(tmp_path / "b")
    doc = yaml.safe_load((study2 / nip.RUNS_FILE).read_text())
    del doc["defaults"]["pfam_hmm"]
    write_runs(study2, doc)
    errors, warnings = check(study2, tmp_path / "b")
    assert errors == [] and any("#193" in w for w in warnings)


def test_missing_queue_config(tmp_path):
    study = make_study(tmp_path)
    (study / "queue.config").unlink()
    errors, _ = check(study, tmp_path)
    assert any("queue_config" in e for e in errors)


def test_input_inside_nextflow_assets_rejected(tmp_path):
    study = make_study(tmp_path)
    inside = tmp_path / "assets" / ".repos" / "x" / "Pfam-A.hmm"
    inside.parent.mkdir(parents=True)
    inside.write_text("HMMER3\n")
    doc = yaml.safe_load((study / nip.RUNS_FILE).read_text())
    doc["defaults"]["pfam_hmm"] = str(inside)
    write_runs(study, doc)
    errors, _ = check(study, tmp_path)
    assert any("inside" in e and "assets" in e for e in errors)


def test_results_dir_from_other_commit_refused(tmp_path):
    study = make_study(tmp_path)
    launch = study / ".nf_launch" / "r1"
    launch.mkdir(parents=True)
    (launch / nip.RUN_RECORD).write_text(json.dumps({"pipeline_commit": "f" * 40}))
    errors, _ = check(study, tmp_path)
    assert any("different commit" in e and "new run name" in e for e in errors)


def test_same_commit_record_passes_resume(tmp_path):
    study = make_study(tmp_path)
    launch = study / ".nf_launch" / "r1"
    launch.mkdir(parents=True)
    (launch / nip.RUN_RECORD).write_text(json.dumps({"pipeline_commit": SHA}))
    (study / "results" / "r1").mkdir(parents=True)
    errors, _ = check(study, tmp_path)
    assert errors == []


def test_results_dir_without_record_refused(tmp_path):
    study = make_study(tmp_path)
    (study / "results" / "r1").mkdir(parents=True)
    errors, _ = check(study, tmp_path)
    assert any("not made by ni" in e for e in errors)


def test_publish_true_needs_publish_yaml(tmp_path):
    study = make_study(tmp_path, publish_yaml=False, publish=True)
    errors, _ = check(study, tmp_path)
    assert any("publish.yaml" in e for e in errors)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: the 16 new tests FAIL with `AttributeError: module 'ni_pangenome' has no attribute ...`; the 9 Task 1 tests pass.

- [ ] **Step 3: Write the implementation**

Change the imports at the top of `lib/ni_pangenome.py` to:

```python
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import yaml

from study_runs import load_publish_target, validate_slug
```

Append to `lib/ni_pangenome.py`:

```python
RUN_RECORD = "ni_run.json"
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
# A leaf label follows "(" or ",". Internal-node labels follow ")", so they are skipped.
_TIP_RE = re.compile(r"[(,]\s*('(?:[^']|'')*'|[^():,;\s\[\]]+)")
DATA_COLUMNS = (("Protein", "pep"), ("DNA", "dna"), ("GFF3", "gff3"))


def launch_dir(spec: RunSpec) -> Path:
    return spec.study_dir / ".nf_launch" / spec.name


def results_dir(spec: RunSpec) -> Path:
    return spec.study_dir / "results" / spec.name


def newick_tips(text: str) -> set[str]:
    tips = set()
    for label in _TIP_RE.findall(text):
        if label.startswith("'"):
            label = label[1:-1].replace("''", "'")
        tips.add(label)
    return tips


def _is_under(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def check_run(spec: RunSpec, *, assets_root: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []

    if not SHA_RE.match(spec.pipeline_commit):
        errors.append(f"pipeline_commit {spec.pipeline_commit!r} is not a full 40-character "
                      "SHA (nextflow -r cannot fetch a short SHA)")

    shorts: list[str] = []
    if not spec.samplesheet.is_file():
        errors.append(f"samplesheet {spec.samplesheet} not found")
    else:
        with open(spec.samplesheet, newline="") as fh:
            reader = csv.DictReader(fh)
            if "Short" not in (reader.fieldnames or []):
                errors.append(f"samplesheet {spec.samplesheet} has no Short column")
            else:
                rows = list(reader)
                shorts = [r["Short"].strip() for r in rows]
                dups = sorted({s for s in shorts if shorts.count(s) > 1})
                if dups:
                    errors.append(f"samplesheet has duplicate Short value(s): {dups}")
                missing = []
                for r in rows:
                    for col, sub in DATA_COLUMNS:
                        fname = (r.get(col) or "").strip()
                        if fname and not (spec.data_dir / sub / fname).is_file():
                            missing.append(f"{sub}/{fname}")
                if missing:
                    errors.append(f"{len(missing)} data file(s) missing under {spec.data_dir}: "
                                  f"{missing[:10]}")

    if spec.species_tree is not None:
        if not spec.species_tree.is_file():
            errors.append(f"species_tree {spec.species_tree} not found")
        elif shorts:
            tips = newick_tips(spec.species_tree.read_text())
            extra = sorted(tips - set(shorts))
            absent = sorted(set(shorts) - tips)
            if extra or absent:
                errors.append(f"species_tree tips differ from samplesheet Short: "
                              f"extra in tree {extra[:20]}, missing from tree {absent[:20]}")

    if spec.pfam_hmm is None:
        warnings.append("pfam_hmm unset: BUILD_ISLANDS will not run (NovInvenio #193)")
    elif not spec.pfam_hmm.is_file():
        errors.append(f"pfam_hmm {spec.pfam_hmm} not found")

    if spec.queue_config is not None and not spec.queue_config.is_file():
        errors.append(f"queue_config {spec.queue_config} not found")

    for label, p in (("samplesheet", spec.samplesheet), ("data_dir", spec.data_dir),
                     ("pfam_hmm", spec.pfam_hmm), ("queue_config", spec.queue_config),
                     ("species_tree", spec.species_tree)):
        if p is not None and _is_under(p, assets_root):
            errors.append(f"{label} {p} is inside the Nextflow assets dir {assets_root}; "
                          "inputs must live in the study or at an absolute path outside it")

    record = launch_dir(spec) / RUN_RECORD
    if record.is_file():
        prev = json.loads(record.read_text()).get("pipeline_commit")
        if prev != spec.pipeline_commit:
            errors.append(f"results/{spec.name} was made at a different commit ({prev}); "
                          "choose a new run name")
    elif results_dir(spec).exists():
        errors.append(f"results/{spec.name} exists but was not made by ni "
                      f"(no {RUN_RECORD}); choose a new run name")

    if spec.publish and load_publish_target(spec.study_dir) is None:
        errors.append(f"publish: true needs {spec.study_dir}/publish.yaml (study: <name>)")

    return errors, warnings
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: 25 passed

- [ ] **Step 5: Commit**

```bash
git add lib/ni_pangenome.py tests/test_ni_pangenome.py
git commit -m "ni pangenome: run checks (commit, samplesheet, data files, tree tips, commit reuse)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Rendering (params YAML, head-job script, run record)

**Files:**
- Modify: `lib/ni_pangenome.py` (append)
- Create: `tests/data/ni_pangenome/expected_head_job.sh`
- Test: `tests/test_ni_pangenome.py` (append)

**Interfaces:**
- Consumes: `RunSpec`, `launch_dir`, `results_dir`, `RUN_RECORD`.
- Produces:
  - `clone_path(spec: RunSpec, assets_root: Path) -> Path`
  - `build_params(spec: RunSpec) -> dict`
  - `render_params_yaml(spec: RunSpec) -> str`
  - `render_head_job(spec: RunSpec, *, clone: Path, nii_root: Path, extra_args: list[str]) -> str`
  - `build_run_record(spec: RunSpec, *, clone: Path, params_text: str, submitted_at: str) -> dict` — keys `run, pipeline, pipeline_commit, clone_path, samplesheet, samplesheet_sha256, params_sha256, submitted_at, slurm_job_id (None), exit_status (None)`
  - `record_exit(record_path: str, status: int) -> None`
  - `PYTHON = "/usr/bin/python3.12"`

- [ ] **Step 1: Write the golden file**

Create `tests/data/ni_pangenome/expected_head_job.sh` (paths are the fixed fake paths the test uses):

```bash
#!/bin/bash
#SBATCH -J nf-S-r1
#SBATCH -p stajichlab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 30-00:00:00
#SBATCH -o /S/.nf_launch/r1/slurm-%j.out
#SBATCH -e /S/.nf_launch/r1/slurm-%j.err
#
# Generated by `bin/ni pangenome run` -- do not edit; rerun that command instead.
# Study: /S   Run: r1   Pipeline: stajichlab/NovInvenio @ 0dddb421720e0eb0a3045231608b27b113d49be6
set -uo pipefail
export PATH="$HOME/.pixi/bin:$HOME/.local/bin:$PATH"
cd /S/.nf_launch/r1
rc=0
nextflow run stajichlab/NovInvenio -r 0dddb421720e0eb0a3045231608b27b113d49be6 -main-script pangenome.nf \
    -params-file /S/.nf_launch/r1/params.yaml \
    -profile slurm \
    -c /A/clones/0dddb421720e0eb0a3045231608b27b113d49be6/conf/ucr_hpcc_slurm.config \
    -c /S/queue.config \
    -with-trace /S/results/r1/trace.txt \
    -resume || rc=$?
/usr/bin/python3.12 -c "import sys; sys.path.insert(0, '/N/lib'); from ni_pangenome import record_exit; record_exit('/S/.nf_launch/r1/ni_run.json', $rc)"
if [ "$rc" -eq 0 ]; then
    /usr/bin/python3.12 /N/bin/ni pangenome stage --study-dir /S --run r1 || rc=$?
fi
exit $rc
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_ni_pangenome.py`:

```python
import dataclasses

GOLDEN = REPO / "tests" / "data" / "ni_pangenome" / "expected_head_job.sh"


def fake_spec(**changes) -> nip.RunSpec:
    spec = nip.RunSpec(
        name="r1", study_dir=Path("/S"), pipeline="stajichlab/NovInvenio",
        pipeline_commit=SHA, samplesheet=Path("/S/config.csv"), data_dir=Path("/S/data_dir"),
        pfam_hmm=Path("/P/Pfam-A.hmm"), queue_config=Path("/S/queue.config"),
        species_tree=None, head_job=dict(nip.DEFAULT_HEAD_JOB),
        params={"pangenome_rescue_enable": False, "pangenome_cluster_backend": "mmseqs"},
        publish=True, current=False,
    )
    return dataclasses.replace(spec, **changes)


def test_clone_path():
    assert nip.clone_path(fake_spec(), Path("/home/u/.nextflow/assets")) == Path(
        f"/home/u/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/{SHA}")


def test_params_yaml_keeps_booleans_and_absolute_paths():
    data = yaml.safe_load(nip.render_params_yaml(fake_spec(species_tree=Path("/S/t.nwk"))))
    assert data["pangenome_rescue_enable"] is False
    assert data["pangenome_samplesheet"] == "/S/config.csv"
    assert data["pangenome_data_dir"] == "/S/data_dir"
    assert data["outdir"] == "/S/results/r1"
    assert data["pangenome_island_pfam_hmm"] == "/P/Pfam-A.hmm"
    assert data["pangenome_species_tree"] == "/S/t.nwk"
    assert "pangenome_project" not in data


def test_params_yaml_omits_unset_pfam_and_tree():
    data = yaml.safe_load(nip.render_params_yaml(fake_spec(pfam_hmm=None)))
    assert "pangenome_island_pfam_hmm" not in data and "pangenome_species_tree" not in data


def test_head_job_matches_golden():
    text = nip.render_head_job(fake_spec(), clone=Path(f"/A/clones/{SHA}"),
                               nii_root=Path("/N"), extra_args=[])
    assert text == GOLDEN.read_text()


def test_head_job_rules():
    text = nip.render_head_job(fake_spec(publish=False, head_job={**nip.DEFAULT_HEAD_JOB,
                                                                  "account": "exfab"}),
                               clone=Path("/A/c"), nii_root=Path("/N"),
                               extra_args=["-stub", "--x", "a b"])
    assert "BASH_SOURCE" not in text
    assert "#SBATCH -A exfab" in text
    assert "pangenome stage" not in text          # publish: false -> no auto-stage
    assert "-resume -stub --x 'a b' || rc=$?" in text


def test_run_record_and_record_exit(tmp_path):
    ss = tmp_path / "config.csv"
    ss.write_text("Short\nA\n")
    rec = nip.build_run_record(fake_spec(samplesheet=ss), clone=Path("/A/c"),
                               params_text="a: 1\n", submitted_at="2026-09-26T00:00:00Z")
    assert rec["pipeline_commit"] == SHA and rec["slurm_job_id"] is None
    assert rec["exit_status"] is None and len(rec["samplesheet_sha256"]) == 64
    path = tmp_path / nip.RUN_RECORD
    path.write_text(json.dumps(rec))
    nip.record_exit(str(path), 3)
    assert json.loads(path.read_text())["exit_status"] == 3
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: the 6 new tests FAIL with `AttributeError`; the 25 earlier tests pass.

- [ ] **Step 4: Write the implementation**

Add `import hashlib` and `import shlex` to the imports of `lib/ni_pangenome.py`. Append:

```python
PYTHON = "/usr/bin/python3.12"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def clone_path(spec: RunSpec, assets_root: Path) -> Path:
    org, repo = spec.pipeline.split("/", 1)
    return Path(assets_root) / ".repos" / org / repo / "clones" / spec.pipeline_commit


def build_params(spec: RunSpec) -> dict:
    params = dict(spec.params)
    params["pangenome_samplesheet"] = str(spec.samplesheet)
    params["pangenome_data_dir"] = str(spec.data_dir)
    params["outdir"] = str(results_dir(spec))
    if spec.pfam_hmm is not None:
        params["pangenome_island_pfam_hmm"] = str(spec.pfam_hmm)
    if spec.species_tree is not None:
        params["pangenome_species_tree"] = str(spec.species_tree)
    return params


def render_params_yaml(spec: RunSpec) -> str:
    header = (f"# Generated by `bin/ni pangenome run` for run {spec.name}. Do not edit;\n"
              f"# change {RUNS_FILE} and rerun instead.\n")
    return header + yaml.safe_dump(build_params(spec), sort_keys=True)


def render_head_job(spec: RunSpec, *, clone: Path, nii_root: Path, extra_args: list[str]) -> str:
    L = launch_dir(spec)
    hj = spec.head_job
    study_name = spec.study_dir.name
    lines = [
        "#!/bin/bash",
        f"#SBATCH -J nf-{study_name}-{spec.name}",
        f"#SBATCH -p {hj['partition']}",
    ]
    if hj.get("account"):
        lines.append(f"#SBATCH -A {hj['account']}")
    lines += [
        f"#SBATCH -c {hj['cpus']}",
        f"#SBATCH --mem={hj['mem']}",
        f"#SBATCH -t {hj['time']}",
        f"#SBATCH -o {L}/slurm-%j.out",
        f"#SBATCH -e {L}/slurm-%j.err",
        "#",
        "# Generated by `bin/ni pangenome run` -- do not edit; rerun that command instead.",
        f"# Study: {spec.study_dir}   Run: {spec.name}   Pipeline: {spec.pipeline} @ {spec.pipeline_commit}",
        "set -uo pipefail",
        'export PATH="$HOME/.pixi/bin:$HOME/.local/bin:$PATH"',
        f"cd {L}",
        "rc=0",
        f"nextflow run {spec.pipeline} -r {spec.pipeline_commit} -main-script pangenome.nf \\",
        f"    -params-file {L}/params.yaml \\",
        "    -profile slurm \\",
        f"    -c {clone}/conf/ucr_hpcc_slurm.config \\",
    ]
    if spec.queue_config is not None:
        lines.append(f"    -c {spec.queue_config} \\")
    extra = "".join(" " + shlex.quote(a) for a in extra_args)
    lines += [
        f"    -with-trace {results_dir(spec)}/trace.txt \\",
        f"    -resume{extra} || rc=$?",
        (f'{PYTHON} -c "import sys; sys.path.insert(0, \'{nii_root}/lib\'); '
         f"from ni_pangenome import record_exit; record_exit('{L}/{RUN_RECORD}', $rc)\""),
    ]
    if spec.publish:
        lines += [
            'if [ "$rc" -eq 0 ]; then',
            f"    {PYTHON} {nii_root}/bin/ni pangenome stage --study-dir {spec.study_dir} "
            f"--run {spec.name} || rc=$?",
            "fi",
        ]
    lines.append("exit $rc")
    return "\n".join(lines) + "\n"


def build_run_record(spec: RunSpec, *, clone: Path, params_text: str, submitted_at: str) -> dict:
    return {
        "run": spec.name,
        "pipeline": spec.pipeline,
        "pipeline_commit": spec.pipeline_commit,
        "clone_path": str(clone),
        "samplesheet": str(spec.samplesheet),
        "samplesheet_sha256": _sha256_bytes(spec.samplesheet.read_bytes()),
        "params_sha256": _sha256_bytes(params_text.encode()),
        "submitted_at": submitted_at,
        "slurm_job_id": None,
        "exit_status": None,
    }


def record_exit(record_path: str, status: int) -> None:
    p = Path(record_path)
    rec = json.loads(p.read_text())
    rec["exit_status"] = int(status)
    p.write_text(json.dumps(rec, indent=2) + "\n")
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: 31 passed. If `test_head_job_matches_golden` fails, diff the output against the golden file. Fix the code, not the golden file, unless the golden file breaks a Global Constraint.

- [ ] **Step 6: Commit**

```bash
git add lib/ni_pangenome.py tests/test_ni_pangenome.py tests/data/ni_pangenome/expected_head_job.sh
git commit -m "ni pangenome: render params.yaml, head-job script and run record

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Orchestration (`run_pipeline`)

**Files:**
- Modify: `lib/ni_pangenome.py` (append)
- Test: `tests/test_ni_pangenome.py` (append)

**Interfaces:**
- Consumes: everything from Tasks 1-3.
- Produces:
  - `PROVISION_MARKER = ".ni_provisioned"`
  - `run_pipeline(spec: RunSpec, *, nii_root: Path, extra_args: list[str], assets_root: Path, dry_run: bool = False, foreground: bool = False, runner=subprocess.run, now: str | None = None, out=sys.stdout) -> int`

`runner` is called as `runner(cmd: list[str], **kwargs)` and must return an object with `.returncode` and, for `sbatch` and `du`, `.stdout`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_pangenome.py`:

```python
import io


class FakeRunner:
    """Records commands. nextflow pull creates the clone; sbatch returns a job id."""

    def __init__(self, assets_root, *, sbatch_rc=0, sbatch_out="12345\n", pull_rc=0):
        self.calls, self.assets_root = [], assets_root
        self.sbatch_rc, self.sbatch_out, self.pull_rc = sbatch_rc, sbatch_out, pull_rc

    def __call__(self, cmd, **kwargs):
        self.calls.append(cmd)

        class R:
            returncode, stdout = 0, ""
        r = R()
        if cmd[:2] == ["nextflow", "pull"]:
            r.returncode = self.pull_rc
            clone = self.assets_root / ".repos" / "stajichlab" / "NovInvenio" / "clones" / cmd[-1]
            (clone / "conf").mkdir(parents=True, exist_ok=True)
            (clone / "pixi.lock").write_text("lock-v1\n")
        elif cmd[0] == "sbatch":
            r.returncode, r.stdout = self.sbatch_rc, self.sbatch_out
        elif cmd[0] == "du":
            r.stdout = "9.1G\t/x\n"
        return r

    def names(self):
        return [c[0] if c[0] != "nextflow" else "nextflow " + c[1] for c in self.calls]


def ready_spec(tmp_path, **run_overrides):
    study = make_study(tmp_path, **run_overrides)
    return nip.load_runs_file(study)["r1"]


def test_run_order_pull_install_write_submit(tmp_path):
    spec = ready_spec(tmp_path)
    fr = FakeRunner(tmp_path / "assets")
    rc = nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                          runner=fr, now="2026-09-26T00:00:00Z", out=io.StringIO())
    assert rc == 0
    assert fr.names() == ["nextflow pull", "pixi", "du", "sbatch"]
    launch = spec.study_dir / ".nf_launch" / "r1"
    rec = json.loads((launch / nip.RUN_RECORD).read_text())
    assert rec["slurm_job_id"] == "12345"
    assert (launch / "params.yaml").is_file()
    script = launch / "submit_nextflow_head.sh"
    assert script.is_file() and script.stat().st_mode & 0o111


def test_install_skipped_when_marker_matches(tmp_path):
    spec = ready_spec(tmp_path)
    assets = tmp_path / "assets"
    fr = FakeRunner(assets)
    nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=assets, runner=fr,
                     out=io.StringIO())
    fr2 = FakeRunner(assets)
    nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=assets, runner=fr2,
                     out=io.StringIO())
    assert "pixi" not in fr2.names()


def test_dry_run_calls_nothing_and_writes_nothing(tmp_path):
    spec = ready_spec(tmp_path)
    fr = FakeRunner(tmp_path / "assets")
    out = io.StringIO()
    rc = nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                          dry_run=True, runner=fr, out=out)
    assert rc == 0 and fr.calls == []
    assert not (spec.study_dir / ".nf_launch").exists()
    assert "pangenome_samplesheet" in out.getvalue() and "nextflow run" in out.getvalue()


def test_foreground_runs_bash_not_sbatch(tmp_path):
    spec = ready_spec(tmp_path)
    fr = FakeRunner(tmp_path / "assets")
    nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                     foreground=True, runner=fr, out=io.StringIO())
    assert "sbatch" not in fr.names() and fr.calls[-1][0] == "bash"


def test_check_failure_stops_before_any_command(tmp_path):
    spec = ready_spec(tmp_path, pipeline_commit="0dddb42")
    fr = FakeRunner(tmp_path / "assets")
    rc = nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                          runner=fr, out=io.StringIO())
    assert rc == 1 and fr.calls == []


def test_pull_failure_stops(tmp_path):
    spec = ready_spec(tmp_path)
    fr = FakeRunner(tmp_path / "assets", pull_rc=1)
    rc = nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                          runner=fr, out=io.StringIO())
    assert rc == 1 and fr.names() == ["nextflow pull"]


def test_sbatch_failure_leaves_job_id_null(tmp_path):
    spec = ready_spec(tmp_path)
    fr = FakeRunner(tmp_path / "assets", sbatch_rc=1, sbatch_out="")
    rc = nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                          runner=fr, out=io.StringIO())
    assert rc != 0
    rec = json.loads((spec.study_dir / ".nf_launch" / "r1" / nip.RUN_RECORD).read_text())
    assert rec["slurm_job_id"] is None


def test_sbatch_garbage_output_is_failure(tmp_path):
    spec = ready_spec(tmp_path)
    fr = FakeRunner(tmp_path / "assets", sbatch_out="Submitted batch job\n")
    rc = nip.run_pipeline(spec, nii_root=REPO, extra_args=[], assets_root=tmp_path / "assets",
                          runner=fr, out=io.StringIO())
    assert rc != 0
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: the 8 new tests FAIL with `AttributeError: ... run_pipeline`.

- [ ] **Step 3: Write the implementation**

Add `import datetime`, `import subprocess`, `import sys` and `import time` to the imports. Append:

```python
PROVISION_MARKER = ".ni_provisioned"


def _now_utc() -> str:
    return datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_pipeline(spec: RunSpec, *, nii_root: Path, extra_args: list[str], assets_root: Path,
                 dry_run: bool = False, foreground: bool = False, runner=subprocess.run,
                 now: str | None = None, out=sys.stdout) -> int:
    errors, warnings = check_run(spec, assets_root=assets_root)
    for w in warnings:
        print(f"WARNING: {w}", file=out)
    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=out)
        return 1

    clone = clone_path(spec, assets_root)
    params_text = render_params_yaml(spec)
    script_text = render_head_job(spec, clone=clone, nii_root=nii_root, extra_args=extra_args)
    if dry_run:
        print(f"== params.yaml ==\n{params_text}\n== submit_nextflow_head.sh ==\n{script_text}",
              file=out)
        return 0

    pull = runner(["nextflow", "pull", spec.pipeline, "-r", spec.pipeline_commit])
    if pull.returncode != 0:
        print(f"ERROR: nextflow pull of {spec.pipeline_commit} failed "
              "(is the commit pushed to GitHub?)", file=out)
        return pull.returncode

    lock = clone / "pixi.lock"
    lock_sha = _sha256_bytes(lock.read_bytes()) if lock.is_file() else "no-lock"
    marker = clone / PROVISION_MARKER
    if not (marker.is_file() and marker.read_text().strip() == lock_sha):
        t0 = time.monotonic()
        inst = runner(["pixi", "install", "--frozen", "--manifest-path", str(clone / "pixi.toml")])
        if inst.returncode != 0:
            print("ERROR: pixi install failed in the pipeline clone", file=out)
            return inst.returncode
        size = runner(["du", "-sh", str(clone / ".pixi")], capture_output=True, text=True)
        print(f"== pixi install: {time.monotonic() - t0:.0f} s; env size "
              f"{(size.stdout or '?').split()[0] if size.stdout else '?'} ==", file=out)
        marker.write_text(lock_sha + "\n")

    L = launch_dir(spec)
    L.mkdir(parents=True, exist_ok=True)
    (L / "params.yaml").write_text(params_text)
    script = L / "submit_nextflow_head.sh"
    script.write_text(script_text)
    script.chmod(0o755)
    record = build_run_record(spec, clone=clone, params_text=params_text,
                              submitted_at=now or _now_utc())
    record_path = L / RUN_RECORD
    record_path.write_text(json.dumps(record, indent=2) + "\n")

    if foreground:
        return runner(["bash", str(script)]).returncode

    sub = runner(["sbatch", "--parsable", str(script)], capture_output=True, text=True)
    job_id = (sub.stdout or "").strip().split(";")[0]
    if sub.returncode != 0 or not job_id.isdigit():
        print(f"ERROR: sbatch failed (rc={sub.returncode}, output={sub.stdout!r})", file=out)
        return sub.returncode or 1
    record["slurm_job_id"] = job_id
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    print(f"== submitted {spec.name} as job {job_id}; log {L}/slurm-{job_id}.out ==", file=out)
    return 0
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: 39 passed

- [ ] **Step 5: Commit**

```bash
git add lib/ni_pangenome.py tests/test_ni_pangenome.py
git commit -m "ni pangenome: run_pipeline (pull, one-time pixi install, head job, sbatch)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: CLI wiring in `bin/ni` (list, check, run, stage, publish)

**Files:**
- Modify: `lib/ni_pangenome.py` (append)
- Modify: `bin/ni`
- Test: `tests/test_ni_pangenome.py` (append)

**Interfaces:**
- Consumes: Tasks 1-4; `study_runs.load_publish_target`.
- Produces:
  - `run_state(spec: RunSpec, repo_root: Path) -> str` → one of `not started`, `submitted <jobid>`, `failed (exit N)`, `done`, `staged`
  - `stage_command(spec: RunSpec, repo_root: Path) -> list[str]`
  - `publish_command(spec: RunSpec, repo_root: Path) -> list[str]` (raises `RunsFileError` if not allowed)
  - `bin/ni pangenome <list|check|run|stage|publish> --study-dir DIR [--run R] [--dry-run] [--foreground] [-- extra]`

`repo_root` is the NII repo root (`BIN.parent` in `bin/ni`). Study `<domain>/<folder>` is `spec.study_dir` relative to `repo_root / "studies"`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_pangenome.py`:

```python
import importlib.machinery
import importlib.util


def load_ni_cli():
    loader = importlib.machinery.SourceFileLoader("ni_cli", str(REPO / "bin" / "ni"))
    spec = importlib.util.spec_from_file_location("ni_cli", REPO / "bin" / "ni", loader=loader)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def repo_study(tmp_path, **kw):
    """make_study puts the study at <tmp>/studies/fungi/s1, so <tmp> acts as repo root."""
    study = make_study(tmp_path, **kw)
    return tmp_path, nip.load_runs_file(study)["r1"]


def test_stage_command(tmp_path):
    root, spec = repo_study(tmp_path, current=True)
    cmd = nip.stage_command(spec, root)
    assert cmd == [nip.PYTHON, str(root / "bin" / "sync_pangenome_report.py"),
                   "--study", "fungi/s1", "--run", "r1", "--current"]


def test_publish_refuses_unstaged_and_unpublishable(tmp_path):
    root, spec = repo_study(tmp_path, publish=True)
    with pytest.raises(nip.RunsFileError, match="stage"):
        nip.publish_command(spec, root)
    root2, spec2 = repo_study(tmp_path / "b", publish=False)
    with pytest.raises(nip.RunsFileError, match="publish: true"):
        nip.publish_command(spec2, root2)


def test_publish_command_when_staged(tmp_path):
    root, spec = repo_study(tmp_path, publish=True)
    run_docs = root / "docs" / "fungi" / "s1" / "r1"
    run_docs.mkdir(parents=True)
    (run_docs / "run.json").write_text("{}")
    assert nip.publish_command(spec, root) == [
        str(root / "bin" / "publish_report_release.sh"), "fungi/s1/r1"]


def test_run_state(tmp_path):
    root, spec = repo_study(tmp_path, publish=True)
    assert nip.run_state(spec, root) == "not started"
    launch = spec.study_dir / ".nf_launch" / "r1"
    launch.mkdir(parents=True)
    rec = {"pipeline_commit": SHA, "slurm_job_id": "77", "exit_status": None}
    (launch / nip.RUN_RECORD).write_text(json.dumps(rec))
    assert nip.run_state(spec, root) == "submitted 77"
    rec["exit_status"] = 2
    (launch / nip.RUN_RECORD).write_text(json.dumps(rec))
    assert nip.run_state(spec, root) == "failed (exit 2)"
    rec["exit_status"] = 0
    (launch / nip.RUN_RECORD).write_text(json.dumps(rec))
    assert nip.run_state(spec, root) == "done"
    (root / "docs" / "fungi" / "s1" / "r1").mkdir(parents=True)
    (root / "docs" / "fungi" / "s1" / "r1" / "run.json").write_text("{}")
    assert nip.run_state(spec, root) == "staged"


def test_cli_check_exit_codes(tmp_path, monkeypatch, capsys):
    ni = load_ni_cli()
    root, spec = repo_study(tmp_path)
    monkeypatch.setattr(sys, "argv", ["ni", "pangenome", "check", "--study-dir",
                                      str(spec.study_dir), "--run", "r1"])
    assert ni.main() == 0
    (spec.study_dir / "queue.config").unlink()
    assert ni.main() == 1
    assert "queue_config" in capsys.readouterr().out


def test_cli_list(tmp_path, monkeypatch, capsys):
    ni = load_ni_cli()
    root, spec = repo_study(tmp_path)
    monkeypatch.setattr(sys, "argv", ["ni", "pangenome", "list", "--study-dir",
                                      str(spec.study_dir)])
    assert ni.main() == 0
    out = capsys.readouterr().out
    assert "r1" in out and SHA[:7] in out and "not started" in out


def test_cli_run_dry_run_forwards_extra_args(tmp_path, monkeypatch, capsys):
    ni = load_ni_cli()
    root, spec = repo_study(tmp_path)
    monkeypatch.setattr(sys, "argv", ["ni", "pangenome", "run", "--study-dir",
                                      str(spec.study_dir), "--run", "r1", "--dry-run",
                                      "--", "-stub"])
    assert ni.main() == 0
    assert "-resume -stub" in capsys.readouterr().out


def test_cli_unknown_run(tmp_path, monkeypatch, capsys):
    ni = load_ni_cli()
    root, spec = repo_study(tmp_path)
    monkeypatch.setattr(sys, "argv", ["ni", "pangenome", "check", "--study-dir",
                                      str(spec.study_dir), "--run", "nope"])
    assert ni.main() == 1
    assert "nope" in capsys.readouterr().out
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py`
Expected: the 8 new tests FAIL.

- [ ] **Step 3: Implement the lib helpers**

Append to `lib/ni_pangenome.py`:

```python
def _study_rel(spec: RunSpec, repo_root: Path) -> str:
    try:
        return str(spec.study_dir.resolve().relative_to((repo_root / "studies").resolve()))
    except ValueError:
        raise RunsFileError(f"{spec.study_dir} is not under {repo_root}/studies") from None


def _run_docs(spec: RunSpec, repo_root: Path) -> Path | None:
    target = load_publish_target(spec.study_dir)
    if target is None:
        return None
    try:
        domain = _study_rel(spec, repo_root).split("/", 1)[0]
    except RunsFileError:
        return None  # study outside repo_root/studies: it cannot have site pages
    return repo_root / "docs" / domain / target.study / spec.name


def run_state(spec: RunSpec, repo_root: Path) -> str:
    docs = _run_docs(spec, repo_root)
    if docs is not None and (docs / "run.json").is_file():
        return "staged"
    record = launch_dir(spec) / RUN_RECORD
    if not record.is_file():
        return "not started"
    rec = json.loads(record.read_text())
    status = rec.get("exit_status")
    if status is None:
        return f"submitted {rec.get('slurm_job_id')}"
    return "done" if status == 0 else f"failed (exit {status})"


def stage_command(spec: RunSpec, repo_root: Path) -> list[str]:
    cmd = [PYTHON, str(repo_root / "bin" / "sync_pangenome_report.py"),
           "--study", _study_rel(spec, repo_root), "--run", spec.name]
    if spec.current:
        cmd.append("--current")
    return cmd


def publish_command(spec: RunSpec, repo_root: Path) -> list[str]:
    if not spec.publish:
        raise RunsFileError(f"run {spec.name} does not set publish: true")
    docs = _run_docs(spec, repo_root)
    if docs is None:
        raise RunsFileError(f"{spec.study_dir}/publish.yaml not found")
    if not (docs / "run.json").is_file():
        raise RunsFileError(f"run {spec.name} is not staged ({docs}/run.json missing); "
                            "run `ni pangenome stage` first")
    domain = _study_rel(spec, repo_root).split("/", 1)[0]
    target = load_publish_target(spec.study_dir)
    return [str(repo_root / "bin" / "publish_report_release.sh"),
            f"{domain}/{target.study}/{spec.name}"]
```

- [ ] **Step 4: Wire the CLI in `bin/ni`**

In `bin/ni`, add to the module docstring usage block:

```
  bin/ni pangenome list    --study-dir studies/<domain>/<set>
  bin/ni pangenome check   --study-dir studies/<domain>/<set> --run <run>
  bin/ni pangenome run     --study-dir studies/<domain>/<set> --run <run> [--dry-run] [--foreground] [-- extra nextflow args]
  bin/ni pangenome stage   --study-dir studies/<domain>/<set> --run <run>
  bin/ni pangenome publish --study-dir studies/<domain>/<set> --run <run>
```

After the `run_ap` parser block (before `args = ap.parse_args()`), add:

```python
    pg_ap = sub.add_parser("pangenome", help="Run pangenome.nf for a study by run name")
    pg_sub = pg_ap.add_subparsers(dest="pg_command", required=True)
    for name in ("list", "check", "run", "stage", "publish"):
        p = pg_sub.add_parser(name)
        p.add_argument("--study-dir", required=True)
        if name != "list":
            p.add_argument("--run", required=True)
        if name == "run":
            p.add_argument("--dry-run", action="store_true")
            p.add_argument("--foreground", action="store_true")
            p.add_argument("extra_args", nargs=argparse.REMAINDER)
```

Before the final `return 1` in `main()`, add:

```python
    if args.command == "pangenome":
        return _pangenome(args)
```

Add this function above `main()`:

```python
def _pangenome(args) -> int:
    import ni_pangenome as nip

    repo_root = BIN.parent
    assets_root = Path.home() / ".nextflow" / "assets"
    try:
        specs = nip.load_runs_file(Path(args.study_dir).resolve())
    except nip.RunsFileError as e:
        print(f"ERROR: {e}")
        return 1
    if args.pg_command == "list":
        for spec in specs.values():
            n = "?"
            if spec.samplesheet.is_file():
                n = str(max(sum(1 for _ in open(spec.samplesheet)) - 1, 0))
            print(f"{spec.name}\t{spec.pipeline_commit[:7]}\t{n} rows\t"
                  f"{nip.run_state(spec, repo_root)}")
        return 0
    spec = specs.get(args.run)
    if spec is None:
        print(f"ERROR: run {args.run!r} not in {nip.RUNS_FILE} (have: {sorted(specs)})")
        return 1
    try:
        if args.pg_command == "check":
            errors, warnings = nip.check_run(spec, assets_root=assets_root)
            for w in warnings:
                print(f"WARNING: {w}")
            for e in errors:
                print(f"ERROR: {e}")
            if not errors:
                print(f"OK: {spec.name}")
            return 1 if errors else 0
        if args.pg_command == "run":
            extra = [a for a in args.extra_args if a != "--"]
            return nip.run_pipeline(spec, nii_root=repo_root, extra_args=extra,
                                    assets_root=assets_root, dry_run=args.dry_run,
                                    foreground=args.foreground)
        if args.pg_command == "stage":
            return subprocess.run(nip.stage_command(spec, repo_root)).returncode
        if args.pg_command == "publish":
            return subprocess.run(nip.publish_command(spec, repo_root)).returncode
    except nip.RunsFileError as e:
        print(f"ERROR: {e}")
        return 1
    return 1
```

- [ ] **Step 5: Run all ni tests**

Run: `pixi run python -m pytest -q tests/test_ni_pangenome.py tests/test_ni_run.py tests/test_ni_resolve.py tests/test_ni_discover.py`
Expected: 47 passed in `test_ni_pangenome.py`, plus the 72 existing ni tests (119 passed in total)

- [ ] **Step 6: Commit**

```bash
git add bin/ni lib/ni_pangenome.py tests/test_ni_pangenome.py
git commit -m "ni pangenome: CLI (list, check, run, stage, publish)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: `sync_pangenome_report.py` records the pipeline commit

**Files:**
- Modify: `bin/sync_pangenome_report.py` (`stage_run` and its call in `main`)
- Test: `tests/test_sync_pangenome_report.py` (append)

**Interfaces:**
- Consumes: `.nf_launch/<run>/ni_run.json` written by Task 4 (keys `pipeline`, `pipeline_commit`).
- Produces: `stage_run(..., today: str, run_record: Path | None = None) -> dict`. When `run_record` is an existing file, `run.json` gains `"pipeline_repo"` and `"pipeline_commit"`. `"pipeline": "pangenome.nf"` stays unchanged, because the gallery reads it.

- [ ] **Step 1: Write the failing tests**

`tests/test_sync_pangenome_report.py` already has `_fake_study(root, run="run_a", ...)` (builds
`<root>/studies/fungi/demo_pangenome/` with a run under `results/run_a/output/pangenome/`) and
`_run(root, *extra)` (calls `main(["--study", "fungi/demo_pangenome", "--repo-root", root, ...])`).
Append:

```python
def test_run_json_gets_pipeline_commit_from_ni_run_record(tmp_path):
    study = _fake_study(tmp_path)
    launch = study / ".nf_launch" / "run_a"
    launch.mkdir(parents=True)
    (launch / "ni_run.json").write_text(json.dumps(
        {"pipeline": "stajichlab/NovInvenio", "pipeline_commit": "a" * 40}))
    assert _run(tmp_path, "--run", "run_a") == 0
    meta = json.loads((tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "run.json").read_text())
    assert meta["pipeline"] == "pangenome.nf"
    assert meta["pipeline_repo"] == "stajichlab/NovInvenio"
    assert meta["pipeline_commit"] == "a" * 40


def test_run_json_without_ni_run_record_is_unchanged(tmp_path):
    _fake_study(tmp_path)
    assert _run(tmp_path, "--run", "run_a") == 0
    meta = json.loads((tmp_path / "docs" / "fungi" / "demo_pangenome" / "run_a" / "run.json").read_text())
    assert "pipeline_commit" not in meta and "pipeline_repo" not in meta
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m pytest -q tests/test_sync_pangenome_report.py`
Expected: `test_run_json_gets_pipeline_commit_from_ni_run_record` FAILS with `KeyError: 'pipeline_repo'`; the other test passes.

- [ ] **Step 3: Implement**

In `bin/sync_pangenome_report.py`, change the `stage_run` signature to:

```python
def stage_run(pangenome_dir: Path, run_docs: Path, domain: str, set_name: str, run: str,
              source_label: str, today: str, run_record: Path | None = None) -> dict:
```

Just before `(run_docs / "run.json").write_text(...)`, insert:

```python
    if run_record is not None and run_record.is_file():
        rec = json.loads(run_record.read_text())
        meta["pipeline_repo"] = rec.get("pipeline")
        meta["pipeline_commit"] = rec.get("pipeline_commit")
```

In `main()`, change the `stage_run(...)` call to pass the record:

```python
    meta = stage_run(pangenome_dir, study_docs / args.run, domain, set_name, args.run,
                     source_label, today,
                     run_record=study_dir / ".nf_launch" / args.run / "ni_run.json")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run python -m pytest -q tests/test_sync_pangenome_report.py tests/test_ni_pangenome.py`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add bin/sync_pangenome_report.py tests/test_sync_pangenome_report.py
git commit -m "sync_pangenome_report: record pipeline commit from ni_run.json in run.json

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Live verification (controller runs this; not a subagent task)

**Files:**
- Create: `studies/fungi/fusarium_FOL/pangenome_runs.yaml`
- Create: `studies/fungi/coccidioides_pangenome/pangenome_runs.yaml`

- [ ] **Step 1: Write the FOL runs file**

Copy the non-rescue params from `studies/fungi/fusarium_FOL/params_norescue.yaml`, and the head-job partition and account from its `.nf_launch/mmseqs_norescue/submit_nextflow_head.sh` (`-p exfab -A exfab`):

```yaml
# pangenome.nf runs for fusarium_FOL, launched with `bin/ni pangenome run`.
# Spec: notes/superpowers/specs/2026-09-26-ni-pangenome-mode-design.md
defaults:
  pipeline_commit: <full SHA of NovInvenio origin/main at run time>
  pfam_hmm: /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm
  queue_config: stajichlab_queue.config
  head_job: {partition: exfab, account: exfab, time: 30-00:00:00, mem: 8G, cpus: 1}
  params:
    pangenome_cluster_backend: mmseqs
    pangenome_pfam_chunk_size: 2000
    pangenome_mash_threshold: 0.0005
runs:
  ni_smoke_norescue:
    samplesheet: config.csv
    params: {pangenome_rescue_enable: false}
```

Fill `<full SHA ...>` with the output of `git -C /bigdata/stajichlab/jstajich/projects/NovInvenio rev-parse origin/main` after `git fetch`. The value must be a real 40-character SHA before this file is committed.

- [ ] **Step 2: Check and dry-run**

Run: `bin/ni pangenome check --study-dir studies/fungi/fusarium_FOL --run ni_smoke_norescue`
Expected: `WARNING`-free and `OK: ni_smoke_norescue`.
Run: `bin/ni pangenome run --study-dir studies/fungi/fusarium_FOL --run ni_smoke_norescue --dry-run`
Expected: `params.yaml` and the head-job script print. The script has no `BASH_SOURCE`, and `stage` is absent because `publish` is false.

- [ ] **Step 3: Real pull, install and submit**

Run: `bin/ni pangenome run --study-dir studies/fungi/fusarium_FOL --run ni_smoke_norescue`
Record: the pixi install time and env size printed by `ni`, and the job ID.

- [ ] **Step 4: Watch the job to completion**

Use a background until-loop on `ni pangenome list` (state leaves `submitted`). On `done`: compare the family count (`wc -l frequency_table.tsv`) and presence-matrix size of `results/ni_smoke_norescue/output/pangenome/` with `results/mmseqs_norescue/output/pangenome/` (pipeline `af6fd68`). Report the numbers as measured. On `failed`: report the Nextflow error from `.nf_launch/ni_smoke_norescue/.nextflow.log`.

- [ ] **Step 5: Cocci runs file (dry-run only)**

Write `studies/fungi/coccidioides_pangenome/pangenome_runs.yaml` with the four planned runs: both reciprocal directions (with the 529-tip species tree, rescue on, `publish: true`, one `current: true`), plus immitis-only and posadasii-only (single-species samplesheets `config_immitis.csv` / `config_posadasii.csv`, no species tree). Leave `pipeline_commit` at the current origin/main SHA. It is bumped to the locus-view commit before the real runs. Run `check` and `--dry-run` on all four and record the output.

- [ ] **Step 6: Commit and push**

```bash
git add studies/fungi/fusarium_FOL/pangenome_runs.yaml studies/fungi/coccidioides_pangenome/pangenome_runs.yaml
git commit -m "fusarium_FOL, coccidioides: pangenome_runs.yaml for ni pangenome

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push origin main
```

On a successful FOL smoke run, remove `studies/fungi/fusarium_FOL/run_pangenome.sh` in a separate commit (spec section 6). Keep `params_*.yaml` if other notes cite them.
