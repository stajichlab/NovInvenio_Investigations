"""ni pangenome: run nf_NovInvenio's pangenome.nf for an NII study by run name.

Spec: notes/superpowers/specs/2026-09-26-ni-pangenome-mode-design.md.
A study declares its runs in studies/<domain>/<set>/pangenome_runs.yaml.
"""
from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import yaml

from study_runs import load_publish_target, validate_slug

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
