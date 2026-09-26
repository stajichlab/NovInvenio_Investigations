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
