"""ni pangenome: run nf_NovInvenio's pangenome.nf for an NII study by run name.

Spec: notes/superpowers/specs/2026-09-26-ni-pangenome-mode-design.md.
A study declares its runs in studies/<domain>/<set>/pangenome_runs.yaml.
"""
from __future__ import annotations

import csv
import datetime
import hashlib
import json
import re
import shlex
import subprocess
import sys
import time
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
        "module load java 2>/dev/null || true",
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

    record_path = launch_dir(spec) / RUN_RECORD
    if record_path.is_file():
        prev_rec = json.loads(record_path.read_text())
        job_id = prev_rec.get("slurm_job_id")
        if isinstance(job_id, str) and job_id.isdigit() and prev_rec.get("exit_status") is None:
            sq = runner(["squeue", "-h", "-j", job_id, "-o", "%T"], capture_output=True, text=True)
            state = (sq.stdout or "").strip()
            if sq.returncode == 0 and state:
                print(f"ERROR: run {spec.name} is SLURM job {job_id} ({state}); scancel it "
                      "first, or wait for it to end", file=out)
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
    if foreground:
        record["slurm_job_id"] = "foreground"
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
