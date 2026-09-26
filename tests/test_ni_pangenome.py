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
