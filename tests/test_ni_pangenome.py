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
