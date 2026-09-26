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
