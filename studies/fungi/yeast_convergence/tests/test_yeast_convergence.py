"""Tests for studies/fungi/yeast_convergence/bin/ (run by `pixi run test`, which
collects studies/). Paths are derived from this file so the suite runs from any
checkout; the scripts themselves default to the cluster path and take overrides."""
import csv
import gzip
import importlib.util
import shutil
from pathlib import Path

import pytest
import yaml

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
NII = FOLDER.parents[2]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, FOLDER / "bin" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


yc = _load("yc_batches")
cc = _load("convergence_compare")
GROWTH = yc.load_growth_forms(NII / "config_support/traits/trait_definitions.yaml")


# ---------- batch specs ----------

def test_committed_batches_pass_check():
    batches = yc.load_batches(FOLDER / "batches")
    assert yc.check(batches, GROWTH) == []


def test_committed_batches_are_novinvenio_compatible():
    """Same required keys build_targeted_configs.py reads."""
    for p in sorted((FOLDER / "batches").glob("*.yaml")):
        spec = yaml.safe_load(p.read_text())
        assert isinstance(spec["batch"], str)
        pools = spec["outgroup_pools"]
        for s in spec["studies"]:
            assert isinstance(s["focal"], str)
            assert s["outgroup_pool"] in pools
            assert s["ingroup_extra"]["mode"] == "explicit"
            assert all(isinstance(m, str) for m in s["ingroup_extra"]["members"])


def _mini_batch(tmp_path, **study_overrides):
    study = {"id": "yc_test", "origin": "x", "contrast": "yeast_vs_filamentous_sister",
             "focal": "Yeast one",
             "ingroup_extra": {"mode": "explicit", "members": [], "reason": "r"},
             "outgroup_pool": "pool"}
    study.update(study_overrides)
    spec = {"batch": "b", "outgroup_pools": {"pool": {"members": ["Fil one"]}},
            "studies": [study],
            "species": {"Yeast one": {"taxon_group": "T", "growth_form": "yeast"},
                        "Fil one": {"taxon_group": "T", "growth_form": "filamentous"},
                        "Dimo hyph": {"taxon_group": "T", "growth_form": "dimorphic",
                                      "dominant_phase": "hyphal"}}}
    d = tmp_path / "batches"
    d.mkdir(exist_ok=True)
    (d / "b.yaml").write_text(yaml.safe_dump(spec))
    return yc.load_batches(d)


def test_check_accepts_minimal(tmp_path):
    assert yc.check(_mini_batch(tmp_path), GROWTH) == []


def test_check_rejects_focal_in_outgroup(tmp_path):
    probs = yc.check(_mini_batch(tmp_path, focal="Fil one"), GROWTH)
    assert any("both IN and OUT" in p for p in probs)


def test_check_rejects_origin_without_yeast(tmp_path):
    probs = yc.check(_mini_batch(tmp_path, focal="Dimo hyph"), GROWTH)
    assert any("origin IN has no yeast_dominant" in p for p in probs)


def test_capacity_origin_accepts_hyphal_dimorph(tmp_path):
    b = _mini_batch(tmp_path, focal="Dimo hyph", contrast="dimorphic_vs_filamentous_sister")
    assert yc.check(b, GROWTH) == []


def test_check_rejects_reversal_without_compare_with(tmp_path):
    b = _mini_batch(tmp_path, focal="Fil one", contrast="reversal_x")
    probs = yc.check(b, GROWTH)
    assert any("compare_with" in p for p in probs)


def test_life_state():
    assert yc.life_state("filamentous", None) == "filamentous"
    assert yc.life_state("yeast", None) == "yeast_dominant"
    assert yc.life_state("dimorphic", "yeast") == "yeast_dominant"
    assert yc.life_state("dimorphic", "hyphal") == "dimorphic_hyphal"


def test_render_prefills_and_protects(tmp_path):
    batches = _mini_batch(tmp_path)
    known = {"Yeast one": {"Short": "Yone", "Strain": "S1", "Protein_Source": "uniprot",
                           "Protein_Accession": "UP0001", "Taxon_ID": "1",
                           "Genome_Source": "ncbi", "Genome_Accession": "GCA_1"}}
    index = yc.render(batches, tmp_path, known, force=False)
    assert index[0]["n_prefilled"] == 1 and index[0]["n_to_resolve"] == 1
    rows = list(csv.DictReader(open(tmp_path / "yc_test/species.csv")))
    assert rows[0]["Short"] == "Yone" and rows[0]["Protein_Accession"] == "UP0001"
    assert rows[1]["Short"] == "Filone" and rows[1]["Protein_Source"] == ""
    assert rows[1]["Group"] == "OUT"
    # a resolved/hand-edited species.csv is not overwritten without --force
    (tmp_path / "yc_test/species.csv").write_text("edited\n")
    yc.render(batches, tmp_path, known, force=False)
    assert (tmp_path / "yc_test/species.csv").read_text() == "edited\n"


def test_committed_species_csv_match_batches(tmp_path):
    """The committed sub-study species.csv files are what render would write,
    unless a row has since been filled by bin/ni resolve (non-blank source)."""
    batches = yc.load_batches(FOLDER / "batches")
    known = yc.known_sources(NII / "studies", exclude=FOLDER)
    yc.render(batches, tmp_path, known, force=True)
    for sub in sorted(tmp_path.glob("yc_*")):
        fresh = list(csv.DictReader(open(sub / "species.csv")))
        committed = list(csv.DictReader(open(FOLDER / sub.name / "species.csv")))
        assert [r["Species"] for r in fresh] == [r["Species"] for r in committed]
        assert [r["Group"] for r in fresh] == [r["Group"] for r in committed]


# ---------- convergence comparison ----------

def test_poisson_binomial_sf():
    assert cc.poisson_binomial_sf([0.5, 0.5], 0) == 1.0
    assert cc.poisson_binomial_sf([0.5, 0.5], 2) == pytest.approx(0.25)
    assert cc.poisson_binomial_sf([0.5, 0.5], 1) == pytest.approx(0.75)
    assert cc.poisson_binomial_sf([0.1, 0.2, 0.3], 3) == pytest.approx(0.006)


def test_bh_monotone():
    q = cc.bh([0.01, 0.04, 0.03, 0.5])
    # sorted p: .01 .03 .04 .5 -> raw p*m/rank: .04 .06 .0533 .5 -> step-up minimum
    assert q == pytest.approx([0.04, 0.04 / 0.75, 0.04 / 0.75, 0.5])


HEADER = ["scope", "level", "vocab", "term", "name", "k_cand_with_term", "N_cand_annotated",
          "n_bg_with_term", "M_bg_annotated", "expected", "fold", "pvalue", "fdr_bh"]


def _table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEADER)
        for term, k, fold, q in rows:
            w.writerow(["whole_set", "family", "pfam", term, f"{term} name", k, 100, 10, 1000,
                        1, fold, q, q])


def _synthetic_folder(tmp_path):
    """3 origins share PF00001 (lost) and PF00172 only in one; the pseudo-origin
    shares PF00002 with the origins; reversal regains PF00001."""
    index = [
        ("yc_o1", "origin", ""), ("yc_o2", "origin", ""), ("yc_o3", "origin", ""),
        ("yc_p1", "pseudo", ""), ("yc_r1", "reversal", "yc_o1"),
    ]
    with open(tmp_path / "study_index.tsv", "w") as fh:
        fh.write("id\tbatch\trole\torigin\tcontrast\tcompare_with\tn_in\tn_out\tn_prefilled\tn_to_resolve\n")
        for i, role, cw in index:
            fh.write(f"{i}\tb\t{role}\tx\tc\t{cw}\t1\t1\t0\t0\n")
    shutil.copy(FOLDER / "focus_terms.tsv", tmp_path / "focus_terms.tsv")
    background = [(f"PF9{n:04d}", 1, 0.5, 1.0) for n in range(40)]
    sig_loss = {"yc_o1": ["PF00001", "PF00002", "PF00172"], "yc_o2": ["PF00001", "PF00002"],
                "yc_o3": ["PF00001", "PF00002"], "yc_p1": ["PF00002"], "yc_r1": []}
    for i, _, _ in index:
        base = tmp_path / i / "analysis/enrichment/outputs"
        loss = [(t, 5, 4.0, 0.001) for t in sig_loss[i]]
        loss += [(t, 1, 0.5, 1.0) for t in ("PF00001", "PF00002", "PF00172") if t not in sig_loss[i]]
        _table(base / "loss/proteome_clustered/enrichment_domains.tsv.gz", loss + background)
        nov = [("PF00001", 4, 3.0, 0.01)] if i == "yc_r1" else [("PF00001", 0, 0.0, 1.0)]
        _table(base / "novelty/proteome_clustered/enrichment_domains.tsv.gz", nov + background)
    return tmp_path


def test_convergence_end_to_end(tmp_path):
    folder = _synthetic_folder(tmp_path)
    out = tmp_path / "out"
    assert cc.main(["--folder", str(folder), "--out-dir", str(out), "--vocab", "pfam"]) == 0
    rec = {(r["direction"], r["term"]): r
           for r in csv.DictReader(open(out / "term_recurrence.tsv"), delimiter="\t")}
    assert rec[("loss", "PF00001")]["class"] == "yeast_convergent"
    assert rec[("loss", "PF00002")]["class"] == "streamlining_generic"
    assert rec[("loss", "PF00172")]["class"] == ""          # one origin only
    assert int(rec[("loss", "PF00001")]["n_origin_sig"]) == 3
    rev = list(csv.DictReader(open(out / "reversal_overlap.tsv"), delimiter="\t"))
    assert [(r["term"], r["pattern"]) for r in rev] == [
        ("PF00001", "lost_in_yeast_origin__gained_in_reversal")]
    focus = list(csv.DictReader(open(out / "focus_terms.tsv"), delimiter="\t"))
    zn = [f for f in focus if f["term"] == "PF00172" and f["direction"] == "loss"]
    assert {f["id"]: f["status"] for f in zn}["yc_o1"] == "significant"
    assert "yeast_convergent" in (out / "summary.md").read_text()


def test_convergence_reports_missing(tmp_path):
    folder = _synthetic_folder(tmp_path)
    shutil.rmtree(folder / "yc_o3/analysis/enrichment/outputs/novelty")
    out = tmp_path / "out"
    assert cc.main(["--folder", str(folder), "--out-dir", str(out), "--vocab", "pfam"]) == 0
    assert "yc_o3/novelty" in (out / "summary.md").read_text()
