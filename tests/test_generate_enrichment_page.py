import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

import generate_enrichment_page as gep  # noqa: E402
from study_runs import inject_run_report_links, rebuild_study_pages  # noqa: E402

REPORT = """<html><body>
    <div class="cards">
      <a class="card-link" href="novelties.html"><span class="card-link-title">N</span></a>
      <a class="card-link" href="core.html"><span class="card-link-title">C</span></a>
      <a class="card-link" href="losses.html"><span class="card-link-title">L</span></a>
    </div>
</body></html>
"""

DOMAIN_HEADER = ("scope\tlevel\tvocab\tterm\tname\tk_cand_with_term\tN_cand_annotated\tn_bg_with_term\t"
                 "M_bg_annotated\texpected\tfold\tpvalue\tfdr_bh\tfdr_bh_all_species\tcandidate_ids\n")
GO_HEADER = DOMAIN_HEADER.rstrip("\n").replace("\tcandidate_ids", "") + (
    "\tgo_evidence\tn_cand_units_non_IEA\tcand_evidence_codes\tcandidate_ids\n")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _outputs(d: Path, root: Path, folder: str) -> None:
    d.mkdir(parents=True)
    (d / "coverage.tsv").write_text(
        "species\tsubset\tn_proteins\tn_in_matrix\tn_pfam\tfrac_pfam\tn_pfam_clan\tfrac_pfam_clan\t"
        "n_interpro\tfrac_interpro\tn_go\tfrac_go\tn_any\tfrac_any\n"
        "ALL_IN\tuniverse\t100\t90\t80\t0.8\t70\t0.7\t85\t0.85\t60\t0.6\t88\t0.88\n"
        "ALL_IN\tcandidates\t10\t10\t4\t0.4\t3\t0.3\t6\t0.6\t2\t0.2\t6\t0.6\n")
    (d / "tests.tsv").write_text(
        "study\tuniverse\tbackground_clustered\tscope\tlevel\tvocab\trole\tn_cand_units\tN_cand_annotated\t"
        "frac_cand_annotated\tn_bg_units\tM_units_annotated\tn_terms_tested\tn_fdr_lt_0.05\t"
        "n_fdr_lt_0.05_k_ge_3\tstatus\n"
        "s\tproteome\tTrue\twhole_set\tfamily\tpfam\tprimary\t10\t4\t0.4\t90\t80\t12\t1\t1\tok\n")
    with gzip.open(d / "enrichment_domains.tsv.gz", "wt") as fh:
        fh.write(DOMAIN_HEADER)
        fh.write("whole_set\tfamily\tpfam\tPF00001\tX: thing\t3\t4\t5\t80\t0.2\t12\t1e-5\t0.01\t\tA::p1;A::p2\n")
        fh.write("whole_set\tfamily\tinterpro\tIPR000001\tDomain:Y\t3\t6\t5\t85\t0.2\t8\t1e-4\t0.2\t\tA::p1\n")
    with gzip.open(d / "enrichment_go.tsv.gz", "wt") as fh:
        fh.write(GO_HEADER)
        fh.write("whole_set\tfamily\tgo\tGO:0000001\tbp:z\t2\t2\t5\t60\t0.1\t9\t1e-3\t0.04\t\tIEA_only\t0\tIEA:2\t"
                 "A::secret_protein\n")
    (d / "controls.tsv").write_text("control\tpfam_ids\tcandidates\tcandidate_ids\nArp2\tPF00022\t1\tA::secret_ctl\n")
    (d / "PROVENANCE.yaml").write_text(
        f"generated: '2026-10-12'\nuniverse: proteome\ncommand: python x.py --results {root}/results/{folder} "
        "--study s\nparameters:\n  cluster_background: true\n")


@pytest.fixture
def repo(tmp_path):
    folder = tmp_path / "studies" / "fungi" / "s1"
    folder.mkdir(parents=True)
    (folder / "publish.yaml").write_text("study: s1\nrun: pw\n")
    _outputs(folder / "analysis" / "enrichment" / "outputs" / "proteome_clustered", tmp_path, "s1")
    run = tmp_path / "docs" / "fungi" / "s1" / "pw"
    run.mkdir(parents=True)
    (run / "report.html").write_text(REPORT)
    (run / "run.json").write_text(json.dumps({"run": "pw", "source_dir": "results/s1",
                                              "report_sha256": _sha(run / "report.html")}))
    return tmp_path


def _main(repo, *args):
    return gep.main(["--repo-root", str(repo), "--study", "fungi/s1", *args])


def test_novelty_page_links_and_sha(repo):
    run = repo / "docs" / "fungi" / "s1" / "pw"
    assert _main(repo) == 0
    page = (run / "enrichment.html").read_text()
    assert "PF00001" in page and "https://www.ebi.ac.uk/interpro/entry/pfam/PF00001/" in page
    assert "GO:0000001" in page and "Low annotation" in page
    assert "PF00022" in page                                   # controls keep pfam_ids
    assert "secret_protein" not in page and "secret_ctl" not in page  # no protein IDs
    assert (run / "enrichment_provenance.yaml").is_file()
    report = (run / "report.html").read_text()
    assert report.count('href="enrichment.html"') == 1
    assert report.index('href="losses.html"') < report.index('href="enrichment.html"')
    assert json.loads((run / "run.json").read_text())["report_sha256"] == _sha(run / "report.html")
    assert not (run / "enrichment_loss.html").exists()
    before = (run / "report.html").read_bytes()
    assert _main(repo) == 0                                    # idempotent
    assert (run / "report.html").read_bytes() == before


COV_HEAD = ("species\tsubset\tn_proteins\tn_in_matrix\tn_pfam\tfrac_pfam\tn_pfam_clan\tfrac_pfam_clan\t"
            "n_interpro\tfrac_interpro\tn_go\tfrac_go\tn_any\tfrac_any\n")


def test_no_vocabulary_gap_no_caveat(repo):
    assert _main(repo) == 0
    page = (repo / "docs" / "fungi" / "s1" / "pw" / "enrichment.html").read_text()
    assert "cover only" not in page and "covers only" not in page


def test_vocabulary_gap_caveat(repo):
    d = repo / "studies" / "fungi" / "s1" / "analysis" / "enrichment" / "outputs" / "proteome_clustered"
    (d / "coverage.tsv").write_text(
        COV_HEAD
        + "Uni\tuniverse\t100\t90\t80\t0.8\t70\t0.7\t85\t0.85\t60\t0.6\t88\t0.88\n"
        + "Uni\tcandidates\t10\t10\t4\t0.4\t3\t0.3\t6\t0.6\t2\t0.2\t6\t0.6\n"
        + "HmmA\tuniverse\t50\t40\t30\t0.6\t20\t0.4\t0\t0.0\t0\t0.0\t30\t0.6\n"
        + "HmmB\tuniverse\t50\t40\t30\t0.6\t20\t0.4\t0\t0.0\t0\t0.0\t30\t0.6\n"
        + "None1\tuniverse\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"           # no proteins: ignored
        + "ALL_IN\tuniverse\t200\t170\t140\t0.7\t110\t0.55\t85\t0.425\t60\t0.3\t148\t0.74\n"
        + "ALL_IN\tcandidates\t10\t10\t4\t0.4\t3\t0.3\t6\t0.6\t2\t0.2\t6\t0.6\n")
    assert _main(repo) == 0
    page = (repo / "docs" / "fungi" / "s1" / "pw" / "enrichment.html").read_text()
    assert "GO and InterPro cover only 1 of 3 species: Uni." in page
    assert "No GO and InterPro annotation for HmmA, HmmB" in page
    assert "Pfam covers only" not in page and "Pfam and" not in page


def test_pfam_only_study_not_available(repo):
    """Every species has Pfam-only HMM annotation: GO and InterPro say 'not available'."""
    d = repo / "studies" / "fungi" / "s1" / "analysis" / "enrichment" / "outputs" / "proteome_clustered"
    (d / "coverage.tsv").write_text(
        COV_HEAD
        + "HmmA\tuniverse\t50\t40\t30\t0.6\t20\t0.4\t0\t0.0\t0\t0.0\t30\t0.6\n"
        + "HmmB\tuniverse\t50\t40\t30\t0.6\t20\t0.4\t0\t0.0\t0\t0.0\t30\t0.6\n"
        + "ALL_IN\tuniverse\t100\t80\t60\t0.6\t40\t0.4\t0\t0.0\t0\t0.0\t60\t0.6\n"
        + "ALL_IN\tcandidates\t10\t10\t6\t0.6\t3\t0.3\t0\t0.0\t0\t0.0\t6\t0.6\n")
    with gzip.open(d / "enrichment_domains.tsv.gz", "wt") as fh:
        fh.write(DOMAIN_HEADER)
        fh.write("whole_set\tfamily\tpfam\tPF00001\tX: thing\t3\t4\t5\t80\t0.2\t12\t1e-5\t0.01\t\tA::p1\n")
    with gzip.open(d / "enrichment_go.tsv.gz", "wt") as fh:
        fh.write(GO_HEADER)
    (d / "PROVENANCE.yaml").write_text(
        (d / "PROVENANCE.yaml").read_text()
        + "annotation_sources:\n  HmmA: 'Pfam-A HMM search (hmmsearch); Pfam only'\n"
          "  HmmB: 'Pfam-A HMM search (hmmsearch); Pfam only'\n")
    assert _main(repo) == 0
    page = (repo / "docs" / "fungi" / "s1" / "pw" / "enrichment.html").read_text()
    assert "GO and InterPro: not available for this study." in page
    assert page.count("<b>Not available</b>") == 2                  # InterPro and GO top-hit sections
    assert "GO: not available for this study" in page
    assert "Pfam-A HMM search (hmmsearch, Pfam gathering thresholds) for every species" in page
    assert "not the UniProt release" not in page
    assert "PF00001" in page
    # 2026-10-16: vocabulary tag, link text, no repeated 0% low-annotation line.
    assert '<meta name="nii-enrichment-vocabularies" content="pfam">' in page
    assert "Low annotation" not in page                      # only InterPro/GO 0% were low
    report = (repo / "docs" / "fungi" / "s1" / "pw" / "report.html").read_text()
    assert "Pfam term enrichment of the novelty candidates" in report
    assert "Pfam, InterPro and GO" not in report


def test_pfam_only_low_pfam_keeps_low_line_without_go_interpro(repo):
    d = repo / "studies" / "fungi" / "s1" / "analysis" / "enrichment" / "outputs" / "proteome_clustered"
    (d / "coverage.tsv").write_text(
        COV_HEAD
        + "HmmA\tuniverse\t50\t40\t30\t0.6\t20\t0.4\t0\t0.0\t0\t0.0\t30\t0.6\n"
        + "ALL_IN\tuniverse\t50\t40\t30\t0.6\t20\t0.4\t0\t0.0\t0\t0.0\t30\t0.6\n"
        + "ALL_IN\tcandidates\t10\t10\t2\t0.2\t1\t0.1\t0\t0.0\t0\t0.0\t2\t0.2\n")
    assert _main(repo) == 0
    page = (repo / "docs" / "fungi" / "s1" / "pw" / "enrichment.html").read_text()
    assert ("Low annotation: the candidate set has less than 50% annotation for Pfam 20.0%, "
            "any vocabulary 20.0%.") in page


def test_full_annotation_page_tag_and_default_link_text(repo):
    assert _main(repo) == 0
    run = repo / "docs" / "fungi" / "s1" / "pw"
    page = (run / "enrichment.html").read_text()
    assert '<meta name="nii-enrichment-vocabularies" content="pfam,interpro,go">' in page
    assert "Low annotation: the candidate set has less than 50% annotation for Pfam 40.0%, GO 20.0%." in page
    assert ("Pfam, InterPro and GO term enrichment of the novelty candidates against the proteome "
            "background.") in (run / "report.html").read_text()


def test_vocabulary_gaps_groups():
    def r(sp, n, pf, cl, ipr, go):
        return {"species": sp, "subset": "universe", "n_proteins": n, "frac_pfam": pf,
                "frac_pfam_clan": cl, "frac_interpro": ipr, "frac_go": go, "frac_any": "0"}
    lines = gep.vocabulary_gaps([r("A", "10", "0.5", "0.4", "0.5", "0.5"), r("B", "10", "0", "0", "0", "0")])
    assert len(lines) == 1 and lines[0].startswith("GO, InterPro, Pfam and Pfam clan cover only 1 of 2 species: A.")
    lines = gep.vocabulary_gaps([r("A", "10", "0.5", "0.4", "0.5", "0"), r("B", "10", "0.5", "0.4", "0", "0.2")])
    assert lines[0].startswith("GO covers only 1 of 2 species: B.")
    assert lines[1].startswith("InterPro covers only 1 of 2 species: A.")


def test_loss_absent_writes_nothing(repo):
    assert _main(repo, "--direction", "loss") == 0
    run = repo / "docs" / "fungi" / "s1" / "pw"
    assert not (run / "enrichment_loss.html").exists()
    assert (run / "report.html").read_text() == REPORT


def test_loss_page_and_both_links(repo):
    enr = repo / "studies" / "fungi" / "s1" / "analysis" / "enrichment" / "outputs"
    _outputs(enr / "loss" / "proteome_clustered", repo, "s1")
    assert _main(repo) == 0
    assert _main(repo, "--direction", "loss") == 0
    run = repo / "docs" / "fungi" / "s1" / "pw"
    report = (run / "report.html").read_text()
    assert 'href="enrichment.html"' in report and 'href="enrichment_loss.html"' in report
    assert (run / "enrichment_loss_provenance.yaml").is_file()


def test_links_only_restores_links_after_resync(repo):
    run = repo / "docs" / "fungi" / "s1" / "pw"
    _main(repo)
    (run / "report.html").write_text(REPORT)                   # sync_reports.sh copies a fresh report
    assert gep.main(["--links-only", "--docs-run-dir", str(run)]) == 0
    assert 'href="enrichment.html"' in (run / "report.html").read_text()


def test_missing_column_stops(repo):
    d = repo / "studies" / "fungi" / "s1" / "analysis" / "enrichment" / "outputs" / "proteome_clustered"
    with gzip.open(d / "enrichment_domains.tsv.gz", "wt") as fh:
        fh.write("scope\tlevel\tvocab\tterm\n")
    with pytest.raises(SystemExit, match="fdr_bh"):
        _main(repo)


def test_source_mismatch_stops(repo):
    run = repo / "docs" / "fungi" / "s1" / "pw"
    (run / "run.json").write_text(json.dumps({"run": "pw", "source_dir": "results/other"}))
    with pytest.raises(SystemExit, match="source_dir"):
        _main(repo)


def test_run_list_links_only_existing_pages(tmp_path):
    for run in ("a", "b"):
        (tmp_path / run).mkdir()
        (tmp_path / run / "run.json").write_text(json.dumps({"run": run, "published": "2026-10-01"}))
    (tmp_path / "a" / "enrichment.html").write_text("x")
    rebuild_study_pages(tmp_path, "fungi", "s")
    html = (tmp_path / "report.html").read_text()
    assert 'href="a/enrichment.html"' in html and "b/enrichment" not in html


def test_inject_without_pages_leaves_report_unchanged(tmp_path):
    (tmp_path / "report.html").write_text(REPORT)
    inject_run_report_links(tmp_path)
    assert (tmp_path / "report.html").read_text() == REPORT
