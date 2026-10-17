"""Tests for bin/pangenome_gain_loss_enrichment.py on a tiny synthetic table."""
import csv
import gzip
import importlib.util
import math
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "pangenome_gain_loss_enrichment.py"
spec = importlib.util.spec_from_file_location("pangenome_gain_loss_enrichment", SCRIPT)
ge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ge)

PFAM_DAT = """# STOCKHOLM 1.0
#=GF ID   DomA
#=GF AC   PF00001.5
#=GF DE   Domain A
#=GF CL   CL0001
//
#=GF ID   DomB
#=GF AC   PF00002.3
#=GF DE   Domain B
//
#=GF ID   DomC
#=GF AC   PF00003.1
#=GF DE   Domain C
#=GF CL   CL0001
//
"""


def exact_sf(k, M, K, N):
    """P(X >= k) by direct summation (independent of scipy)."""
    tot = math.comb(M, N)
    return sum(math.comb(K, i) * math.comb(M - K, N - i) for i in range(k, min(K, N) + 1)) / tot


def _write_gz(path, header, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


@pytest.fixture
def synth(tmp_path):
    # 20 shell/soft_core families: f00-f07 gain, f08-f13 loss, f14-f17 ambiguous,
    # f18 shell without Pfam, f19 cloud (out of the default scope), f20 core.
    gl, ann = [], []
    for i in range(21):
        fam = f"f{i:02d}"
        cls = "cloud" if i == 19 else "core" if i == 20 else ("soft_core" if i % 5 == 0 else "shell")
        d = "gain" if i < 8 else "loss" if i < 14 else "ambiguous"
        if i == 19:
            d = "gain"
        gl.append([fam, cls, d, "conserved" if cls == "core" else d])
        # PF00001 on all 8 gain families + f08; PF00002 on loss f09-f13 + f14,f15; PF00003 on f16,f17
        if i < 9:
            p = "PF00001.5"
        elif i < 16:
            p = "PF00002.3"
        elif i < 18:
            p = "PF00003.1"
        elif i == 19:
            p = "PF00001.5"
        else:
            p = ""
        ann.append([fam, p, "DomA" if p.startswith("PF00001") else "DomB" if p.startswith("PF00002")
                    else "DomC" if p else ""])
    _write_gz(tmp_path / "gl.tsv.gz", ["family", "class", "rule_freq", "direction"], gl)
    _write_gz(tmp_path / "ann.tsv.gz", ["accession", "pfam_ids", "pfam_names"], ann)
    (tmp_path / "Pfam-A.hmm.dat").write_text(PFAM_DAT)
    return tmp_path


def _read(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_run_counts_and_pvalues(synth):
    out = synth / "out"
    ge.main(["--gain-loss-table", str(synth / "gl.tsv.gz"), "--family-annotation", str(synth / "ann.tsv.gz"),
             "--out-dir", str(out), "--pfam-dat", str(synth / "Pfam-A.hmm.dat"), "--label", "t",
             "--min-universe", "2"])
    rows = _read(out / "enrichment_domains.tsv.gz")
    a = {(r["test_set"], r["vocab"], r["term"]): r for r in rows}
    # scope = shell+soft_core = f00..f18 (19 families); with Pfam = 18 (f18 empty) -> M = 18
    r = a[("gain", "pfam", "PF00001")]
    assert (int(r["k"]), int(r["n_set"]), int(r["K"]), int(r["M"])) == (8, 8, 9, 18)
    assert float(r["p"]) == pytest.approx(exact_sf(8, 18, 9, 8), rel=1e-9)
    assert r["name"] == "DomA: Domain A"
    r = a[("loss", "pfam", "PF00002")]
    assert (int(r["k"]), int(r["n_set"]), int(r["K"]), int(r["M"])) == (5, 6, 7, 18)
    assert float(r["p"]) == pytest.approx(exact_sf(5, 18, 7, 6), rel=1e-9)
    # k = 0 -> p = 1
    assert float(a[("gain", "pfam", "PF00002")]["p"]) == 1.0
    # clan CL0001 = PF00001 + PF00003; universe for clans = families with a clan (f00-f08, f16, f17)
    r = a[("gain", "pfam_clan", "CL0001")]
    assert (int(r["k"]), int(r["K"]), int(r["M"])) == (8, 11, 11)
    assert float(r["p"]) == pytest.approx(1.0)
    tests = _read(out / "tests.tsv")
    assert {(t["test_set"], t["vocab"]) for t in tests} == {("gain", "pfam"), ("gain", "pfam_clan"),
                                                          ("loss", "pfam"), ("loss", "pfam_clan")}
    cov = {c["direction"]: c for c in _read(out / "coverage.tsv")}
    assert int(cov["ALL"]["n_families"]) == 19 and int(cov["ALL"]["n_with_pfam"]) == 18
    assert int(cov["ambiguous"]["n_families"]) == 5   # f14-f18
    assert (out / "PROVENANCE.yaml").exists()


def test_min_universe_filters_terms(synth):
    out = synth / "out5"
    ge.main(["--gain-loss-table", str(synth / "gl.tsv.gz"), "--family-annotation", str(synth / "ann.tsv.gz"),
             "--out-dir", str(out), "--pfam-dat", str(synth / "Pfam-A.hmm.dat"), "--min-universe", "5"])
    terms = {r["term"] for r in _read(out / "enrichment_domains.tsv.gz") if r["vocab"] == "pfam"}
    assert terms == {"PF00001", "PF00002"}       # PF00003 has K = 2 < 5


def test_classes_option_changes_scope(synth):
    out = synth / "outall"
    ge.main(["--gain-loss-table", str(synth / "gl.tsv.gz"), "--family-annotation", str(synth / "ann.tsv.gz"),
             "--out-dir", str(out), "--pfam-dat", str(synth / "Pfam-A.hmm.dat"), "--min-universe", "2",
             "--classes", "shell,soft_core,cloud"])
    r = next(r for r in _read(out / "enrichment_domains.tsv.gz")
             if r["test_set"] == "gain" and r["term"] == "PF00001")
    assert (int(r["k"]), int(r["n_set"]), int(r["K"]), int(r["M"])) == (9, 9, 10, 19)


def test_deterministic_output(synth):
    outs = []
    for name in ("d1", "d2"):
        ge.main(["--gain-loss-table", str(synth / "gl.tsv.gz"), "--family-annotation",
                 str(synth / "ann.tsv.gz"), "--out-dir", str(synth / name), "--pfam-dat",
                 str(synth / "Pfam-A.hmm.dat"), "--min-universe", "2"])
        outs.append((synth / name / "enrichment_domains.tsv.gz").read_bytes())
    assert outs[0] == outs[1]


def test_missing_family_in_annotation_stops(synth, tmp_path):
    _write_gz(tmp_path / "ann_short.tsv.gz", ["accession", "pfam_ids", "pfam_names"], [["f00", "", ""]])
    with pytest.raises(SystemExit):
        ge.main(["--gain-loss-table", str(synth / "gl.tsv.gz"), "--family-annotation",
                 str(tmp_path / "ann_short.tsv.gz"), "--out-dir", str(tmp_path / "x"),
                 "--pfam-dat", str(synth / "Pfam-A.hmm.dat")])


def test_bh_and_enrich_function():
    uni = {f"f{i}": {"A"} if i < 6 else {"B"} for i in range(12)}
    rows = ge.enrich(uni, {"f0", "f1", "f2", "f3"}, min_universe=5, report_ids_fdr=1.0)
    a = {r["term"]: r for r in rows}
    assert a["A"]["k"] == 4 and a["A"]["K"] == 6 and a["A"]["M"] == 12 and a["A"]["n_set"] == 4
    assert a["A"]["p"] == pytest.approx(exact_sf(4, 12, 6, 4))
    assert a["B"]["p"] == 1.0
    # BH with 2 tests: q_A = min(p_A * 2 / 1, q_B)
    assert a["A"]["q"] == pytest.approx(min(a["A"]["p"] * 2, 1.0))
    assert a["A"]["family_ids"] == "f0;f1;f2;f3"


def test_go_terms_strip_evidence():
    assert ge.split_terms("GO:0000001:IEA|GO:0000002", "go") == {"GO:0000001", "GO:0000002"}
    assert ge.split_terms("PF00001.5|PF00002", "pfam") == {"PF00001", "PF00002"}
