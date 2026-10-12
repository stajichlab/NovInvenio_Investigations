"""bin/go_enrichment.py: one-sided hypergeometric p and BH over all tested terms."""
import csv
import subprocess
import sys
from pathlib import Path

import pytest
from scipy.stats import hypergeom

pytest.importorskip("goatools")

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
from stats_util import bh_fdr  # noqa: E402

OBO = """format-version: 1.2
ontology: go

[Term]
id: GO:0000001
name: term_a
namespace: molecular_function

[Term]
id: GO:0000002
name: term_b
namespace: molecular_function

[Term]
id: GO:0000003
name: term_c
namespace: molecular_function
"""


def test_one_sided_p_and_bh_over_all_tested_terms(tmp_path):
    # 30 background proteins. term_a: 6 proteins (all candidates). term_b: 12 (2 candidates,
    # under-represented). term_c: 12 (2 candidates, under-represented).
    prots = [f"P{i:02d}" for i in range(30)]
    terms = {p: set() for p in prots}
    for p in prots[:6]:
        terms[p].add("GO:0000001")
    for p in prots[6:18]:
        terms[p].add("GO:0000002")
    for p in prots[18:30]:
        terms[p].add("GO:0000003")
    cands = prots[:6] + [prots[6], prots[7], prots[18], prots[19]]
    # candidates also carry the other terms for the first four
    terms[prots[6]].add("GO:0000002"); terms[prots[7]].add("GO:0000002")
    ann = tmp_path / "ann.tsv"
    with open(ann, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "go_ids"])
        for p in prots:
            w.writerow([p, "|".join(f"{t}:IDA" for t in sorted(terms[p]))])
    (tmp_path / "cand.txt").write_text("\n".join(cands) + "\n")
    (tmp_path / "bg.txt").write_text("\n".join(prots) + "\n")
    (tmp_path / "go.obo").write_text(OBO)
    out = tmp_path / "out.tsv"
    r = subprocess.run(
        [sys.executable, str(REPO / "bin" / "go_enrichment.py"),
         "--annotations", str(ann), "--obo", str(tmp_path / "go.obo"),
         "--candidates", str(tmp_path / "cand.txt"), "--background", str(tmp_path / "bg.txt"),
         "--output", str(out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    rows = {x["GO_id"]: x for x in csv.DictReader(open(out), delimiter="\t")}
    a = rows["GO:0000001"]
    n_cand = int(a["n_candidates_annotatable"]); n_bg = int(a["n_background_annotatable"])
    # one-sided over-representation p for term_a: 6 of 6 background carriers are candidates
    expected_p = hypergeom.sf(6 - 1, n_bg, 6, n_cand)
    assert float(a["pvalue"]) == pytest.approx(expected_p)
    # BH must use every tested term (3 here), not only the over-represented ones
    all_p = []
    for t, nbg in (("GO:0000001", 6), ("GO:0000002", 12), ("GO:0000003", 12)):
        k = {"GO:0000001": 6, "GO:0000002": 4, "GO:0000003": 2}[t]
        all_p.append(hypergeom.sf(k - 1, n_bg, nbg, n_cand))
    expected_q = bh_fdr(all_p)[0]
    assert float(a["fdr_bh"]) == pytest.approx(expected_q)
    assert "3 terms tested" in r.stderr
