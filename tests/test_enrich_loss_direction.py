"""Loss-direction candidate rule used by the study enrichment and loss near-miss scripts.

Scripts under test (study-specific, loaded by absolute path):
  studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py   rule_candidates()
  studies/eukaryotes/zoosporic_holozoa_vs_dikarya/analysis/loss_near_miss/loss_loose.py
      filter_pair_hits(), passes()
The rule mirrors nf_NovInvenio bin/build_presence_matrix.py --query-group OUT:
query (target) count / n_query >= min_frac and other count / n_other <= other_max_frac.
"""
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("scipy")

REPO = Path(__file__).resolve().parent.parent
ENRICH = REPO / "studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py"
LOOSE = REPO / "studies/eukaryotes/zoosporic_holozoa_vs_dikarya/analysis/loss_near_miss/loss_loose.py"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


enrich = _load(ENRICH, "enrich_novelties_v2")
loose = _load(LOOSE, "loss_loose")

# Tiny loss-direction matrix: 4 OUT (target) proteomes, 3 IN (other) proteomes.
# Rows are OUT proteins; target count includes the protein's own proteome.
OUT = ["O1", "O2", "O3", "O4"]
IN = ["I1", "I2", "I3"]
ROWS = {
    "O1::a": ({"O1", "O2", "O3", "O4"}, set()),          # 4/4 OUT, 0 IN -> loss candidate
    "O1::b": ({"O1", "O2", "O3"}, set()),                # 3/4 = 0.75 OUT, 0 IN -> candidate (>=)
    "O2::c": ({"O2", "O3"}, set()),                      # 2/4 OUT -> fails OUT rule
    "O3::d": ({"O1", "O2", "O3", "O4"}, {"I2"}),         # 1 IN hit -> not a strict candidate
    "O4::e": ({"O1", "O2", "O3", "O4"}, {"I1", "I3"}),   # 2 IN hits
}


def _counts():
    return {k: (len(o), len(i)) for k, (o, i) in ROWS.items()}


def test_rule_candidates_loss_strict():
    rule, cands = enrich.rule_candidates(_counts(), len(OUT), len(IN), 0.75, 0.0)
    assert rule == {"O1::a", "O1::b", "O3::d", "O4::e"}
    assert cands == {"O1::a", "O1::b"}


def test_rule_candidates_other_max_frac_relaxed():
    # 1/3 IN <= 0.34 keeps O3::d; 2/3 IN does not.
    _, cands = enrich.rule_candidates(_counts(), len(OUT), len(IN), 0.75, 0.34)
    assert cands == {"O1::a", "O1::b", "O3::d"}


def test_novelty_default_matches_zero_other_hits():
    # Novelty default (other_max_frac 0.0) is the old "other count == 0" test.
    counts = {"x": (3, 0), "y": (3, 1), "z": (2, 0)}
    _, cands = enrich.rule_candidates(counts, 3, 5, 0.75)
    assert cands == {"x"}


@pytest.mark.parametrize("k,expected", [
    (0, {"O1::a", "O1::b"}),
    (1, {"O1::a", "O1::b", "O3::d"}),
    (2, {"O1::a", "O1::b", "O3::d", "O4::e"}),
])
def test_loose_k(k, expected):
    got = {key for key, (o, i) in ROWS.items() if loose.passes(len(o), len(OUT), len(i), 0.75, k)}
    assert got == expected
    # K = 0 agrees with the strict enrichment rule
    if k == 0:
        assert got == enrich.rule_candidates(_counts(), len(OUT), len(IN), 0.75, 0.0)[1]


def test_filter_pair_hits_pipeline_filters():
    paralog_of = {"q1": "q1p", "q2": "q2p"}
    rows = [
        ("q1", "t1", 1e-10),   # paralog q1p hits t1 better (1e-12): dropped (no rescue, 1e-10 > 1e-20)
        ("q1p", "t1", 1e-12),
        ("q2", "t2", 1e-30),   # paralog better (1e-40) but rescued (1e-30 <= 1e-20)
        ("q2p", "t2", 1e-40),
        ("q3", "t3", 1e-5),    # filter 1 is strict: 1e-5 is not < 1e-5 -> dropped
        ("q4", "t4", 1e-6),    # kept
        ("q4", "t5", 1e-8),    # kept, better
        ("q5", "t6", 0.5),     # above parse cutoff
    ]
    kept = loose.filter_pair_hits(rows, paralog_of)
    assert set(kept) == {"q1p", "q2", "q2p", "q4"}
    assert kept["q4"] == (1e-8, "t5")


def test_filter_pair_hits_paralog_on_other_target_does_not_drop():
    # Scope 'target': the paralog beating the query on a different target protein is ignored.
    rows = [("q1", "t1", 1e-10), ("q1p", "t9", 1e-50)]
    assert set(loose.filter_pair_hits(rows, {"q1": "q1p"})) == {"q1", "q1p"}
