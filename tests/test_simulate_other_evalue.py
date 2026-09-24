import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import simulate_other_evalue as soe  # noqa: E402

QUERY = {'A', 'B', 'C', 'D'}
OTHER = {'X', 'Y'}


def _row(pres, evs):
    return {**{s: str(v) for s, v in pres.items()}}, {s: e for s, e in evs.items()}


def test_other_group_threshold_only():
    # present in A,B,C (0.75 query) and in X with E=1e-8
    pres = {'A': 1, 'B': 1, 'C': 1, 'D': 0, 'X': 1, 'Y': 0}
    evs = {'A': '', 'B': '1e-30', 'C': '1e-6', 'X': '1e-8'}
    assert not soe.is_candidate(pres, evs, QUERY, OTHER, 0.75, 0.0, 1e-5)
    assert soe.is_candidate(pres, evs, QUERY, OTHER, 0.75, 0.0, 1e-10)
    # query-group cells are never re-thresholded: C's weak 1e-6 hit still counts
    pres2 = dict(pres, X=0)
    assert soe.is_candidate(pres2, evs, QUERY, OTHER, 0.75, 0.0, 1e-30)


def test_present_cell_without_evalue_stays_present():
    pres = {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'X': 1, 'Y': 0}
    assert not soe.is_candidate(pres, {}, QUERY, OTHER, 0.75, 0.0, 1e-30)
