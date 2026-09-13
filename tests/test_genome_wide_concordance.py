import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import genome_wide_concordance as gwc  # noqa: E402


def test_load_candidates(tmp_path):
    p = tmp_path / 'candidates.txt'
    p.write_text('Ncra::pA1\nAfum::pB1\n')
    assert gwc.load_candidates(p) == {'Ncra::pA1', 'Afum::pB1'}


def test_load_candidates_empty_file(tmp_path):
    p = tmp_path / 'candidates.txt'
    p.write_text('')
    assert gwc.load_candidates(p) == set()


def test_jaccard():
    assert gwc.jaccard({'a', 'b'}, {'b', 'c'}) == 1 / 3
    assert gwc.jaccard(set(), set()) == 1.0  # both empty -> perfect agreement, not 0/0


def test_precision_recall():
    predicted = {'a', 'b', 'c'}
    gold = {'a', 'b', 'd'}
    precision, recall = gwc.precision_recall(predicted, gold)
    assert precision == 2 / 3
    assert recall == 2 / 3


def test_count_inflation():
    baseline = {'a', 'b'}
    refined = {'a', 'b', 'c', 'd', 'e'}
    result = gwc.count_inflation(baseline, refined)
    assert result == {'baseline_count': 2, 'refined_count': 5, 'delta': 3, 'delta_pct': 150.0}
