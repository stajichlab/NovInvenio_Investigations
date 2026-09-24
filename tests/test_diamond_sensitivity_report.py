import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import diamond_sensitivity_report as dsr  # noqa: E402


def test_compare_counts():
    default = {'A::1', 'A::2', 'A::3'}
    mode = {'A::1', 'A::4'}
    ch = {'A::1', 'A::4', 'A::5'}
    r = dsr.compare(mode, default, ch)
    assert r['shared_with_default'] == 1
    assert r['gained_vs_default'] == 1
    assert r['lost_vs_default'] == 2
    # C+H extras vs default: A::4, A::5 -> mode recovers A::4 -> 0.5
    assert r['n_ch_extras_vs_default'] == 2
    assert r['ch_extras_recovered'] == 0.5
    # C+H misses vs default: A::2, A::3 -> both gone from mode -> 1.0
    assert r['ch_misses_resolved'] == 1.0
    assert r['ch_precision_vs_mode'] == round(2 / 3, 4)
    assert r['ch_recall_vs_mode'] == 1.0


def test_compare_empty_sets():
    r = dsr.compare(set(), set(), set())
    assert r['jaccard_vs_default'] == 1.0
    assert r['ch_extras_recovered'] == ''


def test_bins():
    assert dsr.e_bin(1e-60) == 'E<1e-50'
    assert dsr.e_bin(1e-30) == '1e-50<=E<1e-20'
    assert dsr.e_bin(1e-8) == '1e-20<=E<1e-5'
    assert dsr.cov_bin(10) == 'qcov<30'
    assert dsr.cov_bin(45) == '30<=qcov<60'
    assert dsr.cov_bin(100) == 'qcov>=60'


def test_best_other_hits(tmp_path):
    rows = [
        'p1\tt1\t1e-10\t50\t100\t40\t80\t70\t120\t130\n',
        'p1\tt2\t1e-30\t90\t100\t50\t25\t60\t120\t140\n',
        'p1\tt3\t1e-3\t20\t50\t30\t10\t10\t120\t400\n',
        'p2\tt4\t1e-40\t99\t100\t60\t90\t90\t100\t100\n',
    ]
    with gzip.open(tmp_path / 'A_vs_X.diamond.tsv.gz', 'wt') as fh:
        fh.writelines(rows[:3])
    with gzip.open(tmp_path / 'A_vs_Y.diamond.tsv.gz', 'wt') as fh:
        fh.writelines(rows[3:])
    best = dsr.best_other_hits(tmp_path, {'A': {'p1'}}, {'X', 'Y'}, 1e-5)
    assert best == {'p1': (1e-30, 25.0, 'X')}
