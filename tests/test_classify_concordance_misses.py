import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import classify_concordance_misses as ccm  # noqa: E402

QUERY = {'A', 'B', 'C', 'D'}
OTHER = {'X', 'Y'}


def _pres(query_present, other_present):
    p = {s: int(s in query_present) for s in QUERY}
    p.update({s: int(s in other_present) for s in OTHER})
    return p


def _setup():
    member_to_rep = {'a1': 'a1', 'b1': 'a1', 'c1': 'a1', 'a2': 'a1',
                     'solo': 'solo', 'big': 'big',
                     'f_a': 'f_a', 'f_b': 'f_a'}
    rep_members = {'a1': ['a1', 'b1', 'c1', 'a2'], 'solo': ['solo'], 'big': ['big'],
                   'f_a': ['f_a', 'f_b']}
    protein_to_proteome = {'a1': 'A', 'b1': 'B', 'c1': 'C', 'a2': 'A', 'solo': 'A',
                           'big': 'A', 'f_a': 'A', 'f_b': 'B'}
    return member_to_rep, rep_members, protein_to_proteome


def _miss(pid, ch_matrix, paralog_of=None, profiled=('a1', 'f_a', 'big'), oversized=('big',)):
    m2r, rm, p2p = _setup()
    return ccm.classify_miss(pid, ch_matrix, m2r, rm, set(profiled), set(oversized), p2p,
                             paralog_of or {}, QUERY, OTHER, 0.75, 0.0)


def test_unprofiled_and_unclustered():
    assert _miss('solo', {}) == 'unprofiled_family'
    assert _miss('not_clustered', {}) == 'unprofiled_family'


def test_oversized():
    assert _miss('big', {}) == 'oversized_family'


def test_hmm_outgroup_paralog_vs_other():
    ch = {'a1': ('A', _pres({'A', 'B', 'C'}, {'X'}))}
    assert _miss('a1', ch, paralog_of={'a1': 'a2'}) == 'hmm_outgroup_paralog'
    assert _miss('a1', ch, paralog_of={'a1': 'elsewhere'}) == 'hmm_outgroup_other'
    assert _miss('a1', ch) == 'hmm_outgroup_other'


def test_fragmentation_vs_hmm_undercall():
    # f_a family: membership covers A,B only (0.5 < 0.75) -> fragmentation.
    assert _miss('f_a', {'f_a': ('A', _pres({'A', 'B'}, set()))}) == 'cluster_fragmentation'
    # a1 family: membership covers A,B,C (0.75), HMM only calls A,B -> HMM undercall.
    assert _miss('a1', {'a1': ('A', _pres({'A', 'B'}, set()))}) == 'hmm_query_undercall'


def test_unexplained_when_matrix_passes():
    assert _miss('a1', {'a1': ('A', _pres({'A', 'B', 'C'}, set()))}) == 'unexplained'


def test_member_row_falls_back_to_rep_row():
    ch = {'a1': ('A', _pres({'A', 'B', 'C'}, {'Y'}))}
    assert _miss('b1', ch) == 'hmm_outgroup_other'


def test_classify_extra():
    p = {'x': ('A', _pres({'A', 'B', 'C'}, {'X'})),
         'y': ('A', _pres({'A'}, set())),
         'z': ('A', _pres({'A', 'B', 'C'}, set()))}
    assert ccm.classify_extra('x', p, QUERY, OTHER, 0.75, 0.0) == 'p_other_presence'
    assert ccm.classify_extra('y', p, QUERY, OTHER, 0.75, 0.0) == 'p_query_short'
    assert ccm.classify_extra('w', p, QUERY, OTHER, 0.75, 0.0) == 'p_no_hits'
    assert ccm.classify_extra('z', p, QUERY, OTHER, 0.75, 0.0) == 'unexplained'


def test_level_of():
    assert ccm.level_of('unprofiled_family') == 'clustering'
    assert ccm.level_of('hmm_outgroup_other') == 'hmm'
    assert ccm.level_of('p_no_hits') == 'pairwise'
    assert ccm.level_of('unexplained') == 'other'
