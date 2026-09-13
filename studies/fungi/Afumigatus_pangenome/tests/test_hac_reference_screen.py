import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from hac_reference_screen import is_present


def test_is_present_true_ortholog_passes():
    hit = {"pident": 98.6, "align_len": 279}
    assert is_present(hit, qlen=362, min_pident=50.0, min_qcov=0.5) is True


def test_is_present_weak_paralog_cross_hit_fails_on_identity():
    # A different PF11001-family paralog, not the true hrmA ortholog --
    # real data showed a clean bimodal split at ~35% vs ~99-100% identity.
    hit = {"pident": 34.7, "align_len": 350}
    assert is_present(hit, qlen=362, min_pident=50.0, min_qcov=0.5) is False


def test_is_present_high_identity_but_low_coverage_fails():
    # High identity but only a small fragment aligned -- below the coverage
    # floor even though it would pass on identity alone.
    hit = {"pident": 100.0, "align_len": 52}
    assert is_present(hit, qlen=362, min_pident=50.0, min_qcov=0.5) is False
