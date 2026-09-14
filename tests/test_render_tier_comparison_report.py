import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import render_tier_comparison_report as rtcr  # noqa: E402


def test_guard_refuses_hand_written_content_without_force(tmp_path):
    out = tmp_path / 'report.md'
    out.write_text(
        '## Analysis 4: decision framework\n\n'
        '**Bottom line (revised after review):** Tier C is not recommended...\n'
    )
    error = rtcr.check_overwrite_guard(str(out), force=False)
    assert error is not None
    assert 'refusing to overwrite' in error
    assert '--force' in error
    # File must be left untouched.
    assert 'Bottom line' in out.read_text()


def test_guard_allows_placeholder_content_without_force(tmp_path):
    out = tmp_path / 'report.md'
    out.write_text(
        '## Analysis 4: decision framework\n\n'
        '_Fill in by hand after reading the tables above — see the spec...\n'
    )
    error = rtcr.check_overwrite_guard(str(out), force=False)
    assert error is None


def test_guard_allows_missing_file_without_force(tmp_path):
    out = tmp_path / 'does_not_exist.md'
    error = rtcr.check_overwrite_guard(str(out), force=False)
    assert error is None


def test_guard_allows_hand_written_content_with_force(tmp_path):
    out = tmp_path / 'report.md'
    out.write_text(
        '## Analysis 4: decision framework\n\n'
        '**Bottom line (revised after review):** Tier C is not recommended...\n'
    )
    error = rtcr.check_overwrite_guard(str(out), force=True)
    assert error is None
