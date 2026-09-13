import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import trace_cost_report as tcr  # noqa: E402


def test_parse_duration_hours_minutes_seconds():
    assert tcr.parse_duration('2h 8m 50s') == 2 * 3600 + 8 * 60 + 50


def test_parse_duration_minutes_seconds():
    assert tcr.parse_duration('7m 34s') == 7 * 60 + 34


def test_parse_duration_milliseconds():
    assert tcr.parse_duration('346ms') == 0.346


def test_parse_duration_zero():
    assert tcr.parse_duration('0') == 0.0


def test_load_trace(tmp_path):
    p = tmp_path / 'trace.txt'
    p.write_text(
        'task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n'
        '1\ta\tb\tDIAMOND_SEARCH (x)\tCOMPLETED\t0\t2026-01-01\t1h\t45m\t100%\t1\t1\t1\t1\n'
        '2\ta\tb\tHMMSEARCH (y)\tFAILED\t1\t2026-01-01\t1h\t10m\t100%\t1\t1\t1\t1\n'
    )
    rows = tcr.load_trace(p)
    assert rows == [
        {'name': 'DIAMOND_SEARCH (x)', 'status': 'COMPLETED', 'realtime_seconds': 45 * 60.0},
        {'name': 'HMMSEARCH (y)', 'status': 'FAILED', 'realtime_seconds': 10 * 60.0},
    ]


def test_cpu_hours_by_process_group_excludes_failed():
    rows = [
        {'name': 'DIAMOND_SEARCH (x)', 'status': 'COMPLETED', 'realtime_seconds': 3600.0},
        {'name': 'DIAMOND_SEARCH (y)', 'status': 'FAILED', 'realtime_seconds': 3600.0},
        {'name': 'HMMSEARCH (z)', 'status': 'CACHED', 'realtime_seconds': 1800.0},
    ]
    groups = {'diamond': ['DIAMOND_SEARCH'], 'hmm': ['HMMSEARCH']}
    result = tcr.cpu_hours_by_process_group(rows, groups)
    assert result == {'diamond': 1.0, 'hmm': 0.5}
