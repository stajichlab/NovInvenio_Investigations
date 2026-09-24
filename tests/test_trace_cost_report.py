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
        {'name': 'DIAMOND_SEARCH (x)', 'status': 'COMPLETED', 'realtime_seconds': 45 * 60.0,
         'cpu_fraction': 1.0},
        {'name': 'HMMSEARCH (y)', 'status': 'FAILED', 'realtime_seconds': 10 * 60.0,
         'cpu_fraction': 1.0},
    ]


def test_wall_hours_by_process_group_excludes_failed():
    rows = [
        {'name': 'DIAMOND_SEARCH (x)', 'status': 'COMPLETED', 'realtime_seconds': 3600.0},
        {'name': 'DIAMOND_SEARCH (y)', 'status': 'FAILED', 'realtime_seconds': 3600.0},
        {'name': 'HMMSEARCH (z)', 'status': 'CACHED', 'realtime_seconds': 1800.0},
    ]
    groups = {'diamond': ['DIAMOND_SEARCH'], 'hmm': ['HMMSEARCH']}
    result = tcr.wall_hours_by_process_group(rows, groups)
    assert result == {'diamond': 1.0, 'hmm': 0.5}


HEADER = ('task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime'
          '\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n')


def _row(name, status='COMPLETED', realtime='1h', cpu='200.0%'):
    return f'1\ta\tb\t{name}\t{status}\t0\t2026-01-01\t1h\t{realtime}\t{cpu}\t1\t1\t1\t1\n'


def test_parse_cpu_pct():
    assert tcr.parse_cpu_pct('2400.0%') == 24.0
    assert tcr.parse_cpu_pct('-') == 0.0
    assert tcr.parse_cpu_pct('') == 0.0


def test_pairwise_search_re_excludes_self_search():
    assert tcr.PAIRWISE_SEARCH_RE.match('SEARCH:DIAMOND_SEARCH (A_vs_B)')
    assert tcr.PAIRWISE_SEARCH_RE.match('LOSS_SEARCH:DIAMOND_MAKEDB (B)')
    assert not tcr.PAIRWISE_SEARCH_RE.match('SEARCH:DIAMOND_SELF (A)')
    assert not tcr.PAIRWISE_SEARCH_RE.match('PROFILE_SEARCH:FAMILY_HMMSEARCH:HMMSEARCH_CHUNK (1)')


def _run_main(tmp_path, monkeypatch, pairwise_rows, extra_args=()):
    pw = tmp_path / 'pw'
    cl = tmp_path / 'cl'
    pw.mkdir()
    cl.mkdir()
    (pw / 'run-trace.txt').write_text(HEADER + ''.join(pairwise_rows))
    (cl / 'run-trace.txt').write_text(HEADER + _row('PROFILE_SEARCH:BUILD_CHUNK (1)', cpu='400.0%')
                                      + _row('PROFILE_LOSS_SEARCH:BUILD_CHUNK (1)', realtime='2h'))
    out = tmp_path / 'cost.tsv'
    monkeypatch.setattr(sys, 'argv', ['trace_cost_report.py', '--pairwise-trace-dir', str(pw),
                                      '--cluster-trace-dir', str(cl), '--output', str(out),
                                      *extra_args])
    tcr.main()
    import csv
    with open(out) as fh:
        return {r['tier']: r for r in csv.DictReader(fh, delimiter='\t')}


def test_main_complete_when_search_rows_present(tmp_path, monkeypatch):
    rows = _run_main(tmp_path, monkeypatch, [
        _row('SEARCH:DIAMOND_SEARCH (A)', realtime='30m', cpu='1000.0%'),
        _row('SEARCH:DIAMOND_SELF (A)', realtime='30m', cpu='100.0%'),
        _row('LOSS_SEARCH:DIAMOND_SEARCH (B)', status='FAILED'),
    ], extra_args=['--tier-r-diamond-seconds', '36', '--tier-r-loss-seconds', '360'])
    assert rows['P']['note'].startswith('COMPLETE')
    assert float(rows['P']['wall_hours']) == 1.0      # 0.5 search + 0.5 self; FAILED excluded
    assert float(rows['P']['cpu_hours']) == 5.5       # 0.5*10 + 0.5*1
    assert float(rows['C+H']['cpu_hours']) == 4.0
    assert float(rows['C+H_loss']['wall_hours']) == 2.0
    assert float(rows['R']['wall_hours']) == 0.01
    assert float(rows['R_loss']['wall_hours']) == 0.1


def test_main_partial_when_search_rows_absent(tmp_path, monkeypatch):
    rows = _run_main(tmp_path, monkeypatch, [_row('SEARCH:DIAMOND_SELF (A)')])
    assert rows['P']['note'].startswith('PARTIAL')
    assert float(rows['P']['wall_hours']) == 1.0
    assert 'R_loss' not in rows
