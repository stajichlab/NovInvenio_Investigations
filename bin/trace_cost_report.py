#!/usr/bin/env python3
"""Analysis 3: real per-tier compute cost from nextflow trace files already published
in each run's nextflow_log/ -- see notes/superpowers/specs/2026-09-13-cluster-vs-
pairwise-sensitivity-design.md.

Two important caveats discovered while building this report (see task-8-report.md):

1. Tier P's all-vs-all search (DIAMOND_SEARCH / DIAMOND_MAKEDB, or the BLAST/PHMMER
   equivalents) uses Nextflow's `storeDir` directive, not ordinary resume-caching. A
   storeDir cache hit skips the task so completely that Nextflow logs no trace row for
   it (unlike -resume, which still logs CACHED). If the search was already in
   search_cache/ before the trace began (pezizo_set1), its cost is UNMEASURABLE from
   that trace, and the P row is a labeled partial / lower bound. If the search ran in
   the traced run (sordariales_shallow, 2026-09-22), its rows ARE in the trace, and
   the P row is the complete cost. main() checks which case applies from the trace
   rows themselves; it does not assume either one.

2. Tier C+H's real cost is not well captured by a short substring list (the plan's
   original suggestion of ['MMSEQS', 'BUILD_FAMILY_PROFILES', 'HMMSEARCH',
   'BUILD_CHUNK'] misses PROFILE_SEARCH:SEED_PROTEIN_MAP and
   PROFILE_SEARCH:PROFILE_PRESENCE_MATRIX). Every real process name under the
   gain-side family-profile pathway is prefixed `PROFILE_SEARCH:`, so main() matches
   on that prefix instead of a substring list.
"""
import argparse
import csv
import glob
import re


DURATION_RE = re.compile(r'(?:(\d+)h)?\s*(?:(\d+)m(?!s))?\s*(?:(\d+(?:\.\d+)?)s)?\s*(?:(\d+)ms)?')


def parse_duration(s):
    s = s.strip()
    if s in ('', '-', '0'):
        return 0.0
    h, m, sec, ms = DURATION_RE.match(s).groups()
    total = 0.0
    if h:
        total += int(h) * 3600
    if m:
        total += int(m) * 60
    if sec:
        total += float(sec)
    if ms:
        total += int(ms) / 1000
    return total


def parse_cpu_pct(s):
    """'2400.5%' -> 24.005 (cores used on average). Missing/'-' -> 0.0."""
    s = (s or '').strip().rstrip('%')
    if s in ('', '-'):
        return 0.0
    return float(s) / 100


# Tier P's all-vs-all search processes: <SEARCH|LOSS_SEARCH>:<TOOL>_<SEARCH|MAKEDB>.
# DIAMOND_SELF (the self-search) is deliberately excluded -- it is already in the
# partial group below.
PAIRWISE_SEARCH_RE = re.compile(r'^(?:LOSS_)?SEARCH:(?:DIAMOND|BLAST|PHMMER)_(?:SEARCH|MAKEDB)\b')


def load_trace(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            rows.append({
                'name': row['name'],
                'status': row['status'],
                'realtime_seconds': parse_duration(row['realtime']),
                'cpu_fraction': parse_cpu_pct(row.get('%cpu', '')),
            })
    return rows


def wall_hours_by_process_group(rows, groups):
    totals = {name: 0.0 for name in groups}
    for row in rows:
        if row['status'] not in ('COMPLETED', 'CACHED'):
            continue
        for group_name, substrings in groups.items():
            if any(sub in row['name'] for sub in substrings):
                totals[group_name] += row['realtime_seconds'] / 3600
                break
    return {k: round(v, 3) for k, v in totals.items()}


def wall_hours_by_prefix(rows, prefix):
    """Sum wall-clock hours for every COMPLETED/CACHED row whose name starts with `prefix`.

    Used instead of a substring list for Tier C+H's cost: every real process name
    under the gain-side family-profile pathway (MMSEQS_FAMILY_CLUSTER,
    BUILD_FAMILY_PROFILES:*, FAMILY_HMMSEARCH:*, SEED_PROTEIN_MAP,
    PROFILE_PRESENCE_MATRIX) shares the `PROFILE_SEARCH:` prefix, so this is a single
    robust check instead of a brittle, incomplete substring guess.
    """
    total = 0.0
    for row in rows:
        if row['status'] not in ('COMPLETED', 'CACHED'):
            continue
        if row['name'].startswith(prefix):
            total += row['realtime_seconds'] / 3600
    return round(total, 3)


def cost_of(rows, predicate):
    """(wall_hours, cpu_hours, n_rows) over COMPLETED/CACHED rows matching predicate.
    cpu_hours = realtime * average cores used (%cpu / 100)."""
    wall = cpu = 0.0
    n = 0
    for row in rows:
        if row['status'] not in ('COMPLETED', 'CACHED') or not predicate(row['name']):
            continue
        wall += row['realtime_seconds'] / 3600
        cpu += row['realtime_seconds'] * row.get('cpu_fraction', 0.0) / 3600
        n += 1
    return round(wall, 3), round(cpu, 3), n


def latest_trace(nextflow_log_dir):
    traces = sorted(glob.glob(f'{nextflow_log_dir}/*-trace.txt'))
    if not traces:
        raise FileNotFoundError(f'no *-trace.txt files under {nextflow_log_dir}')
    return traces[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pairwise-trace-dir', required=True, dest='pairwise_trace_dir',
                    help="pairwise run's nextflow_log/ directory")
    ap.add_argument('--cluster-trace-dir', required=True, dest='cluster_trace_dir',
                    help="cluster+HMM run's nextflow_log/ directory (Tier R's own new "
                         "within-family diamond step is timed separately, not read "
                         "from this trace -- see --tier-r-diamond-seconds)")
    ap.add_argument('--pairwise-trace-file', default=None, dest='pairwise_trace_file',
                    help='Explicit trace file to use instead of the alphabetically-last '
                         'one under --pairwise-trace-dir. Use this when the newest trace '
                         "file in the directory is a later, unrelated/partial rerun that "
                         "doesn't reflect the run whose published outputs you're costing "
                         "(check output-file mtimes against each trace file's process set "
                         'and status counts before trusting "latest" == "the run that '
                         'matters" -- see task-8-report.md for a real case where this bit).')
    ap.add_argument('--cluster-trace-file', default=None, dest='cluster_trace_file',
                    help='Same override as --pairwise-trace-file, for the cluster+HMM run.')
    ap.add_argument('--tier-r-diamond-seconds', type=float, default=None,
                    dest='tier_r_diamond_seconds',
                    help='Wall-clock of the refine_ambiguous_families.py diamond step '
                         '(time it separately with `time` when running Task 5/6 Step 6 '
                         '-- there is no trace file for a plain subprocess call)')
    ap.add_argument('--tier-r-loss-seconds', type=float, default=None,
                    dest='tier_r_loss_seconds',
                    help='Same as --tier-r-diamond-seconds, for the loss-direction '
                         'refine_ambiguous_families.py run (--query-group OUT)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    pairwise_trace = args.pairwise_trace_file or latest_trace(args.pairwise_trace_dir)
    cluster_trace = args.cluster_trace_file or latest_trace(args.cluster_trace_dir)
    pairwise_rows = load_trace(pairwise_trace)
    cluster_rows = load_trace(cluster_trace)

    partial_subs = ['DIAMOND_SELF', 'PARSE_HITS', 'PARSE_SELF_HITS',
                    'TBLASTN_MAKEDB', 'TBLASTN', 'BUILD_PRESENCE_MATRIX']
    p_part_wall, p_part_cpu, _ = cost_of(
        pairwise_rows,
        lambda n: not PAIRWISE_SEARCH_RE.match(n) and any(sub in n for sub in partial_subs))
    p_search_wall, p_search_cpu, n_search = cost_of(
        pairwise_rows, lambda n: bool(PAIRWISE_SEARCH_RE.match(n)))
    if n_search:
        p_note = (f'COMPLETE: all-vs-all search measured from {n_search} '
                  f'*_SEARCH/*_MAKEDB trace rows ({p_search_wall} wall-h, '
                  f'{p_search_cpu} cpu-h), plus DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, '
                  'TBLASTN(+MAKEDB), BUILD_PRESENCE_MATRIX.')
    else:
        p_note = ('PARTIAL/lower-bound only -- this trace has no *_SEARCH/*_MAKEDB rows '
                  '(storeDir cache hit: Nextflow logs no trace row for it, unlike ordinary '
                  '-resume CACHED rows), so the all-vs-all search cost is unmeasurable '
                  'from this trace. This total covers only DIAMOND_SELF, PARSE_HITS, '
                  "PARSE_SELF_HITS, TBLASTN(+MAKEDB), BUILD_PRESENCE_MATRIX -- do NOT read "
                  "this as Tier P's full cost.")

    # Tier C+H: prefix match captures each direction's whole family-profile pathway
    # (see module docstring point 2).
    ch_wall, ch_cpu, _ = cost_of(cluster_rows, lambda n: n.startswith('PROFILE_SEARCH:'))
    loss_wall, loss_cpu, _ = cost_of(cluster_rows, lambda n: n.startswith('PROFILE_LOSS_SEARCH:'))

    def r_row(tier, seconds, direction):
        if seconds is None:
            return {'tier': tier, 'wall_hours': 0.0, 'cpu_hours': '',
                    'note': f'No {direction} Tier R time given; reported as 0.0, not measured.'}
        return {'tier': tier, 'wall_hours': round(seconds / 3600, 3), 'cpu_hours': '',
                'note': ("Measured directly with `time` around refine_ambiguous_families.py "
                         f"({direction}) -- not read from any trace file (plain subprocess "
                         'call, no trace row). cpu_hours not measured.')}

    rows = [
        {'tier': 'P', 'wall_hours': round(p_part_wall + p_search_wall, 3),
         'cpu_hours': round(p_part_cpu + p_search_cpu, 3), 'note': p_note},
        {'tier': 'C+H', 'wall_hours': ch_wall, 'cpu_hours': ch_cpu,
         'note': 'PROFILE_SEARCH:* (gain-side family-profile pathway), measured directly.'},
        {'tier': 'C+H_loss', 'wall_hours': loss_wall, 'cpu_hours': loss_cpu,
         'note': 'PROFILE_LOSS_SEARCH:* (loss-side equivalent), measured directly.'},
        {'tier': 'C', 'wall_hours': 0.0, 'cpu_hours': 0.0,
         'note': 'Free byproduct of Tier C+H\'s own clustering step.'},
        r_row('R', args.tier_r_diamond_seconds, 'gain'),
    ]
    if args.tier_r_loss_seconds is not None:
        rows.append(r_row('R_loss', args.tier_r_loss_seconds, 'loss'))
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['tier', 'wall_hours', 'cpu_hours', 'note'], delimiter='\t',
                            lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
