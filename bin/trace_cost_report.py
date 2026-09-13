#!/usr/bin/env python3
"""Analysis 3: real per-tier compute cost from nextflow trace files already published
in each run's nextflow_log/ -- see notes/superpowers/specs/2026-09-13-cluster-vs-
pairwise-sensitivity-design.md.

Two important caveats discovered while building this report (see task-8-report.md):

1. Tier P's all-vs-all diamond search (DIAMOND_SEARCH / DIAMOND_MAKEDB) uses
   Nextflow's `storeDir` directive, not ordinary resume-caching. A storeDir cache hit
   means Nextflow skips the task so completely it never logs a trace row for it at
   all (unlike -resume, which still logs CACHED). None of pezizo_set1's trace files
   contain any DIAMOND_SEARCH/DIAMOND_MAKEDB row, because the all-vs-all search was
   already cached in search_cache/ before any trace history we have begins. That
   true cost is therefore UNMEASURABLE from any available trace file. What IS
   measurable for Tier P (DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, TBLASTN,
   TBLASTN_MAKEDB, BUILD_PRESENCE_MATRIX) is reported as a labeled partial /
   lower-bound cost, not the true total -- see the `note` column in the output TSV.

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


def load_trace(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            rows.append({
                'name': row['name'],
                'status': row['status'],
                'realtime_seconds': parse_duration(row['realtime']),
            })
    return rows


def cpu_hours_by_process_group(rows, groups):
    totals = {name: 0.0 for name in groups}
    for row in rows:
        if row['status'] not in ('COMPLETED', 'CACHED'):
            continue
        for group_name, substrings in groups.items():
            if any(sub in row['name'] for sub in substrings):
                totals[group_name] += row['realtime_seconds'] / 3600
                break
    return {k: round(v, 3) for k, v in totals.items()}


def cpu_hours_by_prefix(rows, prefix):
    """Sum CPU-hours for every COMPLETED/CACHED row whose name starts with `prefix`.

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
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    pairwise_trace = args.pairwise_trace_file or latest_trace(args.pairwise_trace_dir)
    cluster_trace = args.cluster_trace_file or latest_trace(args.cluster_trace_dir)
    pairwise_rows = load_trace(pairwise_trace)
    cluster_rows = load_trace(cluster_trace)

    # Tier P: the true all-vs-all DIAMOND_SEARCH/DIAMOND_MAKEDB cost is unmeasurable
    # (storeDir cache hit -> no trace row ever logged, confirmed absent from every
    # pezizo_set1 trace file including the earliest). What we report here is a
    # labeled partial/lower-bound: everything else in Tier P's pipeline that DOES
    # appear in the trace.
    pairwise_groups = {
        'search_partial': ['DIAMOND_SELF', 'PARSE_HITS', 'PARSE_SELF_HITS',
                           'TBLASTN_MAKEDB', 'TBLASTN', 'BUILD_PRESENCE_MATRIX'],
    }
    p_cost = cpu_hours_by_process_group(pairwise_rows, pairwise_groups)

    # Tier C+H: prefix match on 'PROFILE_SEARCH:' captures the entire gain-side
    # family-profile pathway in one robust check (see module docstring point 2).
    ch_cost = cpu_hours_by_prefix(cluster_rows, 'PROFILE_SEARCH:')
    # Loss-side equivalent, reported for completeness (not required, but consistent).
    loss_cost = cpu_hours_by_prefix(cluster_rows, 'PROFILE_LOSS_SEARCH:')

    rows = [
        {'tier': 'P', 'cpu_hours': p_cost['search_partial'],
         'note': ('PARTIAL/lower-bound only -- true all-vs-all DIAMOND_SEARCH cost is '
                   'unmeasurable: storeDir cache hit means Nextflow never logs a trace '
                   'row for it (unlike ordinary -resume CACHED rows). This total covers '
                   'only DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, TBLASTN(+MAKEDB), '
                   'BUILD_PRESENCE_MATRIX -- do NOT read this as Tier P\'s full cost.')},
        {'tier': 'C+H', 'cpu_hours': ch_cost,
         'note': 'PROFILE_SEARCH:* (gain-side family-profile pathway), measured directly.'},
        {'tier': 'C+H_loss', 'cpu_hours': loss_cost,
         'note': 'PROFILE_LOSS_SEARCH:* (loss-side equivalent), measured directly.'},
        {'tier': 'C', 'cpu_hours': 0.0,
         'note': 'Free byproduct of Tier C+H\'s own clustering step.'},
        {'tier': 'R', 'cpu_hours': round((args.tier_r_diamond_seconds or 0.0) / 3600, 3),
         'note': ('Measured directly with `time` around refine_ambiguous_families.py\'s '
                  'diamond step -- not read from any trace file (plain subprocess call, '
                  'no trace row).' if args.tier_r_diamond_seconds is not None
                  else 'No --tier-r-diamond-seconds given; reported as 0.0, not measured.')},
    ]
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['tier', 'cpu_hours', 'note'], delimiter='\t',
                            lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
