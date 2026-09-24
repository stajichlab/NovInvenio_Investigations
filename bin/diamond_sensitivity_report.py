#!/usr/bin/env python3
"""Compare Tier P runs that differ only in diamond sensitivity (default /
--sensitive / --very-sensitive) for one clade, against each other and against
the clade's Tier C+H (cluster + family HMM) run.

Benchmark design: notes/diamond-sensitivity/README.md.

For every mode and direction (gain = candidates.txt, loss = loss_candidates.txt):

  - n candidates; vs the default-mode run: shared, gained, lost, Jaccard.
  - vs Tier C+H: precision, recall, Jaccard (C+H measured against this mode).
  - ch_extras_recovered: of the C+H candidates that the default run rejects, the
    fraction this mode calls. A high value means the mode now finds the seed-group
    homologs that default diamond missed.
  - ch_misses_resolved: of the default-run candidates that C+H rejects, the
    fraction this mode no longer calls (it now agrees with C+H).

For each non-default mode, the candidates it LOSES relative to the default run
are the ones where the more sensitive search now finds an other-group hit. The
lost-candidate evidence table bins that best new other-group hit by E-value and
query coverage (qcovhsp), from the mode's raw search_cache. Many low-coverage or
weak hits would point at new paralog/fragment cross-hits rather than real
orthologs.

Cost per mode: the pairwise search (SEARCH:/LOSS_SEARCH: DIAMOND_SEARCH+MAKEDB),
DIAMOND_SELF, and the whole Tier P total, from the run's trace
(trace_cost_report.py's parser).

Outputs (all small, class 2): --output-prefix.{cost,concordance,lost_evidence}.tsv.
Optional --output-lost-detail lists the lost candidates themselves (class 3).
"""
import argparse
import csv
import gzip
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import trace_cost_report as tcr  # noqa: E402

NII_PIPELINE_LIB = Path('/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
sys.path.insert(0, str(NII_PIPELINE_LIB))
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config  # noqa: E402

DIRECTIONS = {'gain': 'candidates.txt', 'loss': 'loss_candidates.txt'}
E_BINS = [(1e-50, 'E<1e-50'), (1e-20, '1e-50<=E<1e-20'), (1e-5, '1e-20<=E<1e-5')]
COV_BINS = [(30, 'qcov<30'), (60, '30<=qcov<60'), (101, 'qcov>=60')]


def load_candidates(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else 1.0


def ratio(n, d):
    return round(n / d, 4) if d else ''


def compare(mode_set, default_set, ch_set):
    ch_extras = ch_set - default_set          # C+H calls, default P rejects
    ch_misses = default_set - ch_set          # default P calls, C+H rejects
    return {
        'n_candidates': len(mode_set),
        'shared_with_default': len(mode_set & default_set),
        'gained_vs_default': len(mode_set - default_set),
        'lost_vs_default': len(default_set - mode_set),
        'jaccard_vs_default': round(jaccard(mode_set, default_set), 4),
        'n_ch': len(ch_set),
        'ch_precision_vs_mode': ratio(len(ch_set & mode_set), len(ch_set)),
        'ch_recall_vs_mode': ratio(len(ch_set & mode_set), len(mode_set)),
        'jaccard_vs_ch': round(jaccard(mode_set, ch_set), 4),
        'n_ch_extras_vs_default': len(ch_extras),
        'ch_extras_recovered': ratio(len(ch_extras & mode_set), len(ch_extras)),
        'n_ch_misses_vs_default': len(ch_misses),
        'ch_misses_resolved': ratio(len(ch_misses - mode_set), len(ch_misses)),
    }


def e_bin(e):
    for cut, label in E_BINS:
        if e < cut:
            return label
    return 'E>=1e-5'


def cov_bin(q):
    for cut, label in COV_BINS:
        if q < cut:
            return label
    return COV_BINS[-1][1]


def best_other_hits(cache_dir, wanted_by_short, other_shorts, evalue_max):
    """pid -> (evalue, qcovhsp, target_short) for the best hit with E < evalue_max
    against any other-group proteome, from <Q>_vs_<T>.diamond.tsv.gz."""
    best = {}
    for q, pids in wanted_by_short.items():
        for t in other_shorts:
            f = Path(cache_dir) / f'{q}_vs_{t}.diamond.tsv.gz'
            if not f.exists():
                continue
            with gzip.open(f, 'rt') as fh:
                for line in fh:
                    a = line.rstrip('\n').split('\t')
                    if a[0] not in pids:
                        continue
                    e = float(a[2])
                    if e >= evalue_max:
                        continue
                    qcov = float(a[6]) if len(a) > 6 and a[6] else float('nan')
                    if a[0] not in best or e < best[a[0]][0]:
                        best[a[0]] = (e, qcov, t)
    return best


def search_cost(results_dir, trace_file=None):
    trace = trace_file or tcr.latest_trace(str(Path(results_dir) / 'nextflow_log'))
    rows = tcr.load_trace(trace)
    s_wall, s_cpu, n_s = tcr.cost_of(rows, lambda n: bool(tcr.PAIRWISE_SEARCH_RE.match(n)))
    self_wall, self_cpu, _ = tcr.cost_of(rows, lambda n: 'DIAMOND_SELF' in n)
    partial = ['DIAMOND_SELF', 'PARSE_HITS', 'PARSE_SELF_HITS', 'TBLASTN_MAKEDB', 'TBLASTN',
               'BUILD_PRESENCE_MATRIX']
    tot_wall, tot_cpu, _ = tcr.cost_of(
        rows, lambda n: bool(tcr.PAIRWISE_SEARCH_RE.match(n)) or any(p in n for p in partial))
    return {'trace': os.path.basename(trace), 'n_search_rows': n_s,
            'search_wall_h': s_wall, 'search_cpu_h': s_cpu,
            'self_wall_h': self_wall, 'self_cpu_h': self_cpu,
            'tier_p_wall_h': tot_wall, 'tier_p_cpu_h': tot_cpu}


def parse_trace_overrides(pairs):
    """[[mode, trace_file], ...] from --trace -> {mode: trace_file}."""
    return {m: t for m, t in (pairs or [])}


def write_tsv(path, rows):
    if not rows:
        return
    cols = list(rows[0].keys())
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--label', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--mode', action='append', nargs=2, metavar=('NAME', 'RESULTS_DIR'),
                    required=True, dest='modes',
                    help='Repeatable. The first --mode is the baseline (default diamond).')
    ap.add_argument('--trace', action='append', nargs=2, metavar=('NAME', 'TRACE_FILE'),
                    dest='traces', help='Repeatable. Trace file to cost a mode from, instead of '
                    'the newest one in its nextflow_log/ (use when a resumed run left a '
                    'later, partial trace).')
    ap.add_argument('--ch-dir', required=True, dest='ch_dir', help="Tier C+H run's results dir")
    ap.add_argument('--evalue', type=float, default=1e-5,
                    help='Significance cutoff used by the runs (for lost-candidate evidence)')
    ap.add_argument('--output-prefix', required=True, dest='output_prefix')
    ap.add_argument('--output-lost-detail', default=None, dest='output_lost_detail')
    args = ap.parse_args()

    samples = parse_config(args.config)
    ingroup = {s.short for s in samples if s.group in INGROUP_ROLES}
    outgroup = {s.short for s in samples if s.group in OUTGROUP_ROLES}
    base_name, base_dir = args.modes[0]

    cost_rows, conc_rows, ev_rows, detail = [], [], [], []
    traces = parse_trace_overrides(args.traces)
    for name, rdir in args.modes:
        c = search_cost(rdir, traces.get(name))
        cost_rows.append({'label': args.label, 'mode': name, **c})

    for direction, fname in DIRECTIONS.items():
        other = outgroup if direction == 'gain' else ingroup
        base = load_candidates(Path(base_dir) / fname)
        ch = load_candidates(Path(args.ch_dir) / fname)
        for name, rdir in args.modes:
            mode_set = load_candidates(Path(rdir) / fname)
            conc_rows.append({'label': args.label, 'direction': direction, 'mode': name,
                              **compare(mode_set, base, ch)})
            if name == base_name:
                continue
            lost = base - mode_set
            wanted = defaultdict(set)
            for line in lost:
                short, _, pid = line.partition('::')
                wanted[short].add(pid)
            hits = best_other_hits(Path(rdir) / 'search_cache', wanted, other, args.evalue)
            bins = Counter()
            for line in sorted(lost):
                short, _, pid = line.partition('::')
                h = hits.get(pid)
                if h is None:
                    key = ('no significant raw hit', '')
                else:
                    key = (e_bin(h[0]), cov_bin(h[1]) if h[1] == h[1] else 'qcov n/a')
                bins[key] += 1
                detail.append((args.label, direction, name, line,
                               '' if h is None else h[0], '' if h is None else h[1],
                               '' if h is None else h[2]))
            for (eb, cb), n in sorted(bins.items()):
                ev_rows.append({'label': args.label, 'direction': direction, 'mode': name,
                                'n_lost': len(lost), 'evalue_bin': eb, 'qcov_bin': cb, 'n': n,
                                'frac': ratio(n, len(lost))})

    write_tsv(f'{args.output_prefix}.cost.tsv', cost_rows)
    write_tsv(f'{args.output_prefix}.concordance.tsv', conc_rows)
    write_tsv(f'{args.output_prefix}.lost_evidence.tsv', ev_rows)
    if args.output_lost_detail:
        with open(args.output_lost_detail, 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['label', 'direction', 'mode', 'candidate', 'best_other_evalue',
                        'best_other_qcov', 'best_other_proteome'])
            w.writerows(detail)
    print(f'{args.label}: wrote {args.output_prefix}.{{cost,concordance,lost_evidence}}.tsv',
          file=sys.stderr)


if __name__ == '__main__':
    main()
