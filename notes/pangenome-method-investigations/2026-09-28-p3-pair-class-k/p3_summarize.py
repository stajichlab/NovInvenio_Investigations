#!/usr/bin/env python3.12
"""Summarize the #132 priority-3 sweep (run_p3_sweep.sh outputs).

For every dataset and setting: label counts, and against the default
(k=10, physical 0.5, trans 0.05) the share of pairs whose label is unchanged,
plus the flows between physical (unexplained_physical, starship_explained),
ambiguous_linkage and low-linkage (trans, trans_unconfirmed) groups.
Also the linkage_fraction distribution at the default k (how many pairs sit
near each threshold).

Usage: p3_summarize.py --sweep <dir with <dataset>/k*_p*_t*.tsv.zst> --out p3_summary.tsv
"""
import argparse
import collections
import glob
import io
import os
import re
import subprocess

GROUP = {'unexplained_physical': 'physical', 'starship_explained': 'physical',
         'ambiguous_linkage': 'ambiguous', 'trans': 'low', 'trans_unconfirmed': 'low',
         'insufficient_data': 'insufficient'}


def read(path):
    txt = subprocess.run(['zstd', '-dc', path], capture_output=True, check=True, text=True).stdout
    fh = io.StringIO(txt)
    head = fh.readline().rstrip('\n').split('\t')
    ia, ib, ic, il = (head.index(x) for x in ('family_a', 'family_b', 'classification', 'linkage_fraction'))
    lab, lf = {}, {}
    for line in fh:
        p = line.rstrip('\n').split('\t')
        key = (p[ia], p[ib])
        lab[key] = p[ic]
        lf[key] = p[il]
    return lab, lf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sweep', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    rows = []
    for ddir in sorted(glob.glob(os.path.join(args.sweep, '*/'))):
        name = os.path.basename(ddir.rstrip('/'))
        files = sorted(glob.glob(ddir + 'k*_p*_t*.tsv.zst'))
        if not files:
            continue
        base_f = ddir + 'k10_p0.5_t0.05.tsv.zst'
        base, base_lf = read(base_f)
        # linkage_fraction near the thresholds at the default k
        vals = [float(v) for v in base_lf.values() if v not in ('', 'NA', 'nan')]
        near = lambda c, d: sum(1 for v in vals if abs(v - c) <= d)
        print(f'{name}: {len(base)} pairs; linkage_fraction measured for {len(vals)}; '
              f'exactly 0: {sum(1 for v in vals if v == 0)}; within 0.05 of 0.5: {near(0.5, 0.05)}; '
              f'within 0.02 of 0.05: {near(0.05, 0.02)}')
        for f in files:
            k, p, t = re.match(r'.*k(\d+)_p([\d.]+)_t([\d.]+)\.tsv\.zst', f).groups()
            lab, _ = read(f) if f != base_f else (base, None)
            cnt = collections.Counter(lab.values())
            same = sum(1 for key, v in lab.items() if base.get(key) == v) / len(lab)
            flow = collections.Counter((GROUP[base[key]], GROUP[v]) for key, v in lab.items()
                                       if key in base and GROUP[base[key]] != GROUP[v])
            rows.append({
                'dataset': name, 'k': int(k), 'physical_threshold': float(p), 'trans_threshold': float(t),
                'pairs': len(lab),
                **{c: cnt.get(c, 0) for c in ('unexplained_physical', 'starship_explained', 'ambiguous_linkage',
                                              'trans', 'trans_unconfirmed', 'insufficient_data')},
                'same_label_as_default': round(same, 4),
                'group_changes': '; '.join(f'{a}->{b}:{n}' for (a, b), n in flow.most_common()),
            })
    rows.sort(key=lambda r: (r['dataset'], r['physical_threshold'] != 0.5 or r['trans_threshold'] != 0.05,
                             r['k'], r['physical_threshold'], r['trans_threshold']))
    with open(args.out, 'w') as fh:
        fh.write('\t'.join(rows[0]) + '\n')
        for r in rows:
            fh.write('\t'.join(str(v) for v in r.values()) + '\n')


if __name__ == '__main__':
    main()
