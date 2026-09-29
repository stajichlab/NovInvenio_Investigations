#!/usr/bin/env python3.12
"""#132 priority 2: do the core / soft-core / shell cutoffs (0.95 / 0.90 / 0.15)
sit on natural breaks in the family frequency spectrum?

Input: one or more pipeline frequency_table.tsv files (per-group bins, #212).
Families counted in >= 1 representative of the ingroup (strain_count >= 1;
old-format tables also put strain_count 0 in 'singleton') are used; k = strain_count, n = number of
ingroup representatives (recovered as k / frequency).

Per dataset it writes:
  - the spectrum (families per k)
  - integer thresholds each cutoff maps to at this n
  - class counts when each cutoff moves (sensitivity)
  - a zero-truncated binomial mixture (Snipen et al. 2009, BMC Genomics 10:385)
    fitted by EM for K = 2..10 components; BIC choice; for every K the smallest
    k whose most likely component is the highest-p component (the data's own
    core boundary) and the largest k whose most likely component is the
    lowest-p component (the data's own cloud boundary)

Usage:
  p2_frequency_cutoffs.py --table NAME=path/frequency_table.tsv ... --outdir .
"""
import argparse
import csv
import math
import os

import numpy as np

COUNTED = {'singleton', 'cloud', 'shell', 'soft_core', 'core'}
CUTS = {'core': 0.95, 'soft_core': 0.90, 'shell': 0.15}


def read_table(path):
    ks, fs = [], []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['bin'] not in COUNTED:
                continue
            k, f = int(r['strain_count']), float(r['frequency'])
            if k < 1:       # old-format tables put k = 0 (no ingroup representative) in 'singleton'
                continue
            ks.append(k)
            fs.append(f)
    ks, fs = np.array(ks), np.array(fs)
    ratio = ks[fs > 0.2] / fs[fs > 0.2]
    n = int(round(float(np.median(ratio))))
    return ks, n


def log_binom(n):
    k = np.arange(n + 1)
    return np.array([math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) for i in k])


def solve_p(mean_k, n):
    """p such that the zero-truncated binomial mean n p / (1 - (1-p)^n) = mean_k."""
    if mean_k <= 1 + 1e-9:
        return 1e-6
    if mean_k >= n - 1e-9:
        return 1 - 1e-9
    lo, hi = 1e-9, 1 - 1e-12
    for _ in range(100):
        mid = (lo + hi) / 2
        m = n * mid / (1 - (1 - mid) ** n)
        lo, hi = (mid, hi) if m < mean_k else (lo, mid)
    return (lo + hi) / 2


def fit_mixture(counts, n, K, seed=0, iters=2000):
    """EM for a K-component zero-truncated binomial mixture on aggregated counts.

    counts[k] = number of families seen in k strains, k = 1..n.
    Returns (loglik, weights, p) with components sorted by p.
    """
    k = np.arange(1, n + 1)
    c = counts[1:].astype(float)
    lb = log_binom(n)[1:]
    rng = np.random.default_rng(seed)
    p = np.sort(rng.uniform(0.01, 0.99, K))
    p[0], p[-1] = min(p[0], 1.5 / n), max(p[-1], 1 - 0.5 / n)
    w = np.full(K, 1 / K)
    prev = -np.inf
    for _ in range(iters):
        logf = (lb[:, None] + k[:, None] * np.log(p)[None, :]
                + (n - k)[:, None] * np.log1p(-p)[None, :]
                - np.log1p(-(1 - p) ** n)[None, :])
        logw = logf + np.log(w)[None, :]
        mx = logw.max(1, keepdims=True)
        lse = mx[:, 0] + np.log(np.exp(logw - mx).sum(1))
        ll = float((c * lse).sum())
        r = np.exp(logw - lse[:, None]) * c[:, None]
        tot = r.sum(0)
        w = np.maximum(tot / c.sum(), 1e-12)
        p = np.array([solve_p(float((r[:, j] * k).sum() / tot[j]), n) if tot[j] > 1e-9 else p[j]
                      for j in range(K)])
        p = np.clip(p, 1e-9, 1 - 1e-9)
        if ll - prev < 1e-7 * abs(ll):
            break
        prev = ll
    o = np.argsort(p)
    return ll, w[o], p[o], (logw - lse[:, None])[:, o]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', action='append', required=True, help='NAME=path')
    ap.add_argument('--outdir', default='.')
    ap.add_argument('--kmax', type=int, default=10)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    summary, mix_rows, sens_rows, spectra = [], [], [], {}
    for spec in args.table:
        name, path = spec.split('=', 1)
        ks, n = read_table(path)
        counts = np.bincount(ks, minlength=n + 1)[:n + 1]
        spectra[name] = (n, counts)
        N = int(counts.sum())
        thr = {c: math.ceil(v * n - 1e-9) for c, v in CUTS.items()}
        # class counts under the pipeline rule (singleton k=1 first)
        core = int(counts[thr['core']:].sum())
        soft = int(counts[thr['soft_core']:thr['core']].sum())
        shell = int(counts[max(thr['shell'], 2):thr['soft_core']].sum())
        cloud = int(counts[2:max(thr['shell'], 2)].sum())
        single = int(counts[1])
        # a trough: smoothed spectrum, minimum over 0.1..0.9, and its depth
        wlen = max(1, round(0.05 * n))
        sm = np.convolve(counts[1:].astype(float), np.ones(wlen) / wlen, mode='same')
        f = np.arange(1, n + 1) / n
        mid = (f >= 0.1) & (f <= 0.9)
        imin = np.flatnonzero(mid)[np.argmin(sm[mid])]
        band = (f >= 0.2) & (f <= 0.8)
        flat = float(sm[band].max() / max(sm[band].min(), 1e-9)) if band.sum() > 1 else float('nan')
        summary.append({
            'dataset': name, 'n_reps': n, 'families': N,
            'k_core>=': thr['core'], 'k_soft>=': thr['soft_core'], 'k_shell>=': max(thr['shell'], 2),
            'singleton': single, 'cloud': cloud, 'shell': shell, 'soft_core': soft, 'core': core,
            'trough_freq': round(float(f[imin]), 3),
            'max/min_0.2-0.8': round(flat, 2),
            'frac_in_k=n': round(counts[n] / N, 3),
            'frac_in_k=1': round(counts[1] / N, 3),
        })
        # sensitivity: class sizes as each cutoff moves
        for label, grid in (('core', [0.90, 0.93, 0.95, 0.97, 0.99, 1.0]),
                            ('shell', [0.05, 0.10, 0.15, 0.20, 0.25])):
            for g in grid:
                t = math.ceil(g * n - 1e-9)
                v = int(counts[t:].sum()) if label == 'core' else int(counts[max(t, 2):thr['soft_core']].sum())
                sens_rows.append({'dataset': name, 'class': label, 'cutoff': g, 'k_threshold': max(t, 1 if label == 'core' else 2), 'families': v})
        # mixture
        fits = []
        for K in range(2, args.kmax + 1):
            best = None
            for seed in range(5):
                res = fit_mixture(counts, n, K, seed)
                if best is None or res[0] > best[0]:
                    best = res
            ll, w, p, post = best
            bic = -2 * ll + (2 * K - 1) * math.log(N)
            arg = post.argmax(1)                 # most likely component per k = 1..n
            top = np.flatnonzero(arg == K - 1)
            low = np.flatnonzero(arg == 0)
            core_k = int(top.min() + 1) if len(top) else None
            cloud_k = int(low.max() + 1) if len(low) else None
            fits.append((bic, K))
            mix_rows.append({'dataset': name, 'K': K, 'loglik': round(ll, 1), 'BIC': round(bic, 1),
                             'p_top': round(float(p[-1]), 4), 'w_top': round(float(w[-1]), 4),
                             'p_low': round(float(p[0]), 4), 'w_low': round(float(w[0]), 4),
                             'core_boundary_k': core_k,
                             'core_boundary_freq': round(core_k / n, 3) if core_k else None,
                             'low_component_max_k': cloud_k,
                             'low_component_max_freq': round(cloud_k / n, 3) if cloud_k else None,
                             'p_all': ' '.join(f'{x:.3g}' for x in p)})
        summary[-1]['BIC_K'] = min(fits)[1]
        print(f'{name}: n={n}, {N} families, BIC K={min(fits)[1]}', flush=True)

    def write(fn, rows):
        with open(os.path.join(args.outdir, fn), 'w', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
            wr.writeheader()
            wr.writerows(rows)
    write('p2_summary.tsv', summary)
    write('p2_sensitivity.tsv', sens_rows)
    write('p2_mixture.tsv', mix_rows)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    m = len(spectra)
    cols = 4 if m > 4 else m
    rows_ = math.ceil(m / cols)
    fig, axes = plt.subplots(rows_, cols, figsize=(4 * cols, 3.1 * rows_), squeeze=False)
    bic_k = {s['dataset']: s['BIC_K'] for s in summary}
    for ax, (name, (n, counts)) in zip(axes.flat, spectra.items()):
        f = np.arange(1, n + 1) / n
        ax.bar(f, counts[1:], width=1 / n, color='#2b6f8e', linewidth=0)
        ax.set_yscale('log')
        for v, ls in ((0.15, ':'), (0.90, '--'), (0.95, '-')):
            ax.axvline(v, color='#b3452c', ls=ls, lw=1)
        row = [r for r in mix_rows if r['dataset'] == name and r['K'] == bic_k[name]][0]
        if row['core_boundary_freq']:
            ax.axvline(row['core_boundary_freq'], color='#3d8a5a', lw=1.4)
        ax.set_title(f'{name} (n={n})', fontsize=10)
        ax.set_xlabel('frequency (share of representatives)', fontsize=8)
        ax.set_ylabel('families', fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[m:]:
        ax.axis('off')
    fig.suptitle('Family frequency spectra. Red: 0.15 (dotted), 0.90 (dashed), 0.95 (solid). '
                 'Green: mixture core boundary at the BIC-chosen K', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, 'p2_spectra.png'), dpi=130)


if __name__ == '__main__':
    main()
