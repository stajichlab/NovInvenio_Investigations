#!/usr/bin/env python3
"""#132 priority 4: Leiden resolution x seed stability for the trans-pair module step.

Uses the pipeline's own functions (bin/pangenome_detect_trans_modules.py at the
pinned commit): load_trans_edges, build_graph, detect_modules
(RBConfigurationVertexPartition, weight = -log10 fdr_q clipped at 300).

Per resolution, over seeds 0..S-1:
  modules (count), modules with >= 10 families, largest module share of families,
  mean pairwise ARI and NMI between seeds (partition stability), and
  ARI / NMI against the next lower resolution at seed 0 (plateau detection).

Usage (pinned env python):
  p4_leiden_sweep.py --code <clone> --pairs NAME=pair_classification.tsv.zst ... --out p4_leiden.tsv
"""
import argparse
import itertools
import sys

import numpy as np


def contingency(a, b):
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    key = ia.astype(np.int64) * len(ub) + ib
    _, cnt = np.unique(key, return_counts=True)
    rows = np.bincount(ia).astype(float)
    cols = np.bincount(ib).astype(float)
    return cnt.astype(float), rows, cols, len(a)


def ari(a, b):
    n_ij, a_i, b_j, n = contingency(a, b)
    c2 = lambda x: x * (x - 1) / 2
    s_ij, s_a, s_b, s_n = c2(n_ij).sum(), c2(a_i).sum(), c2(b_j).sum(), c2(n)
    exp = s_a * s_b / s_n
    mx = (s_a + s_b) / 2
    return (s_ij - exp) / (mx - exp) if mx != exp else 1.0


def nmi(a, b):
    n_ij, a_i, b_j, n = contingency(a, b)
    ent = lambda x: -(x / n * np.log(x / n)).sum()
    ha, hb = ent(a_i), ent(b_j)
    # mutual information from the sparse contingency (recompute pairs)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    key = ia.astype(np.int64) * len(ub) + ib
    k, cnt = np.unique(key, return_counts=True)
    pi = a_i[k // len(ub)] / n
    pj = b_j[k % len(ub)] / n
    pij = cnt / n
    mi = (pij * np.log(pij / (pi * pj))).sum()
    return 2 * mi / (ha + hb) if ha + hb > 0 else 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True, help='NovInvenio clone (bin/ and lib/)')
    ap.add_argument('--pairs', action='append', required=True, help='NAME=pair_classification.tsv.zst')
    ap.add_argument('--resolutions', default='0.25,0.5,1,1.5,2,3,4,6,8,12,16')
    ap.add_argument('--seeds', type=int, default=10)
    ap.add_argument('--out', required=True)
    ap.add_argument('--labels_dir', help='write <dataset>.labels.npz (family names + labels[res, seed]) for AMI')
    args = ap.parse_args()
    sys.path.insert(0, args.code + '/bin')
    sys.path.insert(0, args.code + '/lib')
    import pangenome_detect_trans_modules as m

    res_list = [float(x) for x in args.resolutions.split(',')]
    out = open(args.out, 'w')
    out.write('dataset\tfamilies\tedges\tresolution\tmodules_mean\tmodules_min\tmodules_max\t'
              'modules_ge10_mean\tlargest_share_mean\tseed_ARI_mean\tseed_ARI_min\t'
              'seed_NMI_mean\tARI_vs_prev_res\tNMI_vs_prev_res\n')
    for spec in args.pairs:
        name, path = spec.split('=', 1)
        edges = m.load_trans_edges(path)
        g = m.build_graph(edges)
        names = g.vs['name']
        print(f'{name}: {g.vcount()} families, {g.ecount()} trans edges', flush=True)
        prev = None
        all_labs = []
        for r in res_list:
            labs = []
            for s in range(args.seeds):
                fm = m.detect_modules(g, resolution=r, seed=s)
                labs.append(np.array([fm[x] for x in names]))
            nmod = [len(np.unique(l)) for l in labs]
            ge10 = [int((np.bincount(l) >= 10).sum()) for l in labs]
            big = [np.bincount(l).max() / len(l) for l in labs]
            pa = [ari(a, b) for a, b in itertools.combinations(labs, 2)]
            pn = [nmi(a, b) for a, b in itertools.combinations(labs, 2)]
            vs_a = f'{ari(prev, labs[0]):.3f}' if prev is not None else ''
            vs_n = f'{nmi(prev, labs[0]):.3f}' if prev is not None else ''
            prev = labs[0]
            all_labs.append(np.stack(labs))
            out.write(f'{name}\t{g.vcount()}\t{g.ecount()}\t{r}\t{np.mean(nmod):.1f}\t{min(nmod)}\t{max(nmod)}\t'
                      f'{np.mean(ge10):.1f}\t{np.mean(big):.3f}\t{np.mean(pa):.3f}\t{min(pa):.3f}\t'
                      f'{np.mean(pn):.3f}\t{vs_a}\t{vs_n}\n')
            out.flush()
            print(f'  r={r}: modules {np.mean(nmod):.1f} (>=10: {np.mean(ge10):.1f}), largest {np.mean(big):.3f}, '
                  f'seed ARI {np.mean(pa):.3f} (min {min(pa):.3f}), NMI {np.mean(pn):.3f}', flush=True)
        if args.labels_dir:
            np.savez_compressed(f'{args.labels_dir}/{name}.labels.npz', names=np.array(names),
                                resolutions=np.array(res_list), labels=np.stack(all_labs))


if __name__ == '__main__':
    main()
