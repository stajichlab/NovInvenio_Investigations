#!/usr/bin/env python3.12
"""Lay out a splits network: SANS splits on a NeighborNet circular ordering.

1. Circular ordering: SplitsPy's NeighborNet cycle (Bryant & Moulton 2004) on a
   distance matrix. Only the ordering is used; SplitsPy's least-squares split
   weights are too slow in Python for 529 taxa.
2. Splits: a SANS splits file (weight, then taxa on one side). A split is drawn
   only if its taxa are contiguous on the cycle (circular split). The fraction of
   splits and of total weight that is kept is recorded: that is the part of the
   signal a planar network can show.
3. Layout: SplitsPy's outline algorithm (Huson et al. 2021).

Usage:
  build_network.py --splits dna.weakly.splits --dist jaccard.npy --labels strains.txt \
     [--weight sqrt] --title "..." --out net.json
  --dist may be a .npy square matrix (rows in --labels order) or a mash-style
  square TSV (header '#query', names ending .dna.fa are stripped).
"""
import argparse
import gzip
import json
import math

import numpy as np
from splitspy.nnet import nnet_cycle
from splitspy.outlines import outline_algo
from splitspy.splits.basic_split import Split


def read_dist(path, labels):
    if path.endswith('.npy'):
        return np.load(path)
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        head = [h.replace('.dna.fa', '') for h in fh.readline().rstrip('\n').split('\t')[1:]]
        rows = {}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            rows[f[0].replace('.dna.fa', '')] = [float(x) for x in f[1:]]
    col = {h: i for i, h in enumerate(head)}
    return np.array([[rows[a][col[b]] for b in labels] for a in labels])


def is_circular(side_pos, n):
    p = sorted(side_pos)
    gaps = sum(1 for a, b in zip(p, p[1:]) if b - a > 1) + (1 if p[0] + n - p[-1] > 1 else 0)
    return gaps <= 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--splits', required=True)
    ap.add_argument('--dist', required=True)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--weight', choices=['raw', 'sqrt', 'log'], default='sqrt',
                    help='edge length transform for display (raw weights are kept in the JSON)')
    ap.add_argument('--title', default='')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    labels = [l.strip() for l in open(args.labels) if l.strip()]
    n = len(labels)
    tax = {s: i + 1 for i, s in enumerate(labels)}     # SplitsPy taxa are 1-based
    D = read_dist(args.dist, labels)
    cycle = nnet_cycle.compute(labels, D.tolist())     # cycle[1..n]
    pos = {cycle[i]: i for i in range(1, n + 1)}

    raw = []
    for line in open(args.splits):
        f = line.rstrip('\n').split('\t')
        if len(f) < 2:
            continue
        side = [tax[s] for s in f[1:] if s in tax]
        if 0 < len(side) < n:
            raw.append((float(f[0]), side))
    total_w = sum(w for w, _ in raw)
    tf = {'raw': lambda w: w, 'sqrt': math.sqrt, 'log': lambda w: math.log1p(w)}[args.weight]

    kept, kept_w, info = [], 0.0, []
    all_t = set(range(1, n + 1))
    for w, side in raw:
        if not is_circular([pos[t] for t in side], n):
            continue
        kept.append(Split(set(side), all_t - set(side), tf(w)))
        info.append({'w': w, 'k': len(side), 'taxa': [t - 1 for t in side]})
        kept_w += w
    print(f'{len(kept)} of {len(raw)} splits circular on the NeighborNet ordering, '
          f'{kept_w / total_w:.3f} of total weight')

    graph, _ = outline_algo.compute(labels, cycle, kept, use_wts=True)
    nodes, nid = [], {}
    for v in graph.nodes():
        nid[v.id()] = len(nodes)
        nodes.append({'x': round(v.pos[0], 5), 'y': round(v.pos[1], 5),
                      'taxa': [tax[s] - 1 for s in v.label.split(',')] if v.label else []})
    edges = []
    for e in graph.edges():
        s = e.info if isinstance(e.info, int) else -1
        edges.append([nid[e.src().id()], nid[e.tar().id()], s])
    json.dump({
        'title': args.title,
        'n_splits_input': len(raw), 'n_splits_drawn': len(kept),
        'weight_frac_drawn': round(kept_w / total_w, 4) if total_w else 0,
        'weight_transform': args.weight,
        'cycle': [cycle[i] - 1 for i in range(1, n + 1)],
        'nodes': nodes, 'edges': edges, 'splits': info,
    }, open(args.out, 'w'), separators=(',', ':'))


if __name__ == '__main__':
    main()
