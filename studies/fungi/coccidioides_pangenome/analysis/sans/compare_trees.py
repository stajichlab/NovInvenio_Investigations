#!/usr/bin/env python3.12
"""Compare split sets of several trees on their shared strains.

For every (query, reference) pair, all trees are restricted to the strains
that are in both. A split is a non-trivial bipartition (both sides >= 2).

Reported per pair:
  q_splits, r_splits  non-trivial splits of each restricted tree
  shared              splits in both
  precision           shared / q_splits
  recall              shared / r_splits
  recall_sup95        recall over reference splits with support >= 95
                      (only if the reference has node support values)
  rf_norm             (q + r - 2*shared) / (q + r)

Optional --families: a presence matrix. Each family's carrier set,
restricted to the shared strains of the reference tree, is classed as
  clade        equals one side of a reference split (or one strain side = tip)
  compatible   does not conflict with any reference split, but is not a clade
  conflict     crosses at least one reference split
This counts how many gene families fit the tree as one gain or one loss.

Usage:
  compare_trees.py --tree NAME=path.nwk ... --rename OLD=NEW ... \
      [--families presence_matrix.tsv --families-ref NAME]
"""
import argparse
import itertools
from Bio import Phylo


def load(path, rename):
    tree = Phylo.read(path, 'newick')
    for tip in tree.get_terminals():
        tip.name = rename.get(tip.name, tip.name)
    return tree


def splits(tree, keep):
    """Return {frozenset(side): support} for non-trivial splits on `keep`."""
    keep = frozenset(keep)
    n = len(keep)
    anchor = min(keep)          # canonical side: the side without `anchor`
    out = {}
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            continue
        side = frozenset(t.name for t in clade.get_terminals()) & keep
        if anchor in side:
            side = keep - side
        if 2 <= len(side) <= n - 2:
            sup = clade.confidence
            if side not in out or (sup is not None and (out[side] is None or sup > out[side])):
                out[side] = sup
    return out


def fitch_lengths(tree, keep, pres):
    """Unweighted small-parsimony length per family (Hartigan's rule at polytomies).

    pres: {strain: bool numpy array over families}. Returns an int array: the
    minimum number of presence/absence changes of each family on `tree`,
    restricted to `keep`. Unrooted length is the same as rooted for this rule.
    """
    import numpy as np
    n_fam = len(next(iter(pres.values())))
    cost = np.zeros(n_fam, dtype=np.int32)

    def post(clade):
        # returns (has0, has1) boolean arrays = the Fitch state set of this node,
        # or None if no kept strain below
        if clade.is_terminal():
            if clade.name not in keep:
                return None
            p = pres[clade.name]
            return (~p, p)
        kids = [s for s in (post(c) for c in clade.clades) if s is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        c0 = sum(k[0].astype(np.int32) for k in kids)
        c1 = sum(k[1].astype(np.int32) for k in kids)
        m = np.maximum(c0, c1)
        cost[:] += len(kids) - m
        return (c0 == m, c1 == m)

    import sys
    sys.setrecursionlimit(100000)
    post(tree.root)
    return cost


def parsimony_report(trees, tips, matrix, present):
    import numpy as np
    common = frozenset.intersection(*(frozenset(t) for t in tips.values()))
    with open(matrix) as fh:
        header = fh.readline().rstrip('\n').split('\t')[1:]
        rows = [line.rstrip('\n').split('\t')[1:] for line in fh]
    cols = {s: i for i, s in enumerate(header) if s in common}
    arr = np.array([[r[i] in present for i in cols.values()] for r in rows], dtype=bool)
    pres = {s: arr[:, j] for j, s in enumerate(cols)}
    minor = np.minimum(arr.sum(1), len(cols) - arr.sum(1))
    informative = minor >= 2
    bands = [(2, 2), (3, 5), (6, 20), (21, 100), (101, 10**9)]
    print(f'# parsimony on {len(cols)} strains common to all trees; '
          f'{int(informative.sum())} informative families (minor side >= 2)')
    print('tree\ttotal_changes\tmean_per_family\tfrac_one_change\t'
          + '\t'.join(f'mean_minor{a}-{b if b < 10**9 else "max"}' for a, b in bands))
    for name, tree in trees.items():
        c = fitch_lengths(tree, cols, pres)[informative]
        mi = minor[informative]
        band = '\t'.join(f'{c[(mi >= a) & (mi <= b)].mean():.2f}' for a, b in bands)
        print(f'{name}\t{int(c.sum())}\t{c.mean():.2f}\t{(c == 1).mean():.3f}\t{band}')
    counts = '\t'.join(str(int(((minor >= a) & (minor <= b)).sum())) for a, b in bands)
    print(f'n_families_per_band\t\t\t\t{counts}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tree', action='append', required=True, help='NAME=path')
    ap.add_argument('--rename', action='append', default=[], help='OLD=NEW tip name')
    ap.add_argument('--families')
    ap.add_argument('--families-ref', action='append', default=[])
    ap.add_argument('--present', default='present,genome_only')
    ap.add_argument('--parsimony', action='store_true',
                    help='with --families: parsimony length of every family on every tree')
    ap.add_argument('--shuffle', action='append', default=[],
                    help='add NAME_shuffled: NAME with tip labels permuted (seed 1), a no-signal control')
    args = ap.parse_args()

    rename = dict(r.split('=', 1) for r in args.rename)
    trees = {}
    for spec in args.tree:
        name, path = spec.split('=', 1)
        trees[name] = load(path, rename)
    for name in args.shuffle:
        import copy
        import random
        t = copy.deepcopy(trees[name])
        leaves = t.get_terminals()
        labels = [x.name for x in leaves]
        random.Random(1).shuffle(labels)
        for x, lab in zip(leaves, labels):
            x.name = lab
        trees[name + '_shuffled'] = t
    tips = {k: {t.name for t in v.get_terminals()} for k, v in trees.items()}
    for k, v in tips.items():
        print(f'# {k}: {len(v)} tips')

    print('query\treference\tn_shared_tips\tq_splits\tr_splits\tshared\tprecision\trecall\trecall_sup95\trf_norm')
    for q, r in itertools.permutations(trees, 2):
        keep = tips[q] & tips[r]
        qs, rs = splits(trees[q], keep), splits(trees[r], keep)
        sh = set(qs) & set(rs)
        sup = [s for s, v in rs.items() if v is not None and v >= 95]
        r95 = (f'{sum(s in qs for s in sup) / len(sup):.3f} ({len(sup)})' if sup else '-')
        rf = (len(qs) + len(rs) - 2 * len(sh)) / (len(qs) + len(rs))
        print(f'{q}\t{r}\t{len(keep)}\t{len(qs)}\t{len(rs)}\t{len(sh)}\t'
              f'{len(sh) / len(qs):.3f}\t{len(sh) / len(rs):.3f}\t{r95}\t{rf:.3f}')

    if args.families and args.parsimony:
        parsimony_report(trees, tips, args.families, set(args.present.split(',')))

    if args.families:
        present = set(args.present.split(','))
        for ref in args.families_ref:
            keep = frozenset(tips[ref])
            anchor = min(keep)
            rs = list(splits(trees[ref], keep))
            rset = set(rs)
            counts = {'clade': 0, 'compatible': 0, 'conflict': 0, 'trivial': 0}
            with open(args.families) as fh:
                header = fh.readline().rstrip('\n').split('\t')[1:]
                idx = [i for i, s in enumerate(header) if s in keep]
                names = [header[i] for i in idx]
                for line in fh:
                    cells = line.rstrip('\n').split('\t')[1:]
                    car = frozenset(names[j] for j, i in enumerate(idx) if cells[i] in present)
                    side = keep - car if anchor in car else car
                    if len(side) <= 1 or len(side) >= len(keep) - 1:
                        counts['trivial'] += 1
                    elif side in rset:
                        counts['clade'] += 1
                    elif all(not (side & s and side - s and s - side and (keep - side) & (keep - s))
                             for s in rs):
                        counts['compatible'] += 1
                    else:
                        counts['conflict'] += 1
            nt = sum(counts.values()) - counts['trivial']
            print(f'# families vs {ref} ({len(keep)} strains): ' + ', '.join(
                f'{k} {v}' for k, v in counts.items())
                + f'; of {nt} non-trivial: clade {counts["clade"] / nt:.3f}, '
                  f'compatible {counts["compatible"] / nt:.3f}, conflict {counts["conflict"] / nt:.3f}')


if __name__ == '__main__':
    main()
