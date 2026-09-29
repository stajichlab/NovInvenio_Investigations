#!/usr/bin/env python3.12
"""Build the data file for the pangenome tree/heatmap/network prototype.

Writes one JSON (gzip-free, embedded in the page later) with:
  strains      name, species
  trees        per tree: nodes (parent, x, y, n_tips, support, preorder), tip order
  families     id, pfam, packed presence bits over all strains (base64)
  per tree     family order, best-clade node, jaccard, carrier side, fitch changes
  networks     per distance: NeighborNet outline nodes/edges (from --network-json)

Usage:
  build_viz_data.py --matrix presence_matrix.rescued.tsv --config config.csv \
     --pfam pfam.domtblout --tree snp=... --tree busco=... \
     --rename Ref_CimmRS=RS --rename Coahuilla_2=Coahuila_2 \
     [--network gene_content=net_gc.json --network mash=net_mash.json] --out viz_data.json
"""
import argparse
import base64
import csv
import json
import sys

import numpy as np
from Bio import Phylo

sys.path.insert(0, __file__.rsplit('/', 2)[0])      # analysis/sans, for compare_trees
from compare_trees import fitch_lengths          # noqa: E402

PRESENT = {'present', 'genome_only'}


def load_tree(path, rename, strains_set, species):
    tree = Phylo.read(path, 'newick')
    for t in tree.get_terminals():
        t.name = rename.get(t.name, t.name)
    for t in list(tree.get_terminals()):
        if t.name not in strains_set:
            tree.prune(t)
    # root between the species when the tree is unrooted (SNP trees)
    cp = [t for t in tree.get_terminals() if species.get(t.name) == 'Cp']
    ci = [t for t in tree.get_terminals() if species.get(t.name) == 'Ci']
    if len(tree.root.clades) != 2 and cp and ci:
        try:
            tree.root_with_outgroup(cp)
        except ValueError:
            tree.root_with_outgroup(ci)
    tree.ladderize()
    return tree


def layout(tree):
    """Cladogram with aligned tips: x = depth from root in edges scaled to tips."""
    nodes, index = [], {}
    tips = []

    def height(c):
        return 0 if c.is_terminal() else 1 + max(height(k) for k in c.clades)

    sys.setrecursionlimit(100000)
    H = height(tree.root)

    def walk(c, parent, depth):
        i = len(nodes)
        index[id(c)] = i
        nodes.append({'p': parent, 'd': depth, 'sup': c.confidence, 'name': c.name if c.is_terminal() else None})
        if c.is_terminal():
            nodes[i]['y'] = len(tips)
            tips.append(c.name)
            nodes[i]['n'] = 1
            nodes[i]['lo'] = nodes[i]['hi'] = nodes[i]['y']
        else:
            kids = [walk(k, i, depth + 1) for k in c.clades]
            nodes[i]['y'] = sum(nodes[k]['y'] for k in kids) / len(kids)
            nodes[i]['n'] = sum(nodes[k]['n'] for k in kids)
            nodes[i]['lo'] = min(nodes[k]['lo'] for k in kids)
            nodes[i]['hi'] = max(nodes[k]['hi'] for k in kids)
        return i

    walk(tree.root, -1, 0)
    for nd in nodes:          # tips aligned at x = 1
        nd['x'] = 1.0 if nd['name'] else nd['d'] / H
    return nodes, tips


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--pfam')
    ap.add_argument('--pfam-evalue', type=float, default=1e-5)
    ap.add_argument('--tree', action='append', required=True, help='NAME=path')
    ap.add_argument('--label', action='append', default=[], help='NAME=display label')
    ap.add_argument('--rename', action='append', default=[])
    ap.add_argument('--network', action='append', default=[], help='NAME=json from build_network.py')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    rename = dict(r.split('=', 1) for r in args.rename)
    labels = dict(r.split('=', 1) for r in args.label)

    species = {}
    for r in csv.DictReader(open(args.config)):
        sp = r['Species']
        species[r['Short']] = 'Ci' if 'immitis' in sp else 'Cp' if 'posadasii' in sp else sp

    with open(args.matrix) as fh:
        strains = fh.readline().rstrip('\n').split('\t')[1:]
        fam_ids, rows = [], []
        for line in fh:
            f = line.rstrip('\n').split('\t')
            fam_ids.append(f[0])
            rows.append([v in PRESENT for v in f[1:]])
    A = np.array(rows, dtype=bool)                 # families x strains
    n_all = A.shape[1]
    k_all = A.sum(1)
    keep_fam = np.minimum(k_all, n_all - k_all) >= 2
    print(f'{A.shape[0]} families, {int(keep_fam.sum())} with minor side >= 2 over all strains')
    fam_idx = np.flatnonzero(keep_fam)
    A = A[fam_idx]
    fam_ids = [fam_ids[i] for i in fam_idx]
    sidx = {s: i for i, s in enumerate(strains)}

    pfam = {}
    if args.pfam:
        best = {}
        for line in open(args.pfam):
            if line.startswith('#'):
                continue
            p = line.split()
            name, q, iev = p[0], p[3], float(p[12])
            if iev <= args.pfam_evalue:
                best.setdefault(q, {})
                best[q][name] = min(iev, best[q].get(name, 1.0))
        for q, d in best.items():
            pfam[q] = [n for n, _ in sorted(d.items(), key=lambda kv: kv[1])][:4]

    out = {
        'strains': [{'name': s, 'sp': species.get(s, '?')} for s in strains],
        'families': {
            'id': fam_ids,
            'pfam': [pfam.get(f, []) for f in fam_ids],
            'count': A.sum(1).tolist(),
            'bits': base64.b64encode(np.packbits(A, axis=1).tobytes()).decode(),
            'bytes_per_family': int(np.packbits(A[:1], axis=1).shape[1]),
        },
        'trees': {},
        'networks': {},
    }

    for spec in args.tree:
        name, path = spec.split('=', 1)
        tree = load_tree(path, rename, set(strains), species)
        nodes, tips = layout(tree)
        tip_cols = np.array([sidx[t] for t in tips])
        sub = A[:, tip_cols]                        # families x tips (tree order)
        n = len(tips)
        k = sub.sum(1)
        informative = np.minimum(k, n - k) >= 2
        loss_side = k > n / 2                       # minor side is the absent set
        side = np.where(loss_side[:, None], ~sub, sub).astype(np.float32)
        # clade membership over tips, for internal nodes and tips
        C = np.zeros((n, len(nodes)), dtype=np.float32)
        for j, nd in enumerate(nodes):
            C[nd['lo']:nd['hi'] + 1, j] = 1.0
        inter = side @ C                            # families x nodes
        sizes = side.sum(1)[:, None] + C.sum(0)[None, :] - inter
        jac = np.where(sizes > 0, inter / np.maximum(sizes, 1), 0)
        jac[:, 0] = np.where(sub.all(1) | (~sub).all(1), 1, 0)  # root clade is uninformative
        best = jac.argmax(1)
        bestj = jac[np.arange(len(best)), best]
        pres = {t: sub[:, i] for i, t in enumerate(tips)}
        changes = fitch_lengths(tree, set(tips), pres)
        fams = np.flatnonzero(informative)
        order = sorted(fams, key=lambda f: (int(best[f]), -float(bestj[f]), int(k[f])))
        out['trees'][name] = {
            'label': labels.get(name, name),
            'path': path,
            'tips': tips,
            'tip_strain_index': tip_cols.tolist(),
            'nodes': [{kk: nd[kk] for kk in ('p', 'x', 'y', 'n', 'sup', 'lo', 'hi')} for nd in nodes],
            'order': [int(f) for f in order],
            'best': [int(best[f]) for f in order],
            'jac': [round(float(bestj[f]), 3) for f in order],
            'loss': [bool(loss_side[f]) for f in order],
            'changes': [int(changes[f]) for f in order],
            'k': [int(k[f]) for f in order],
        }
        print(f'{name}: {n} tips, {len(nodes)} nodes, {len(order)} informative families, '
              f'{int((changes[fams] == 1).sum())} with one change')

    for spec in args.network:
        name, path = spec.split('=', 1)
        out['networks'][name] = json.load(open(path))

    with open(args.out, 'w') as fh:
        json.dump(out, fh, separators=(',', ':'))


if __name__ == '__main__':
    main()
