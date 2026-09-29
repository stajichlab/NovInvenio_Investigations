#!/usr/bin/env python3
"""Build the Coccidioides pangenome species tree for the Dollo step.

Source: existing 50-BUSCO-locus CDS IQ-TREE consensus (.contree, UFboot supports).
Steps: map tip names -> pangenome Short IDs; root on the Ci/Cp split; prune to pangenome
strains; graft strains absent from the CDS tree as zero-length sisters of their nearest
mash neighbour; add U. reesii as the root outgroup; write full + support-collapsed trees,
with and without U. reesii.
"""
import argparse, copy, csv, re, sys
from collections import Counter
from pathlib import Path
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

ap = argparse.ArgumentParser()
ap.add_argument('--cds_tree', required=True)
ap.add_argument('--mash_dist', required=True)
ap.add_argument('--config', required=True)
ap.add_argument('--outdir', required=True)
ap.add_argument('--collapse_below', type=float, default=70.0)
ap.add_argument('--outgroup_short', default='Uree')
ap.add_argument('--outgroup_branch', type=float, default=0.05)
a = ap.parse_args()
out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

cfg = list(csv.DictReader(open(a.config)))
species = {r['Short']: r['Species'] for r in cfg}
strain2short = {}
for r in cfg:
    strain2short[r['Strain']] = r['Short']; strain2short[r['Short']] = r['Short']
CI, CP = 'Coccidioides immitis', 'Coccidioides posadasii'
cocci = sorted(s for s in species if species[s] in (CI, CP))

# ---- mash matrix, keyed by Short
rows = [l.rstrip('\n').split('\t') for l in open(a.mash_dist)]
dna2short = {Path(r['DNA']).name: r['Short'] for r in cfg}
mnames = [dna2short[n] for n in rows[0][1:]]
D = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
midx = {s: i for i, s in enumerate(mnames)}
assert set(mnames) == set(species), 'mash matrix and config disagree'

# ---- name map
# Spelling variants accepted only if the mash data agree (checked below).
SPELLING = {'Coccidioides_posadasii_Coahuilla_2': 'Coahuila_2'}
t = Phylo.read(a.cds_tree, 'newick')
map_rows, drop = [], []
tipmap = {}
conflicts = []
for tip in t.get_terminals():
    n = tip.name
    m = re.match(r'^Coccidioides_(immitis|posadasii)_(.+)$', n)
    short, how = None, ''
    if n in SPELLING:
        short, how = SPELLING[n], 'spelling_variant(mash-checked)'
    elif m and m.group(2) in strain2short:
        short, how = strain2short[m.group(2)], 'exact(prefix stripped)'
    if short:
        if short in tipmap.values():
            sys.exit(f'duplicate mapping to {short}')
        tipmap[n] = short
        treesp = 'Coccidioides ' + m.group(1)
        if species[short] != treesp:
            how += f';SPECIES_LABEL_CONFLICT(cds_tree={treesp},config={species[short]})'
            conflicts.append((n, short))
    map_rows.append((n, short or '', how or 'UNMAPPED'))
    if not short:
        drop.append(n)
if len(t.get_terminals()) != sum(1 for n in [x.name for x in t.get_terminals()] if n.startswith(('Coccidioides_', 'Cimmitis', 'Cposadasii'))):
    sys.exit('non-Coccidioides tips present; rooting logic assumes none')

# ---- prune unmapped tips first (they are not pangenome strains)
for n in drop:
    t.prune(n)

# ---- root on the Ci/Cp split (the CDS tree has no non-Coccidioides tips)
ci_tips = [x for x in t.get_terminals() if x.name in tipmap and species[tipmap[x.name]] == CI]
print('species-label conflicts:', conflicts)
for n, s_ in conflicts:
    i = midx[s_]
    print('  mash median to config-Ci', np.median([D[i, midx[x]] for x in cocci if species[x]==CI and x!=s_]), 'to config-Cp', np.median([D[i, midx[x]] for x in cocci if species[x]==CP and x!=s_]))
# Rooting: first root on one Ci tip so the Cp set is a proper clade, then root on the Cp
# clade. Bio.Phylo leaves a trifurcating root here (Ci | Cp-part1 | Cp-part2); every
# non-Ci root child is pure Cp, so they are wrapped into one Cp clade. That is the
# same unrooted topology, rooted on the Ci/Cp branch.
t.root_with_outgroup(ci_tips[0])
cp_tips = [x for x in t.get_terminals() if species[tipmap[x.name]] == CP]
cp_mrca = t.common_ancestor(*cp_tips)
if len(cp_mrca.get_terminals()) != len(cp_tips):
    sys.exit('config-labelled Cp is not a clade in the CDS tree')
cp_support = cp_mrca.confidence
t.root_with_outgroup(cp_mrca)
kids = t.root.clades
ci_kids = [c for c in kids if all(species[tipmap[x.name]] == CI for x in c.get_terminals())]
cp_kids = [c for c in kids if all(species[tipmap[x.name]] == CP for x in c.get_terminals())]
if len(ci_kids) + len(cp_kids) != len(kids) or len(ci_kids) != 1:
    sys.exit('root children are not a clean Ci / Cp partition')
if len(cp_kids) > 1:
    t.root.clades = [ci_kids[0], Clade(clades=cp_kids, branch_length=0.0, confidence=cp_support)]
print('root children (tips):', [len(c.get_terminals()) for c in t.root.clades], 'Cp-node support', cp_support,
      'Ci-node support', t.root.clades[0].confidence)
for n, s_ in conflicts:
    tip = next(x for x in t.get_terminals() if x.name == n)
    par = t.get_path(tip)[-2]
    print('  conflict tip', n, 'parent clade (config species):', dict(Counter(species[tipmap[x.name]] for x in par.get_terminals())))
if len(t.root.clades) != 2:
    sys.exit('root not bifurcating after rooting on Ci')
for tip in t.get_terminals():
    tip.name = tipmap[tip.name]
in_tree = {x.name for x in t.get_terminals()}

# verify spelling variant with mash: Coahuila_2's mash nearest neighbours vs tree neighbours
def mash_nn(s, pool):
    i = midx[s]; best = min((D[i, midx[p]], p) for p in pool if p != s)
    return best
spell_check = []
for n, s in SPELLING.items():
    d, nn = mash_nn(s, in_tree)
    parent = t.get_path(s)[-2] if len(t.get_path(s)) > 1 else t.root
    sibs = {x.name for x in parent.get_terminals()}
    spell_check.append((n, s, nn, d, nn in sibs, len(sibs)))

# ---- grafts
missing = [s for s in cocci if s not in in_tree]
graft_rows = []
for s in missing:
    d, sister = mash_nn(s, in_tree - {a.outgroup_short})
    dall, nn_all = mash_nn(s, set(cocci))
    if species[sister] != species[s]:
        sys.exit(f'graft {s} nearest {sister} is other species')
    leaf = next(x for x in t.get_terminals() if x.name == sister)
    # leaf becomes an internal node holding the sister and the grafted strain at length 0.
    leaf.clades = [Clade(name=sister, branch_length=0.0), Clade(name=s, branch_length=0.0)]
    leaf.name = None; leaf.confidence = None
    graft_rows.append((s, species[s], sister, f'{d:.6g}', nn_all, f'{dall:.6g}'))
    # grafts attach only to original CDS-tree strains (in_tree not updated on purpose)

cocci_tree = t
def species_monophyly(tree):
    res = {}
    for sp_ in (CI, CP):
        tips = [x for x in tree.get_terminals() if species.get(x.name) == sp_]
        mrca = tree.common_ancestor(*tips)
        res[sp_] = Counter(species[x.name] for x in mrca.get_terminals())
    return res

def with_outgroup(tree):
    tr = copy.deepcopy(tree)
    old = tr.root; old.branch_length = a.outgroup_branch
    tr.root = Clade(clades=[old, Clade(name=a.outgroup_short, branch_length=a.outgroup_branch)])
    tr.rooted = True
    return tr

def collapse(tree):
    tr = copy.deepcopy(tree)
    n0 = len(tr.get_nonterminals())
    # never collapse the root or its direct children (the Ci/Cp split and the Cocci node)
    protect = {id(tr.root)} | {id(c) for c in tr.root.clades}
    targets = [c for c in tr.get_nonterminals()
               if id(c) not in protect and c.confidence is not None and c.confidence < a.collapse_below]
    for c in targets:
        tr.collapse(c)
    return tr, n0, len(tr.get_nonterminals()), len(targets)

outputs = {}
genus = with_outgroup(cocci_tree)
for label, tree in [('genus_with_Uree', genus), ('coccidioides_only', cocci_tree)]:
    full = copy.deepcopy(tree)
    col, n0, n1, nc = collapse(tree)
    for kind, tr in [('full', full), (f'collapsed_ufboot{int(a.collapse_below)}', col)]:
        p = out / f'coccidioides_{label}.{kind}.nwk'
        tr.rooted = True
        Phylo.write(tr, p, 'newick')
        outputs[p.name] = (tr, n0 if kind == 'full' else n1)
    print(f'{label}: internal nodes {n0} -> {n1} after collapsing {nc} nodes with support < {a.collapse_below}')

with open(out / 'tip_name_map.tsv', 'w') as f:
    f.write('cds_tree_tip\tpangenome_short\tmapping\n')
    for r in map_rows: f.write('\t'.join(r) + '\n')
with open(out / 'grafts.tsv', 'w') as f:
    f.write('grafted_strain\tspecies\tsister_in_cds_tree\tmash_dist_to_sister\tnearest_any_pangenome_strain\tmash_dist_nearest_any\n')
    for r in graft_rows: f.write('\t'.join(r) + '\n')
print('CDS tips', len(map_rows), 'mapped', len(tipmap), 'unmapped', drop)
print('grafted', len(graft_rows))
print('spelling check (tree tip, short, mash NN among tree strains, dist, NN is in same parent clade, parent size):', spell_check)
for n, (tr, k) in outputs.items():
    print(n, 'tips', len(tr.get_terminals()), 'root children', len(tr.root.clades), 'monophyly', {k2: dict(v) for k2, v in species_monophyly(tr).items()})
