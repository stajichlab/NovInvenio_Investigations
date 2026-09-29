import sys, random, time
sys.dont_write_bytecode = True
sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
from Bio import Phylo
from ancestral_states import dollo_polarize
from collections import Counter
import csv
T = sys.argv[1] if len(sys.argv) > 1 else '.'
R = '/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/results'
species = {r['Short']: r['Species'] for r in csv.DictReader(open('/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/config_genus_vs_ureesii.csv'))}
runs = {'rescue_structural_genus_vs_ureesii': 'coccidioides_genus_with_Uree',
        'rescue_freqpol_immitis_in_posadasii_out': 'coccidioides_coccidioides_only',
        'rescue_freqpol_posadasii_in_immitis_out': 'coccidioides_coccidioides_only'}
ok = True
for run, stem in runs.items():
    for kind in ('full', 'collapsed_ufboot70'):
        tr = Phylo.read(f'{T}/{stem}.{kind}.nwk', 'newick')   # same parser as bin/pangenome_cooccurrence.py:533
        tips = [x.name for x in tr.get_terminals()]
        dup = [k for k, v in Counter(tips).items() if v > 1]
        for m in ('presence_matrix.tsv', 'presence_matrix.rescued.tsv'):
            cols = set(open(f'{R}/{run}/output/pangenome/{m}').readline().rstrip('\n').split('\t')[1:])
            miss, extra = cols - set(tips), set(tips) - cols
            good = not miss and not extra and not dup and len(tr.root.clades) == 2
            ok &= good
            mono = {}
            for sp_ in ('Coccidioides immitis', 'Coccidioides posadasii'):
                ts = [x for x in tr.get_terminals() if species[x.name] == sp_]
                mono[sp_[13:]] = len(tr.common_ancestor(*ts).get_terminals()) == len(ts)
            ok &= all(mono.values())
            print(f'{run} {m} vs {stem}.{kind}: tips={len(tips)} matrix={len(cols)} missing={len(miss)} extra={len(extra)} dups={len(dup)} root_children={len(tr.root.clades)} monophyletic={mono}')
print('ALL CHECKS PASS' if ok else 'CHECK FAILURE')
# presence = present or genome_only, as lib/pangenome_matrix.py:113-114 is_present()
# Dollo smoke test on real families: loss-event counts, full vs collapsed tree
full = Phylo.read(f'{T}/coccidioides_genus_with_Uree.full.nwk', 'newick')
col = Phylo.read(f'{T}/coccidioides_genus_with_Uree.collapsed_ufboot70.nwk', 'newick')
f = open(f'{R}/rescue_structural_genus_vs_ureesii/output/pangenome/presence_matrix.rescued.tsv')
hdr = f.readline().rstrip('\n').split('\t')[1:]
rows = [l.rstrip('\n').split('\t') for l in f]
random.seed(1)
acc = [r for r in rows if 0 < sum(v in ('present', 'genome_only') for v in r[1:]) < len(hdr)]
print('families total', len(rows), 'accessory (not in all 530)', len(acc))
sample = random.sample(acc, min(1000, len(acc)))
t0 = time.time(); diffs = []; errs = 0
for r in sample:
    carriers = {h for h, v in zip(hdr, r[1:]) if v in ('present', 'genome_only')}
    try:
        a = dollo_polarize(full, carriers).n_loss_events; b = dollo_polarize(col, carriers).n_loss_events
        diffs.append((a, b))
    except Exception as e:
        errs += 1; print('ERR', e)
dt = time.time() - t0
more = sum(b > a for a, b in diffs); less = sum(b < a for a, b in diffs); same = sum(a == b for a, b in diffs)
import statistics as st
print(f'dollo on {len(diffs)} random accessory families x 2 trees: {dt:.1f}s, errors={errs}')
print(f'n_loss_events collapsed vs full: higher={more} lower={less} equal={same}; median full={st.median(a for a,_ in diffs)} median collapsed={st.median(b for _,b in diffs)}; sum full={sum(a for a,_ in diffs)} sum collapsed={sum(b for _,b in diffs)}')
