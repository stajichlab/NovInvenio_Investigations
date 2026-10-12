import csv, sys, math, collections
import numpy as np
from scipy.stats import fisher_exact
P='/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/FSSC_ambrosia/results/mmseqs_norescue/output/pangenome/'
inv={r['Short']:r for r in csv.DictReader(open(P+'strain_inventory.tsv'),delimiter='\t')}
ss={r['Short']:r for r in csv.DictReader(open(P+'samplesheet.with_clades.csv'))}
reps=[s for s,r in ss.items() if r['GROUP']=='IN' and inv.get(s,{}).get('is_representative')=='1']
print('ingroup representatives:', len(reps))
freq={r['family']:r for r in csv.DictReader(open(P+'frequency_table.tsv'),delimiter='\t')}
fams=[]; rows=[]
with open(P+'presence_matrix.tsv') as fh:
    hdr=fh.readline().rstrip('\n').split('\t'); idx=[hdr.index(s) for s in reps]
    for line in fh:
        f=line.rstrip('\n').split('\t')
        v=np.array([f[i] in ('present','genome_only') for i in idx])
        k=v.sum()
        if k>=5 and freq.get(f[0],{}).get('bin') in ('shell','cloud','soft_core') and k < len(reps):
            fams.append(f[0]); rows.append(v)
X=np.array(rows,dtype=np.uint8); F,n=X.shape
print('eligible families (5 <= k < n, shell/cloud/soft_core):',F, 'pairs', F*(F-1)//2)
ks=X.sum(1)
print('carrier count distribution (k):', dict(sorted(collections.Counter(ks.tolist()).items())))
# species groups from species name
sp=[ss[s]['Species'] for s in reps]
spc=collections.Counter(sp); print('species among reps:', len(spc), dict(spc.most_common(6)))
spi=np.array([sorted(spc).index(x) for x in sp])
nsp=np.array([len(set(spi[r.astype(bool)])) for r in X])
print('species spanned per eligible family:', dict(sorted(collections.Counter(nsp.tolist()).items())))
# p lookup over (ka,kb,both)
cache={}
def pv(ka,kb,b):
    key=(ka,kb,b)
    if key not in cache:
        cache[key]=fisher_exact([[b,ka-b],[kb-b,n-ka-kb+b]],alternative='greater')[1]
    return cache[key]
m=F*(F-1)//2
# histogram of -log10 p in bins
edges=np.arange(0,12.5,0.5); hist=np.zeros(len(edges),dtype=np.int64)
pmin_all=1.0
Xi=X.astype(np.int32)
for i0 in range(0,F,500):
    B=Xi[i0:i0+500]@Xi.T
    for a in range(B.shape[0]):
        i=i0+a; js=np.arange(i+1,F)
        if not len(js): continue
        b=B[a,i+1:]; ka=ks[i]; kb=ks[i+1:]
        keys=list(zip([int(ka)]*len(js), kb.tolist(), b.tolist()))
        p=np.array([pv(*k) for k in keys])
        lp=np.clip(-np.log10(np.maximum(p,1e-300)),0,12.25)
        np.add.at(hist,(lp/0.5).astype(int),1)
        pmin_all=min(pmin_all,p.min())
print('smallest pair p:', pmin_all)
cum=np.cumsum(hist[::-1])[::-1]
print('pairs with p <= 10^-x  (and BH needs N >= p*m/0.05):')
for e,c in zip(edges,cum):
    if e>=4 and c>0:
        need=(10**-e)*m/0.05
        print(f'  x={e:4.1f}: {c:>12,}   need {need:>12,.0f}   {"PASS" if c>=need else ""}')
