#!/usr/bin/env python3
"""Near-miss table for a pairwise (diamond/blastp) NovInvenio novelty run.

A novelty candidate passes two rules (bin/build_presence_matrix.py in nf_NovInvenio):
  - present in >= --ingroup fraction of IN proteomes (the protein's own proteome counts)
  - present in 0 OUT proteomes
A "near miss" passes the same ingroup rule but is present in 1..--max-out OUT
proteomes. This script lists those proteins. It does not change any call.

"Present" means exactly what the pipeline's presence_matrix.tsv says: a hit with
E < 1e-5 that survived the paralog-competition filter (and the coverage floor, if
the run used one). This script reads that matrix; it does not re-score hits.

Inputs (all under --results unless stated):
  presence_matrix.tsv                    0/1 presence per (IN protein, proteome)
  presence_matrix.evalues.tsv            best kept-hit E-value per cell
  presence_matrix.targets.tsv            best kept-hit target ID per cell
  presence_matrix.other_evidence.tsv.gz  per (protein, OUT proteome) cell: n_hits,
                                         max_qcov, status (kept / paralog_filtered)
  candidates.txt                         cross-check of the candidate count
  presence_matrix.uniprot.tsv            optional UniProt/Pfam annotation columns
  search_cache/<IN>_vs_<OUT>.diamond.tsv.gz
                                         optional raw diamond rows (outfmt 6 qseqid sseqid
                                         evalue bitscore length pident qcovhsp scovhsp qlen
                                         slen) for percent identity and qcov of the best
                                         kept hit
  --config CSV                           GROUP, Short, Protein and a lineage column
                                         (TaxonGroup by default)
  --proteins-dir                         protein FASTAs named as in config Protein column
                                         (default: <config dir>/data_dir/pep); used only
                                         for family clustering

Families: near-miss proteins are clustered with mmseqs easy-cluster using the same
parameters as the pipeline's MMSEQS_CLUSTER on candidates (--min-seq-id 0.3 -c 0.8
--cov-mode 0). Family IDs are the mmseqs representative IDs. They are NOT the same
clusters as the candidates' family_id in novelties.<Short>.tsv (that clustering used
candidates only). --family-method none sets family_id = protein_id.

Outputs (prefix = --out-prefix):
  <prefix>.tsv.gz                  protein-level near-miss table
  <prefix>_families.tsv            one row per near-miss family
  <prefix>_pfam_vs_candidates.tsv  Pfam counts, near misses vs candidates (plain counts)
  <prefix>_out_species.tsv         per OUT species: proteins/families it alone removes, etc.
  <prefix>_summary.txt             summary numbers
"""
import argparse
import csv
import gzip
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ANNOT_COLS = ['uniprot_accession', 'uniprot_description', 'uniprot_pfam_ids',
              'uniprot_pfam_names', 'uniprot_interpro_ids', 'uniprot_ec_numbers',
              'uniprot_go_ids']
DIAMOND_COLS = ['qseqid', 'sseqid', 'evalue', 'bitscore', 'length', 'pident',
                'qcovhsp', 'scovhsp', 'qlen', 'slen']


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def load_config(path, lineage_col):
    cfg = pd.read_csv(path, dtype=str).fillna('')
    cfg.columns = [c.strip() for c in cfg.columns]
    for need in ('GROUP', 'Short'):
        if need not in cfg.columns:
            sys.exit(f'ERROR: config {path} has no {need} column')
    if lineage_col not in cfg.columns:
        alt = [c for c in ('TaxonGroup', 'Group', 'Lineage', 'Taxon') if c in cfg.columns]
        if not alt:
            sys.exit(f'ERROR: config has no lineage column ({lineage_col}); pass --lineage-col')
        log(f'lineage column {lineage_col} not found; using {alt[0]}')
        lineage_col = alt[0]
    cfg['GROUP'] = cfg['GROUP'].str.strip().str.upper()
    cfg['Short'] = cfg['Short'].str.strip()
    in_ids = cfg.loc[cfg['GROUP'] == 'IN', 'Short'].tolist()
    out_ids = cfg.loc[cfg['GROUP'] == 'OUT', 'Short'].tolist()
    lineage = dict(zip(cfg['Short'], cfg[lineage_col].str.strip()))
    protein_file = dict(zip(cfg['Short'], cfg.get('Protein', pd.Series([''] * len(cfg)))))
    return in_ids, out_ids, lineage, protein_file


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 'NA'
    if isinstance(x, float):
        return f'{x:.{nd}g}'
    return str(x)


def read_raw_hits(results, in_short, out_short, keep_pairs):
    """Best raw diamond row per (query, target) in keep_pairs for one IN->OUT file."""
    path = Path(results) / 'search_cache' / f'{in_short}_vs_{out_short}.diamond.tsv.gz'
    if not path.exists():
        return None
    best = {}
    with open_text(path) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            key = (c[0], c[1])
            if key not in keep_pairs:
                continue
            ev = float(c[2])
            if key not in best or ev < best[key][0]:
                rec = dict(zip(DIAMOND_COLS, c))
                best[key] = (ev, rec)
    return {k: v[1] for k, v in best.items()}


def read_fasta_subset(path, keep):
    seqs, cur = {}, None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                pid = line[1:].split()[0]
                cur = pid if pid in keep else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def mmseqs_families(proteins, source_of, protein_file, proteins_dir, mmseqs, threads):
    """protein_id -> representative id, from mmseqs easy-cluster on the near-miss set."""
    by_src = defaultdict(set)
    for p in proteins:
        by_src[source_of[p]].add(p)
    seqs = {}
    for src, ids in by_src.items():
        fa = Path(proteins_dir) / protein_file[src]
        if not fa.exists():
            sys.exit(f'ERROR: protein FASTA not found for {src}: {fa}')
        got = read_fasta_subset(fa, ids)
        missing = ids - set(got)
        if missing:
            sys.exit(f'ERROR: {len(missing)} near-miss IDs not found in {fa} '
                     f'(e.g. {sorted(missing)[:3]})')
        seqs.update(got)
    base = os.environ.get('SCRATCH') or os.environ.get('TMPDIR') or None
    tmp = tempfile.mkdtemp(prefix='near_miss_mmseqs_', dir=base)
    try:
        fa = Path(tmp) / 'near_miss.fa'
        with open(fa, 'w') as fh:
            # Numeric headers: mmseqs rewrites UniProt-style 'tr|ACC|NAME' IDs in
            # *_cluster.tsv, so map back through an index instead.
            for i, p in enumerate(proteins):
                fh.write(f'>nm{i}\n{seqs[p]}\n')
        idx = {f'nm{i}': p for i, p in enumerate(proteins)}
        cmd = [mmseqs, 'easy-cluster', str(fa), str(Path(tmp) / 'nm'), str(Path(tmp) / 'work'),
               '--threads', str(threads), '--min-seq-id', '0.3', '-c', '0.8',
               '--cov-mode', '0', '-v', '1']
        log('running: ' + ' '.join(cmd))
        subprocess.run(cmd, check=True)
        fam = {}
        with open(Path(tmp) / 'nm_cluster.tsv') as fh:
            for line in fh:
                rep, mem = line.rstrip('\n').split('\t')[:2]
                fam[idx[mem]] = idx[rep]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    missing = set(proteins) - set(fam)
    if missing:
        # mmseqs can rewrite some header styles; fall back to self for any it lost.
        log(f'WARNING: {len(missing)} proteins not in mmseqs output; family_id = protein_id')
        for p in missing:
            fam[p] = p
    return fam


def split_ids(s):
    if not isinstance(s, str) or not s:
        return []
    return [x for x in s.replace(';', '|').replace(',', '|').split('|') if x]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True, help='pipeline results/<study> directory')
    ap.add_argument('--config', required=True, help='study config.csv (GROUP, Short, lineage)')
    ap.add_argument('--out-prefix', required=True, help='output path prefix')
    ap.add_argument('--max-out', type=int, default=2,
                    help='max number of OUT proteomes with a kept hit (default 2)')
    ap.add_argument('--ingroup', type=float, default=0.75,
                    help='min fraction of IN proteomes present, own proteome included '
                         '(default 0.75, as the pipeline)')
    ap.add_argument('--lineage-col', default='TaxonGroup')
    ap.add_argument('--proteins-dir', default=None,
                    help='dir of protein FASTAs (default <config dir>/data_dir/pep)')
    ap.add_argument('--family-method', choices=['mmseqs', 'none'], default='mmseqs')
    ap.add_argument('--mmseqs', default=shutil.which('mmseqs') or 'mmseqs')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_ON_NODE', 2)))
    ap.add_argument('--no-raw-hits', action='store_true',
                    help='skip search_cache diamond files (no pident / per-hit qcov)')
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--min-pfam-nm', type=int, default=5,
                    help='min near-miss proteins for a Pfam to be ranked by ratio (default 5)')
    args = ap.parse_args()
    t0 = time.time()
    res = Path(args.results)

    in_ids, out_ids, lineage, protein_file = load_config(args.config, args.lineage_col)
    log(f'IN ({len(in_ids)}): {",".join(in_ids)}')
    log(f'OUT ({len(out_ids)}): {",".join(out_ids)}')

    pm = pd.read_csv(res / 'presence_matrix.tsv', sep='\t', dtype={'protein_id': str,
                                                                   'source_proteome': str})
    missing = [c for c in in_ids + out_ids if c not in pm.columns]
    if missing:
        sys.exit(f'ERROR: config species not in presence_matrix.tsv: {missing}')
    pm = pm[pm['source_proteome'].isin(in_ids)].copy()
    pm['n_in'] = pm[in_ids].sum(axis=1)
    pm['n_out'] = pm[out_ids].sum(axis=1)
    pass_in = pm['n_in'] / len(in_ids) >= args.ingroup
    cand = pm[pass_in & (pm['n_out'] == 0)]
    nm = pm[pass_in & (pm['n_out'] >= 1) & (pm['n_out'] <= args.max_out)].copy()

    # Cross-check: our ingroup/outgroup rule must reproduce candidates.txt exactly.
    cand_file = res / 'candidates.txt'
    check_msg = 'candidates.txt not found; candidate rule not cross-checked'
    if cand_file.exists():
        ref = {l.strip() for l in open(cand_file) if l.strip()}
        ours = set(cand['source_proteome'] + '::' + cand['protein_id'])
        check_msg = (f'candidate cross-check: rule gives {len(ours)}, candidates.txt has '
                     f'{len(ref)}, symmetric difference {len(ours ^ ref)}')
        log(check_msg)
    log(f'matrix IN rows {len(pm)}; pass ingroup rule {int(pass_in.sum())}; candidates '
        f'{len(cand)}; near misses (1..{args.max_out} OUT) {len(nm)}')

    keep_ids = set(nm['protein_id'])
    nm_key = set(zip(nm['source_proteome'], nm['protein_id']))
    ev = pd.read_csv(res / 'presence_matrix.evalues.tsv', sep='\t', dtype=str)
    ev = ev[ev['protein_id'].isin(keep_ids)].set_index(['source_proteome', 'protein_id'])
    tg = pd.read_csv(res / 'presence_matrix.targets.tsv', sep='\t', dtype=str)
    tg = tg[tg['protein_id'].isin(keep_ids)].set_index(['source_proteome', 'protein_id'])

    # other_evidence: per OUT cell n_hits / max_qcov; also cells removed by filter 2.
    oe = {}
    oe_path = res / 'presence_matrix.other_evidence.tsv.gz'
    if oe_path.exists():
        with open_text(oe_path) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['protein_id'] in keep_ids and (r['source_proteome'], r['protein_id']) in nm_key:
                    oe[(r['source_proteome'], r['protein_id'], r['other_proteome'])] = r
    else:
        log('other_evidence not found; qcov from other_evidence not reported')

    # Raw diamond rows for the best kept target per OUT cell.
    hit_rows = []
    for (src, pid), row in nm.set_index(['source_proteome', 'protein_id']).iterrows():
        for o in out_ids:
            if row[o] == 1:
                hit_rows.append((src, pid, o, tg.at[(src, pid), o], float(ev.at[(src, pid), o])))
    raw = {}
    if not args.no_raw_hits:
        pairs_by_file = defaultdict(set)
        for src, pid, o, t, _ in hit_rows:
            pairs_by_file[(src, o)].add((pid, t))
        for (src, o), pairs in sorted(pairs_by_file.items()):
            got = read_raw_hits(res, src, o, pairs)
            if got is None:
                log(f'WARNING: no search_cache file for {src}_vs_{o}; pident missing')
                continue
            for (pid, t), rec in got.items():
                raw[(src, pid, o)] = rec
        n_match = 0
        for src, pid, o, t, e in hit_rows:
            rec = raw.get((src, pid, o))
            if rec is not None and math.isclose(float(rec['evalue']), e, rel_tol=1e-6, abs_tol=1e-300):
                n_match += 1
        log(f'raw diamond rows found for {len(raw)}/{len(hit_rows)} kept OUT cells; '
            f'E-value equal to presence_matrix.evalues.tsv for {n_match}')

    # Annotations.
    annot = None
    up = res / 'presence_matrix.uniprot.tsv'
    if up.exists():
        hdr = pd.read_csv(up, sep='\t', nrows=0).columns
        use = ['protein_id', 'source_proteome'] + [c for c in ANNOT_COLS if c in hdr]
        annot = pd.read_csv(up, sep='\t', dtype=str, usecols=use).fillna('')
        annot = annot.drop_duplicates(['source_proteome', 'protein_id'])

    # Families.
    source_of = dict(zip(nm['protein_id'], nm['source_proteome']))
    if args.family_method == 'mmseqs' and len(nm):
        pdir = args.proteins_dir or str(Path(args.config).resolve().parent / 'data_dir' / 'pep')
        fam = mmseqs_families(sorted(keep_ids), source_of, protein_file, pdir,
                              args.mmseqs, args.threads)
    else:
        fam = {p: p for p in keep_ids}

    # Protein-level table.
    rows = []
    for _, r in nm.iterrows():
        src, pid = r['source_proteome'], r['protein_id']
        hit_out = [o for o in out_ids if r[o] == 1]
        hit_out.sort(key=lambda o: float(ev.at[(src, pid), o]))
        lins = sorted({lineage.get(o, '') for o in hit_out})
        evs, pids, qcovs, alns, tgts, maxq, nh = [], [], [], [], [], [], []
        for o in hit_out:
            evs.append(f'{float(ev.at[(src, pid), o]):.3g}')
            tgts.append(tg.at[(src, pid), o])
            rec = raw.get((src, pid, o))
            pids.append(rec['pident'] if rec else 'NA')
            qcovs.append(rec['qcovhsp'] if rec else 'NA')
            alns.append(rec['length'] if rec else 'NA')
            e = oe.get((src, pid, o))
            maxq.append(e['max_qcov'] if e and e['max_qcov'] != '' else 'NA')
            nh.append(e['n_hits'] if e else 'NA')
        filtered = sorted(o for o in out_ids
                          if (src, pid, o) in oe and oe[(src, pid, o)]['status'] != 'kept')
        rows.append({
            'protein_id': pid, 'source_proteome': src, 'family_id': fam[pid],
            'n_in_species_with_hit': int(r['n_in']), 'n_in_total': len(in_ids),
            'n_out_species_hit': int(r['n_out']),
            'out_species': ','.join(hit_out),
            'out_lineages': ','.join(lins),
            'n_out_lineages': len(lins),
            'same_out_lineage': len(lins) == 1,
            'out_best_evalue': ','.join(evs),
            'out_pident': ','.join(pids),
            'out_qcov_best_hit': ','.join(qcovs),
            'out_aln_length': ','.join(alns),
            'out_max_qcov_any_hit': ','.join(maxq),
            'out_n_hits': ','.join(nh),
            'out_best_target': ','.join(tgts),
            'min_out_evalue': min(float(x) for x in evs),
            'out_paralog_filtered_species': ','.join(filtered),
        })
    prot = pd.DataFrame(rows)
    if annot is not None and len(prot):
        prot = prot.merge(annot, on=['protein_id', 'source_proteome'], how='left')
    fsize = prot['family_id'].map(prot['family_id'].value_counts()) if len(prot) else []
    if len(prot):
        prot.insert(3, 'family_size', fsize)
        prot = prot.sort_values(['n_out_species_hit', 'family_size', 'family_id', 'protein_id'],
                                ascending=[True, False, True, True])
    out_tsv = f'{args.out_prefix}.tsv.gz'
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    prot.to_csv(out_tsv, sep='\t', index=False, compression='gzip')

    # Family-level table.
    frows = []
    for fid, g in prot.groupby('family_id', sort=False):
        union = sorted({o for s in g['out_species'] for o in s.split(',')})
        lins = sorted({lineage.get(o, '') for o in union})
        best_ev = {}
        for s, e in zip(g['out_species'], g['out_best_evalue']):
            for o, x in zip(s.split(','), e.split(',')):
                best_ev[o] = min(best_ev.get(o, math.inf), float(x))
        sets = Counter(g['out_species'])
        pf = Counter()
        desc = Counter()
        if 'uniprot_pfam_ids' in g.columns:
            for s in g['uniprot_pfam_ids'].fillna(''):
                pf.update(set(split_ids(s)))
            desc.update(d for d in g['uniprot_description'].fillna('') if d)
        frows.append({
            'family_id': fid, 'n_proteins': len(g),
            'n_in_species_in_family': g['source_proteome'].nunique(),
            'in_species_in_family': ','.join(sorted(g['source_proteome'].unique())),
            'n_out_species_union': len(union),
            'out_species_union': ','.join(union),
            'out_lineages_union': ','.join(lins),
            'same_out_lineage': len(lins) == 1,
            'min_n_out_per_protein': int(g['n_out_species_hit'].min()),
            'max_n_out_per_protein': int(g['n_out_species_hit'].max()),
            'n_distinct_out_sets': len(sets),
            'most_common_out_set': sets.most_common(1)[0][0],
            'best_evalue_per_out': ','.join(f'{o}:{best_ev[o]:.3g}' for o in union),
            'min_out_evalue': min(best_ev.values()),
            'n_annotated': int((g.get('uniprot_accession', pd.Series(dtype=str)).fillna('') != '').sum()),
            'pfam_ids': ','.join(f'{k}({v})' for k, v in pf.most_common()),
            'top_description': desc.most_common(1)[0][0] if desc else '',
            'members': ','.join(g['protein_id']),
        })
    fams = pd.DataFrame(frows)
    if len(fams):
        fams = fams.sort_values(['n_out_species_union', 'n_proteins'], ascending=[True, False])
    fams.to_csv(f'{args.out_prefix}_families.tsv', sep='\t', index=False)

    # Per OUT species table.
    orows = []
    for o in out_ids:
        p1 = prot[prot['out_species'] == o]
        pall = prot[prot['out_species'].str.split(',').apply(lambda L: o in L)]
        f1 = fams[fams['out_species_union'] == o] if len(fams) else fams
        fall = fams[fams['out_species_union'].str.split(',').apply(lambda L: o in L)] if len(fams) else fams
        orows.append({'out_species': o, 'lineage': lineage.get(o, ''),
                      'proteins_sole_out_hit': len(p1), 'proteins_any': len(pall),
                      'families_sole_out_hit': len(f1), 'families_any': len(fall)})
    outsp = pd.DataFrame(orows).sort_values(['families_sole_out_hit', 'proteins_sole_out_hit'],
                                            ascending=False)
    outsp.to_csv(f'{args.out_prefix}_out_species.tsv', sep='\t', index=False)

    # Pfam counts: near misses vs candidates, among proteins with UniProt annotation.
    pfam_lines = []
    if annot is not None:
        a = annot.set_index(['source_proteome', 'protein_id'])
        def pfam_counts(df):
            keys = list(zip(df['source_proteome'], df['protein_id']))
            sub = a.reindex(keys)
            annotated = sub['uniprot_accession'].fillna('') != ''
            cnt = Counter()
            for s in sub.loc[annotated, 'uniprot_pfam_ids'].fillna(''):
                cnt.update(set(split_ids(s)))
            return cnt, int(annotated.sum())
        nm_cnt, nm_ann = pfam_counts(nm)
        cd_cnt, cd_ann = pfam_counts(cand)
        names = {}
        if 'uniprot_pfam_names' in annot.columns:
            for i, n in zip(annot['uniprot_pfam_ids'], annot['uniprot_pfam_names']):
                for k, v in zip(split_ids(i), split_ids(n)):
                    names.setdefault(k, v)
        prow = []
        for k in set(nm_cnt) | set(cd_cnt):
            fn = nm_cnt[k] / nm_ann if nm_ann else 0
            fc = cd_cnt[k] / cd_ann if cd_ann else 0
            prow.append({'pfam_id': k, 'pfam_name': names.get(k, ''),
                         'near_miss_proteins': nm_cnt[k], 'candidate_proteins': cd_cnt[k],
                         'near_miss_frac': fn, 'candidate_frac': fc,
                         'log2_ratio_pc0.5': math.log2(((nm_cnt[k] + 0.5) / max(nm_ann, 1)) /
                                                      ((cd_cnt[k] + 0.5) / max(cd_ann, 1)))})
        pf = pd.DataFrame(prow).sort_values(['near_miss_proteins'], ascending=False)
        pf.to_csv(f'{args.out_prefix}_pfam_vs_candidates.tsv', sep='\t', index=False)
        pfam_lines.append(f'Pfam comparison uses UniProt-annotated proteins only: near misses '
                          f'{nm_ann}/{len(nm)}, candidates {cd_ann}/{len(cand)}.')
        top = pf[pf['near_miss_proteins'] >= args.min_pfam_nm].sort_values(
            'log2_ratio_pc0.5', ascending=False).head(args.top)
        pfam_lines.append(f'Top {args.top} Pfam by log2((nm+0.5)/nm_annotated / '
                          f'((cand+0.5)/cand_annotated)), Pfam with >= {args.min_pfam_nm} '
                          'near-miss proteins. Descriptive ratio only; no test.')
        pfam_lines.append(top[['pfam_id', 'pfam_name', 'near_miss_proteins', 'candidate_proteins',
                               'log2_ratio_pc0.5']].to_string(index=False))
        pfam_lines.append(f'Top {args.top} Pfam by near-miss protein count:')
        pfam_lines.append(pf.head(args.top)[['pfam_id', 'pfam_name', 'near_miss_proteins',
                                             'candidate_proteins']].to_string(index=False))

    # Summary.
    S = []
    S.append(f'results: {res}')
    S.append(f'config: {args.config}')
    S.append(f'rule: IN present >= {args.ingroup} of {len(in_ids)} (own proteome counts); '
             f'near miss = 1..{args.max_out} OUT proteomes present')
    S.append(check_msg)
    S.append(f'IN matrix rows: {len(pm)}; pass ingroup rule: {int(pass_in.sum())}; '
             f'candidates (0 OUT): {len(cand)}')
    for k in range(1, args.max_out + 1):
        pk = prot[prot['n_out_species_hit'] == k]
        S.append(f'near-miss proteins with {k} OUT species: {len(pk)}; same lineage: '
                 f'{int(pk["same_out_lineage"].sum())}')
    S.append(f'near-miss proteins total: {len(prot)}; families (mmseqs on near misses): '
             f'{len(fams)}')
    if len(fams):
        for k in sorted(fams['n_out_species_union'].unique()):
            fk = fams[fams['n_out_species_union'] == k]
            S.append(f'families with {k} OUT species in union: {len(fk)}; same lineage: '
                     f'{int(fk["same_out_lineage"].sum())}')
        multi = fams[fams['n_out_species_union'] >= 2]
        S.append(f'single-lineage families with >= 2 OUT species: '
                 f'{int(multi["same_out_lineage"].sum())} of {len(multi)}')
        S.append(f'families with >= 2 proteins: {int((fams["n_proteins"] >= 2).sum())}; '
                 f'with >= 2 IN species: {int((fams["n_in_species_in_family"] >= 2).sum())}')
    S.append('per OUT species (sorted by families where it is the only OUT hit):')
    S.append(outsp.to_string(index=False))
    S.extend(pfam_lines)
    S.append(f'runtime_seconds: {time.time() - t0:.1f}')
    txt = '\n'.join(S) + '\n'
    with open(f'{args.out_prefix}_summary.txt', 'w') as fh:
        fh.write(txt)
    print(txt)


if __name__ == '__main__':
    main()
