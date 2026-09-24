#!/usr/bin/env python3
"""Analysis 2 follow-up: explain every disagreement between Tier C+H and Tier P.

For each Tier P candidate that Tier C+H misses, and each Tier C+H candidate that
Tier P rejects, assign one cause class. The classes follow the failure modes in
notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md
(anchor-tracing: family membership vs. presence_matrix disagreement).

Misses (in P, not in C+H), checked in this order:

  clustering-level
    unprofiled_family     the protein's mmseqs cluster is not in families.tsv
                          (below the family-size floor, e.g. a singleton), so no HMM
                          exists for it (failure mode 3, LAH-like). A protein that is
                          in no cluster at all is also counted here.
    oversized_family      the cluster is in oversized_families.tsv (not profiled).
    cluster_fragmentation cluster membership covers < --min-frac of the query group,
                          and the HMM does not recover the missing species either.
  HMM-level
    hmm_outgroup_paralog  the HMM calls the family present in the other group, AND
                          the protein's registered within-genome paralog (from
                          --self-hits) is a co-member of its family (failure mode 1,
                          HEX1-like: a merged paralog carries the HMM hit).
    hmm_outgroup_other    the HMM calls the family present in the other group, with
                          no co-member paralog (failure mode 2, ADA1/HAM5-like:
                          the profile turns weak signal into "present").
    hmm_query_undercall   membership covers >= --min-frac of the query group, but the
                          HMM presence call does not (coverage floor / E-value cut).
  other
    unexplained           C+H's matrix row passes both fractions, but the protein is
                          not in C+H's candidates.txt.

Extras (in C+H, not in P):

    p_other_presence      P's search finds the protein in the other group (P found a
                          homolog the family HMM missed).
    p_query_short         P finds it in < --min-frac of the query group.
    p_no_hits             the protein has no P matrix row. build_presence_matrix.py
                          builds rows from surviving hits only, so this means P's
                          search found no significant hit in any other proteome,
                          while the family HMM found homologs in the query group.
    unexplained           P's row passes both fractions, but it is not in P's
                          candidates.txt.

--query-group IN = gains (seed group = ingroup); OUT = losses (seed = outgroup).

Outputs: --output-summary (counts per direction/set/class; small, class 2) and an
optional --output-detail (one row per protein; candidate-bearing, class 3 -- keep it
gitignored).
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

NII_PIPELINE_LIB = Path('/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
sys.path.insert(0, str(NII_PIPELINE_LIB))
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config  # noqa: E402

CLUSTERING = {'unprofiled_family', 'oversized_family', 'cluster_fragmentation'}
HMM = {'hmm_outgroup_paralog', 'hmm_outgroup_other', 'hmm_query_undercall'}


def load_candidates(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_matrix(path):
    """protein_id -> (source_proteome, {proteome: 0/1})."""
    out = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        cols = [c for c in r.fieldnames if c not in ('protein_id', 'source_proteome')]
        for row in r:
            out[row['protein_id']] = (row['source_proteome'],
                                      {c: int(row[c]) for c in cols})
    return out


def load_clusters(path):
    """(member -> rep, rep -> [members])."""
    member_to_rep, rep_members = {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            rep, member = line.rstrip('\n').split('\t')[:2]
            member_to_rep[member] = rep
            rep_members[rep].append(member)
    return member_to_rep, dict(rep_members)


def load_first_column(path, skip_header=True):
    if not path:
        return set()
    with open(path) as fh:
        if skip_header:
            fh.readline()
        return {line.split('\t')[0].strip() for line in fh if line.strip()}


def load_profiled_reps(path):
    reps = set()
    with open(path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1]:
                reps.add(parts[1])
    return reps


def load_paralog_map(paths):
    paralog_of = {}
    for p in paths or []:
        with open(p) as fh:
            fh.readline()
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    paralog_of[parts[0]] = parts[1]
    return paralog_of


def frac(presence, group):
    return sum(presence.get(s, 0) for s in group) / len(group) if group else 0.0


def classify_miss(pid, ch_matrix, member_to_rep, rep_members, profiled, oversized,
                  protein_to_proteome, paralog_of, query, other, min_frac, other_max):
    rep = member_to_rep.get(pid)
    if rep is None:
        return 'unprofiled_family'
    if rep in oversized:
        return 'oversized_family'
    if rep not in profiled:
        return 'unprofiled_family'
    members = rep_members[rep]
    mem_species = {protein_to_proteome.get(m) for m in members}
    mem_q = len(mem_species & query) / len(query) if query else 0.0
    row = ch_matrix.get(pid) or ch_matrix.get(rep)
    if row is None:
        # Profiled family with no matrix row: the pipeline dropped it downstream.
        return 'unexplained'
    presence = row[1]
    hmm_q, hmm_o = frac(presence, query), frac(presence, other)
    if hmm_o > other_max:
        p = paralog_of.get(pid)
        return 'hmm_outgroup_paralog' if p and p != pid and p in set(members) else 'hmm_outgroup_other'
    if hmm_q < min_frac:
        return 'cluster_fragmentation' if mem_q < min_frac else 'hmm_query_undercall'
    return 'unexplained'


def classify_extra(pid, p_matrix, query, other, min_frac, other_max):
    row = p_matrix.get(pid)
    if row is None:
        return 'p_no_hits'
    presence = row[1]
    if frac(presence, other) > other_max:
        return 'p_other_presence'
    if frac(presence, query) < min_frac:
        return 'p_query_short'
    return 'unexplained'


def level_of(cls):
    if cls in CLUSTERING:
        return 'clustering'
    if cls in HMM:
        return 'hmm'
    if cls.startswith('p_'):
        return 'pairwise'
    return 'other'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--label', required=True, help='clade label written to every row')
    ap.add_argument('--direction', required=True, choices=['gain', 'loss'])
    ap.add_argument('--query-group', required=True, choices=['IN', 'OUT'], dest='query_group')
    ap.add_argument('--config', required=True)
    ap.add_argument('--p-candidates', required=True, dest='p_candidates')
    ap.add_argument('--p-matrix', required=True, dest='p_matrix')
    ap.add_argument('--ch-candidates', required=True, dest='ch_candidates')
    ap.add_argument('--ch-matrix', required=True, dest='ch_matrix')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv')
    ap.add_argument('--families', required=True)
    ap.add_argument('--oversized-families', default=None, dest='oversized_families')
    ap.add_argument('--self-hits', nargs='*', default=[], dest='self_hits')
    ap.add_argument('--min-frac', type=float, default=0.75, dest='min_frac')
    ap.add_argument('--other-max-frac', type=float, default=0.0, dest='other_max_frac')
    ap.add_argument('--output-summary', required=True, dest='output_summary')
    ap.add_argument('--output-detail', default=None, dest='output_detail')
    args = ap.parse_args()

    samples = parse_config(args.config)
    ingroup = {s.short for s in samples if s.group in INGROUP_ROLES}
    outgroup = {s.short for s in samples if s.group in OUTGROUP_ROLES}
    query, other = (ingroup, outgroup) if args.query_group == 'IN' else (outgroup, ingroup)

    p_cand = load_candidates(args.p_candidates)
    ch_cand = load_candidates(args.ch_candidates)
    p_matrix = load_matrix(args.p_matrix)
    ch_matrix = load_matrix(args.ch_matrix)
    member_to_rep, rep_members = load_clusters(args.cluster_tsv)
    profiled = load_profiled_reps(args.families)
    oversized = load_first_column(args.oversized_families)
    paralog_of = load_paralog_map(args.self_hits)
    protein_to_proteome = {pid: v[0] for pid, v in p_matrix.items()}
    protein_to_proteome.update({pid: v[0] for pid, v in ch_matrix.items()})
    # candidates.txt lines are SHORT::protein_id; fill any proteome still unknown.
    for line in p_cand | ch_cand:
        short, _, pid = line.partition('::')
        protein_to_proteome.setdefault(pid, short)

    detail = []
    for line in sorted(p_cand - ch_cand):
        pid = line.partition('::')[2]
        cls = classify_miss(pid, ch_matrix, member_to_rep, rep_members, profiled, oversized,
                            protein_to_proteome, paralog_of, query, other,
                            args.min_frac, args.other_max_frac)
        detail.append(('miss', line, cls))
    for line in sorted(ch_cand - p_cand):
        pid = line.partition('::')[2]
        cls = classify_extra(pid, p_matrix, query, other, args.min_frac, args.other_max_frac)
        detail.append(('extra', line, cls))

    counts = Counter((kind, cls) for kind, _, cls in detail)
    totals = Counter(kind for kind, _, _ in detail)
    with open(args.output_summary, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['label', 'direction', 'set', 'level', 'class', 'n', 'frac_of_set',
                    'n_p_candidates', 'n_ch_candidates', 'n_shared'])
        for (kind, cls), n in sorted(counts.items(), key=lambda kv: (kv[0][0], -kv[1])):
            w.writerow([args.label, args.direction, kind, level_of(cls), cls, n,
                        round(n / totals[kind], 4), len(p_cand), len(ch_cand),
                        len(p_cand & ch_cand)])
    if args.output_detail:
        with open(args.output_detail, 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['set', 'candidate', 'class'])
            w.writerows(detail)
    print(f'{args.label} {args.direction}: {totals["miss"]} misses, {totals["extra"]} extras',
          file=sys.stderr)


if __name__ == '__main__':
    main()
