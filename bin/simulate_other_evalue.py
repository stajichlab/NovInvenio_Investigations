#!/usr/bin/env python3
"""Simulate a stricter E-value cutoff for OTHER-group presence on an existing Tier P run,
without re-running the pipeline (notes/diamond-sensitivity/README.md, follow-up 2).

Input is a finished run's presence_matrix.tsv + presence_matrix.evalues.tsv (or the
loss_ pair). The evalues file holds the E-value of every hit that counted as presence
(after the run's own --evalue and paralog filters). For each threshold T, an
other-group cell stays present only if its E-value < T; query-group cells are left as
the run called them. A protein is a candidate when it is present in >= --min-frac of
the query group and in <= --other-max-frac of the other group -- the same keep rule
as build_presence_matrix.py. At T equal to the run's own --evalue the simulation must
reproduce the run's candidates file exactly; main() checks this.

Per threshold it reports the candidate count, overlap with a reference candidate list,
and how many candidates have a TBLASTN hit in any other-group genome (from one or
more --tblastn summaries; a candidate in none of them is counted as unknown). A
TBLASTN hit is independent evidence against a candidate being truly absent from the
other group. With --write-matrix-dir it also writes the thresholded presence matrix
per T, for score_controls.py.
"""
import argparse
import csv
import sys
from pathlib import Path

NII_PIPELINE_LIB = Path('/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
sys.path.insert(0, str(NII_PIPELINE_LIB))
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config  # noqa: E402


def other_present(s, pres, evs, threshold):
    """An other-group cell counts when the run called it present and, if it has a
    recorded E-value, that E-value is below the threshold."""
    if str(pres.get(s, '0')) != '1':
        return False
    e = evs.get(s, '')
    return True if e in ('', None) else float(e) < threshold


def is_candidate(pres, evs, query, other, min_frac, other_max, threshold):
    q = sum(str(pres.get(s, '0')) == '1' for s in query) / len(query)
    o = (sum(other_present(s, pres, evs, threshold) for s in other) / len(other)) if other else 0.0
    return q >= min_frac and o <= other_max


def load_tblastn(paths):
    """protein_id -> True/False (any other-group genome hit), merged over files."""
    out = {}
    for p in paths:
        with open(p) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            cols = [c for c in r.fieldnames if c != 'protein_id']
            for row in r:
                out[row['protein_id']] = out.get(row['protein_id'], False) or any(
                    row[c] == '1' for c in cols)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--label', required=True)
    ap.add_argument('--direction', required=True, choices=['gain', 'loss'])
    ap.add_argument('--run-dir', required=True, type=Path, dest='run_dir')
    ap.add_argument('--config', required=True)
    ap.add_argument('--run-evalue', type=float, default=1e-5, dest='run_evalue')
    ap.add_argument('--thresholds', nargs='+', type=float, required=True)
    ap.add_argument('--min-frac', type=float, default=0.75, dest='min_frac')
    ap.add_argument('--other-max-frac', type=float, default=0.0, dest='other_max')
    ap.add_argument('--reference', nargs=2, action='append', default=[],
                    metavar=('NAME', 'CANDIDATES'), help='Repeatable reference candidate lists')
    ap.add_argument('--tblastn', nargs='+', default=[], help='tblastn summary TSVs')
    ap.add_argument('--write-matrix-dir', type=Path, default=None, dest='matrix_dir')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()

    pre = '' if a.direction == 'gain' else 'loss_'
    samples = parse_config(a.config)
    ingroup = {s.short for s in samples if s.group in INGROUP_ROLES}
    outgroup = {s.short for s in samples if s.group in OUTGROUP_ROLES}
    query, other = (ingroup, outgroup) if a.direction == 'gain' else (outgroup, ingroup)

    with open(a.run_dir / f'{pre}presence_matrix.tsv') as fh:
        mr = csv.DictReader(fh, delimiter='\t')
        header = mr.fieldnames
        matrix = list(mr)
    with open(a.run_dir / f'{pre}presence_matrix.evalues.tsv') as fh:
        evalues = {(r['source_proteome'], r['protein_id']): r for r in csv.DictReader(fh, delimiter='\t')}
    tb = load_tblastn(a.tblastn)
    refs = {}
    for name, path in a.reference:
        with open(path) as fh:
            refs[name] = {ln.strip() for ln in fh if ln.strip()}

    run_cands = {ln.strip() for ln in open(a.run_dir / f'{pre}candidates.txt') if ln.strip()}
    rows = []
    for t in sorted(set(a.thresholds + [a.run_evalue]), reverse=True):
        cands = set()
        kept_rows = []
        for r in matrix:
            evs = evalues.get((r['source_proteome'], r['protein_id']), {})
            if is_candidate(r, evs, query, other, a.min_frac, a.other_max, t):
                cands.add(f"{r['source_proteome']}::{r['protein_id']}")
            if a.matrix_dir:
                kept_rows.append({c: (('1' if other_present(c, r, evs, t) else '0')
                                      if c in other else r[c]) for c in header})
        if t == a.run_evalue and cands != run_cands:
            sys.exit(f'ERROR: simulation at the run evalue {t} does not reproduce '
                     f'{pre}candidates.txt ({len(cands)} vs {len(run_cands)}); aborting')
        pids = [c.split('::', 1)[1] for c in cands]
        tb_yes = sum(1 for p in pids if tb.get(p) is True)
        tb_no = sum(1 for p in pids if tb.get(p) is False)
        row = {'label': a.label, 'direction': a.direction, 'other_evalue': f'{t:g}',
               'n_candidates': len(cands), 'tblastn_other_hit': tb_yes,
               'tblastn_no_hit': tb_no, 'tblastn_unknown': len(cands) - tb_yes - tb_no,
               'tblastn_hit_frac_of_known': round(tb_yes / (tb_yes + tb_no), 4) if tb_yes + tb_no else ''}
        for name, ref in refs.items():
            row[f'shared_{name}'] = len(cands & ref)
        rows.append(row)
        if a.matrix_dir:
            a.matrix_dir.mkdir(parents=True, exist_ok=True)
            with open(a.matrix_dir / f'{a.label}.{a.direction}.E{t:g}.presence_matrix.tsv', 'w',
                      newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=header, delimiter='\t', lineterminator='\n')
                w.writeheader()
                w.writerows(kept_rows)
    with open(a.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{a.label} {a.direction}: reproduced run candidates at E={a.run_evalue:g} '
          f'({len(run_cands)})', file=sys.stderr)


if __name__ == '__main__':
    main()
