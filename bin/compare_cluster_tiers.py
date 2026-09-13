#!/usr/bin/env python3
"""Analysis 1: run pairwise (Tier P), cluster+HMM (Tier C+H), cluster-membership-only
(Tier C), and refined-cluster (Tier R) scoring through nf_NovInvenio/bin/score_controls.py
against one controls CSV, and assemble one per-control + one per-tier-summary table --
see notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

SCORE_CONTROLS = '/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/score_controls.py'


def parse_summary_tsv(path):
    with open(path) as fh:
        r = csv.reader(fh, delimiter='\t')
        next(r)  # header
        return {row[0]: row[1] for row in r if row}


def parse_per_control_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def run_score_controls(*, matrix, controls, config, output, cluster_tsv=None,
                       families=None, presence_mode='hmm', busco_map=None, cpus=1):
    cmd = [sys.executable, SCORE_CONTROLS, '--controls', controls, '--matrix', matrix,
           '--config', config, '--output', str(output), '--presence-mode', presence_mode,
           '--cpus', str(cpus)]
    if cluster_tsv:
        cmd += ['--cluster-tsv', cluster_tsv, '--families', families]
    if busco_map:
        cmd += ['--busco-map', busco_map]
    subprocess.run(cmd, check=True)


def build_tier_comparison(tier_summaries, tier_per_control):
    per_control_rows = []
    for tier, rows in tier_per_control.items():
        for row in rows:
            per_control_rows.append({**row, 'tier': tier})
    summary_rows = [{'tier': tier, **summary} for tier, summary in tier_summaries.items()]
    return per_control_rows, summary_rows


def write_tsv(rows, path):
    if not rows:
        Path(path).write_text('')
        return
    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--controls', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--pairwise-matrix', required=True, dest='pairwise_matrix')
    ap.add_argument('--cluster-hmm-matrix', required=True, dest='cluster_hmm_matrix')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv')
    ap.add_argument('--families', required=True)
    ap.add_argument('--refined-cluster-tsv', required=True, dest='refined_cluster_tsv')
    ap.add_argument('--refined-families', required=True, dest='refined_families')
    ap.add_argument('--busco-map', default=None, dest='busco_map')
    ap.add_argument('--out-dir', required=True, dest='out_dir')
    ap.add_argument('--label', required=True, help='output filename prefix, e.g. pezizo_set1')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tiers = {
        'P': dict(matrix=args.pairwise_matrix, cluster_tsv=None, families=None, presence_mode='hmm'),
        'C+H': dict(matrix=args.cluster_hmm_matrix, cluster_tsv=args.cluster_tsv,
                    families=args.families, presence_mode='hmm'),
        'C': dict(matrix=args.cluster_hmm_matrix, cluster_tsv=args.cluster_tsv,
                  families=args.families, presence_mode='cluster_membership'),
        'R': dict(matrix=args.cluster_hmm_matrix, cluster_tsv=args.refined_cluster_tsv,
                  families=args.refined_families, presence_mode='cluster_membership'),
    }

    tier_summaries, tier_per_control = {}, {}
    for tier, cfg in tiers.items():
        output = out_dir / f'{args.label}.{tier.replace("+", "")}.controls_scored.tsv'
        run_score_controls(matrix=cfg['matrix'], controls=args.controls, config=args.config,
                           output=output, cluster_tsv=cfg['cluster_tsv'],
                           families=cfg['families'], presence_mode=cfg['presence_mode'],
                           busco_map=args.busco_map)
        summary_path = output.with_suffix('').with_suffix('.summary.tsv')
        tier_summaries[tier] = parse_summary_tsv(summary_path)
        tier_per_control[tier] = parse_per_control_tsv(output)

    per_control_rows, summary_rows = build_tier_comparison(tier_summaries, tier_per_control)
    write_tsv(per_control_rows, out_dir / f'{args.label}.per_control.tsv')
    write_tsv(summary_rows, out_dir / f'{args.label}.tier_summary.tsv')
    print(f'Wrote {out_dir}/{args.label}.per_control.tsv and .tier_summary.tsv',
          file=sys.stderr)


if __name__ == '__main__':
    main()
