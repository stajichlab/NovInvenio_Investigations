#!/usr/bin/env python3
"""Assemble Tasks 5-8's TSV outputs into one markdown report with the Analysis 4
decision-framework recommendation -- see notes/superpowers/specs/2026-09-13-
cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def render_table(rows):
    if not rows:
        return '_(no rows)_\n'
    cols = list(rows[0].keys())
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for row in rows:
        lines.append('| ' + ' | '.join(str(row.get(c, '')) for c in cols) + ' |')
    return '\n'.join(lines) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tier-comparison-dir', required=True, dest='dir')
    ap.add_argument('--label', default='pezizo_set1')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    d = Path(args.dir)
    sections = []
    sections.append(f'# Cluster-vs-pairwise sensitivity: {args.label}\n')

    sections.append('## Analysis 1: controls recall/FP by tier\n')
    sections.append(render_table(read_tsv(d / f'{args.label}.tier_summary.tsv')))
    sections.append('\n### Per-control outcomes\n')
    sections.append(render_table(read_tsv(d / f'{args.label}.per_control.tsv')))

    gain_conc = d / f'{args.label}.gain.concordance.tsv'
    if gain_conc.exists():
        sections.append('\n## Analysis 2: genome-wide concordance (gains)\n')
        sections.append(render_table(read_tsv(gain_conc)))
    loss_conc = d / f'{args.label}.loss.concordance.tsv'
    if loss_conc.exists():
        sections.append('\n## Analysis 2: genome-wide concordance (losses)\n')
        sections.append(render_table(read_tsv(loss_conc)))

    cost = d / f'{args.label}.cost.tsv'
    if cost.exists():
        sections.append('\n## Analysis 3: real compute cost (CPU-hours)\n')
        sections.append(render_table(read_tsv(cost)))

    sections.append('\n## Analysis 4: decision framework\n')
    sections.append(
        '_Fill in by hand after reading the tables above — see the spec\'s Analysis 4 '
        'for the exact framework (escalate Tier R near-misses to Tier P; treat '
        'candidate-count inflation and HEX1-style clustering-inherent misses as the '
        'two conditions that would invalidate the framework)._\n'
    )

    Path(args.output).write_text('\n'.join(sections))


if __name__ == '__main__':
    main()
