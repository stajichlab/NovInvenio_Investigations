#!/usr/bin/env python3
"""Assemble Tasks 5-8's TSV outputs into one markdown report with the Analysis 4
decision-framework recommendation -- see notes/superpowers/specs/2026-09-13-
cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv
import sys
from pathlib import Path


PLACEHOLDER = '_Fill in by hand after reading the tables above'


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def check_overwrite_guard(output_path, force):
    """Refuse to clobber a hand-written report.

    Returns None if it's safe to proceed (no existing file, existing file still
    has the unfilled Analysis 4 placeholder, or --force was passed). Returns an
    error message string if the write should be refused.
    """
    out = Path(output_path)
    if force or not out.exists():
        return None
    existing = out.read_text()
    if PLACEHOLDER in existing:
        return None
    return (
        f'refusing to overwrite {output_path}: it already exists and does not '
        f'contain the unfilled Analysis 4 placeholder ("{PLACEHOLDER}..."), so it '
        'looks hand-edited (e.g. a filled-in Analysis 4 conclusion or a hand-merged '
        'subsection). Pass --force to regenerate it anyway and discard those '
        'hand edits.'
    )


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
    ap.add_argument('--force', action='store_true',
                     help='Overwrite --output even if it already exists and does not '
                          'contain the unfilled Analysis 4 placeholder (i.e. it looks '
                          'hand-edited). Without this, such an existing file is left '
                          'untouched and the script exits with an error.')
    args = ap.parse_args()

    guard_error = check_overwrite_guard(args.output, args.force)
    if guard_error is not None:
        print(guard_error, file=sys.stderr)
        raise SystemExit(1)

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
        sections.append('\n## Analysis 3: real compute cost (wall-hours)\n')
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
