#!/usr/bin/env python3
"""Analysis 2: genome-wide precision/recall/Jaccard of Tier C/C+H/R's novelty (or
loss) candidate lists against Tier P's, plus Tier R's candidate-count-inflation
check (the over-splitting signature the 16/2-control set can't see) -- see
notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv


def load_candidates(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def precision_recall(predicted, gold):
    if not predicted:
        precision = 0.0 if gold else 1.0
    else:
        precision = len(predicted & gold) / len(predicted)
    recall = (len(predicted & gold) / len(gold)) if gold else 1.0
    return precision, recall


def count_inflation(baseline, refined):
    delta = len(refined) - len(baseline)
    delta_pct = (delta / len(baseline) * 100) if baseline else (0.0 if not refined else float('inf'))
    return {
        'baseline_count': len(baseline),
        'refined_count': len(refined),
        'delta': delta,
        'delta_pct': round(delta_pct, 1),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gold', required=True, help="Tier P's candidates.txt/loss_candidates.txt")
    ap.add_argument('--tier', action='append', nargs=2, metavar=('NAME', 'PATH'),
                    dest='tiers', required=True,
                    help='Repeatable: --tier C path/to/candidates.txt --tier "C+H" ... '
                         '--tier R path/to/refined_candidates.txt')
    ap.add_argument('--baseline-tier', default='C', dest='baseline_tier',
                    help='Which --tier name is the count_inflation baseline for Tier R '
                         '(default: C, i.e. is R inflated relative to unrefined clustering)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    gold = load_candidates(args.gold)
    tier_candidates = {name: load_candidates(path) for name, path in args.tiers}

    rows = []
    baseline = tier_candidates.get(args.baseline_tier)
    for name, candidates in tier_candidates.items():
        precision, recall = precision_recall(candidates, gold)
        row = {
            'tier': name,
            'n_candidates': len(candidates),
            'precision_vs_P': round(precision, 4),
            'recall_vs_P': round(recall, 4),
            'jaccard_vs_P': round(jaccard(candidates, gold), 4),
        }
        if baseline is not None and name != args.baseline_tier:
            row.update({f'{k}_vs_{args.baseline_tier}': v
                       for k, v in count_inflation(baseline, candidates).items()})
        rows.append(row)

    fieldnames = sorted({k for row in rows for k in row})
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
