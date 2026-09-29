#!/usr/bin/env python3.12
"""Turn a pangenome presence matrix into a SANS splits file.

Each family's carrier set is one bipartition of the strains. The split
weight is the number of families with that bipartition (a pattern and its
complement are the same split). Families carried by 0, 1, n-1 or n strains
are trivial splits and are skipped.

Output format is what `SANS -s` reads: weight, then the strain IDs on the
smaller side, tab-separated, sorted by weight descending.

Usage:
  gene_content_splits.py --matrix presence_matrix.rescued.tsv \
      --present present,genome_only --out gene_content.splits
"""
import argparse
import collections
import gzip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--present', default='present,genome_only',
                    help='cell values counted as present (comma-separated)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--strains-out', help='write the strain list (one per line)')
    args = ap.parse_args()

    present = set(args.present.split(','))
    opener = gzip.open if args.matrix.endswith('.gz') else open
    counts = collections.Counter()
    n_fam = n_trivial = 0
    with opener(args.matrix, 'rt') as fh:
        strains = fh.readline().rstrip('\n').split('\t')[1:]
        n = len(strains)
        full = frozenset(range(n))
        for line in fh:
            cells = line.rstrip('\n').split('\t')[1:]
            carriers = frozenset(i for i, v in enumerate(cells) if v in present)
            n_fam += 1
            k = len(carriers)
            if k <= 1 or k >= n - 1:
                n_trivial += 1
                continue
            other = full - carriers
            # canonical side: the smaller one; ties broken by sorted index tuple
            side = min((carriers, other), key=lambda s: (len(s), sorted(s)))
            counts[side] += 1

    with open(args.out, 'w') as out:
        for side, w in sorted(counts.items(), key=lambda kv: (-kv[1], sorted(kv[0]))):
            out.write(str(w) + '\t' + '\t'.join(strains[i] for i in sorted(side)) + '\n')
    if args.strains_out:
        with open(args.strains_out, 'w') as out:
            out.write('\n'.join(strains) + '\n')
    print(f'{n_fam} families, {n_trivial} trivial (<=1 or >=n-1 carriers), '
          f'{len(counts)} distinct non-trivial splits from {n_fam - n_trivial} families; '
          f'{n} strains')


if __name__ == '__main__':
    main()
