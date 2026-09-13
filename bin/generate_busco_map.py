#!/usr/bin/env python3
"""Build the 2-column busco_id<TAB>protein_id map nf_NovInvenio/bin/score_controls.py's
--busco-map expects, from one species' BUSCO full_table.tsv, restricted to the busco
anchor ids an actual controls CSV references -- fixes the gap where score_controls.py
was run without --busco-map at all, silently reporting every busco-anchored negative
control as unresolved instead of scoring it (see notes/superpowers/specs/2026-09-13-
cluster-vs-pairwise-sensitivity-design.md, Analysis 1).
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
from busco import parse_busco_full_table  # noqa: E402


def filter_wanted(rows, wanted_ids):
    """First Complete row per busco_id, restricted to wanted_ids. rows: iterable of
    (busco_id, species, protein_id, length)."""
    mapping = {}
    for busco_id, _species, protein_id, _length in rows:
        if busco_id in wanted_ids and busco_id not in mapping:
            mapping[busco_id] = protein_id
    return mapping


def load_wanted_busco_ids(controls_csv):
    wanted = set()
    with open(controls_csv) as fh:
        for row in csv.DictReader(fh):
            if (row.get('anchor_type') or '').strip() == 'busco':
                wanted.add((row.get('anchor') or '').strip())
    return wanted


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full-table', required=True, dest='full_table',
                    help="one species' BUSCO full_table.tsv")
    ap.add_argument('--species', required=True, help='label passed through to parse_busco_full_table')
    ap.add_argument('--controls', required=True, help='controls CSV to restrict ids to')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    wanted = load_wanted_busco_ids(args.controls)
    rows = list(parse_busco_full_table(args.full_table, args.species))
    mapping = filter_wanted(rows, wanted)

    missing = wanted - mapping.keys()
    if missing:
        print(f'WARNING: {len(missing)} busco id(s) not Complete in {args.species}: '
              f'{sorted(missing)}', file=sys.stderr)

    with open(args.output, 'w') as fh:
        for busco_id in sorted(mapping):
            fh.write(f'{busco_id}\t{mapping[busco_id]}\n')


if __name__ == '__main__':
    main()
