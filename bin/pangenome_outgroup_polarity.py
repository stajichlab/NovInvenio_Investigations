#!/usr/bin/env python3
"""Per frequency bin: outgroup presence and gain/loss calls, for one pangenome run.

Used for notes/pangenome-method-investigations/2026-09-24-rescued-rerun-and-outgroup-polarity.md.
"present" = present or genome_only (same as NovInvenio PresenceMatrix.is_present).

Usage:
  bin/pangenome_outgroup_polarity.py --matrix presence_matrix.rescued.tsv \
      --freq frequency_table.tsv --pairs cooccurring_pairs.tsv.zst --outgroup Uree

Inputs: presence matrix (family x strain, states present/genome_only/absent),
frequency_table (family, frequency, strain_count, bin), cooccurring_pairs
(family_a, ..., direction_a), outgroup strain names.
"""
import argparse
import collections
import subprocess
import sys


def open_maybe_zst(path):
    if path.endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        return p.stdout
    return open(path)


ap = argparse.ArgumentParser()
ap.add_argument("--matrix", required=True)
ap.add_argument("--freq", required=True)
ap.add_argument("--pairs", required=True)
ap.add_argument("--outgroup", nargs="+", required=True)
args = ap.parse_args()

fam_bin = {}
with open(args.freq) as fh:
    next(fh)
    for line in fh:
        f, _, _, b = line.rstrip("\n").split("\t")
        fam_bin[f] = b

# outgroup states per family, plus genome_only totals per strain
out_state = {}
go_per_strain = collections.Counter()
with open(args.matrix) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {s: i for i, s in enumerate(header)}
    oi = [idx[s] for s in args.outgroup]
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        out_state[parts[0]] = tuple(parts[i] for i in oi)
        for s, i in zip(args.outgroup, oi):
            if parts[i] == "genome_only":
                go_per_strain[s] += 1

# direction per family_a (first occurrence; direction_a depends only on family_a)
fam_dir = {}
rows = collections.Counter()
with open_maybe_zst(args.pairs) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    c = header.index("direction_a")
    for line in fh:
        parts = line.rstrip("\n").split("\t", c + 1)
        d = parts[c]
        rows[d] += 1
        if parts[0] not in fam_dir:
            fam_dir[parts[0]] = d

tot = sum(rows.values())
print("pair rows: " + ", ".join(f"{k} {v} ({v / tot * 100:.2f}%)" for k, v in sorted(rows.items())))
print("genome_only (rescued) cells in outgroup columns: " + ", ".join(f"{s} {go_per_strain[s]}" for s in args.outgroup))
print()
print("bin\tn_fam\tout_all_present%\tout_any_present%\tout_any_genome_only%\tn_polarised\tgain\tloss\tambiguous")
bins = ["core", "soft_core", "shell", "cloud", "singleton"]
for b in bins:
    fams = [f for f, fb in fam_bin.items() if fb == b and f in out_state]
    n = len(fams)
    if not n:
        continue
    allp = sum(all(s in ("present", "genome_only") for s in out_state[f]) for f in fams)
    anyp = sum(any(s in ("present", "genome_only") for s in out_state[f]) for f in fams)
    anygo = sum(any(s == "genome_only" for s in out_state[f]) for f in fams)
    d = collections.Counter(fam_dir[f] for f in fams if f in fam_dir)
    npol = sum(d.values())
    print(f"{b}\t{n}\t{allp / n * 100:.1f}\t{anyp / n * 100:.1f}\t{anygo / n * 100:.1f}\t{npol}\t"
          f"{d['gain']}\t{d['loss']}\t{d['ambiguous']}")
