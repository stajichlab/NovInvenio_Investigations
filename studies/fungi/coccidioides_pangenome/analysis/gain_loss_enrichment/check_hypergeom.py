#!/usr/bin/env python3
"""Independent check of bin/pangenome_gain_loss_enrichment.py output.

For --n random Pfam terms per test set (seeded), recompute k, N, K, M directly from the
raw inputs (family_gain_loss.tsv.gz + family annotation table; no code shared with the
enrichment script) and P(X >= k) by exact summation with math.comb / fractions (no scipy).
Compare with enrichment_domains.tsv.gz. Writes a TSV to stdout; exits 1 on any mismatch.

Usage: check_hypergeom.py --gain-loss-table X --family-annotation Y --enrichment Z
           --classes shell,soft_core [--n 5] [--seed 1]
"""
import argparse
import csv
import gzip
import math
import random
import sys
from fractions import Fraction


def read(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def exact_sf(k, M, K, N):
    num = sum(math.comb(K, i) * math.comb(M - K, N - i) for i in range(k, min(K, N) + 1))
    return float(Fraction(num, math.comb(M, N)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gain-loss-table", required=True)
    ap.add_argument("--family-annotation", required=True)
    ap.add_argument("--enrichment", required=True)
    ap.add_argument("--classes", default="shell,soft_core")
    ap.add_argument("--direction-column", default="rule_freq")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    classes = set(a.classes.split(","))
    scope = {r["family"]: r[a.direction_column] for r in read(a.gain_loss_table) if r["class"] in classes}
    pf = {}
    for r in read(a.family_annotation):
        if r["accession"] in scope and r["pfam_ids"]:
            pf[r["accession"]] = {x.split(".")[0] for x in r["pfam_ids"].split("|")}
    M = len(pf)
    enr = [r for r in read(a.enrichment) if r["vocab"] == "pfam"]
    rng = random.Random(a.seed)
    bad = 0
    print("test_set\tterm\tk\tN\tK\tM\tp_script\tp_exact\trel_diff\tcounts_match")
    for ts in sorted({r["test_set"] for r in enr}):
        rows = sorted((r for r in enr if r["test_set"] == ts and int(r["k"]) > 0), key=lambda r: r["term"])
        for r in rng.sample(rows, min(a.n, len(rows))):
            t = r["term"]
            N = sum(1 for f in pf if scope[f] == ts)
            K = sum(1 for f in pf if t in pf[f])
            k = sum(1 for f in pf if t in pf[f] and scope[f] == ts)
            pe = exact_sf(k, M, K, N)
            ps = float(r["p"])
            rel = abs(ps - pe) / pe if pe else abs(ps)
            ok = (k, N, K, M) == (int(r["k"]), int(r["n_set"]), int(r["K"]), int(r["M"]))
            if not ok or rel > 1e-6:
                bad += 1
            print(f"{ts}\t{t}\t{k}\t{N}\t{K}\t{M}\t{ps:.6g}\t{pe:.6g}\t{rel:.2e}\t{ok}")
    print(f"# mismatches: {bad}", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
