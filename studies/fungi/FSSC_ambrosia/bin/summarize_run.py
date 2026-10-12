#!/usr/bin/env python3
"""Print measured summary numbers for one pangenome.nf output directory.

Usage: summarize_run.py <.../output/pangenome> [<samplesheet.csv>]

Reads only pipeline outputs. Reports: family count; frequency-bin counts
(frequency_table.tsv); families by group occupancy (presence_matrix.tsv or
presence_matrix.rescued.tsv, IN/OUT from the samplesheet); significant
islands; pair-classification counts; Leiden module count; assembly-quality
correlations; openness fit.
"""
from __future__ import annotations

import csv
import io
import subprocess
import sys
from collections import Counter
from pathlib import Path


def rows(path: Path):
    if path.suffix == ".zst":
        p = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, text=True, check=True)
        return list(csv.DictReader(io.StringIO(p.stdout), delimiter="\t"))
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    d = Path(sys.argv[1])
    sheet = Path(sys.argv[2]) if len(sys.argv) > 2 else d / "samplesheet.with_clades.csv"
    group = {r["Short"]: r["GROUP"] for r in csv.DictReader(open(sheet))}

    ft = rows(d / "frequency_table.tsv")
    print(f"families_in_frequency_table\t{len(ft)}")
    for b, n in sorted(Counter(r["bin"] for r in ft).items()):
        print(f"bin_{b}\t{n}")
    # the pipeline's 'singleton' bin also holds families with 0 ingroup carriers
    print(f"bin_singleton_ingroup_count0\t{sum(1 for r in ft if r['bin'] == 'singleton' and int(r['strain_count']) == 0)}")
    print(f"bin_singleton_ingroup_count1\t{sum(1 for r in ft if r['bin'] == 'singleton' and int(r['strain_count']) == 1)}")

    mfile = d / "presence_matrix.rescued.tsv"
    if not mfile.exists():
        mfile = d / "presence_matrix.tsv"
    occ = Counter()
    per_strain = Counter()
    with open(mfile) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")[1:]
        g = [group.get(s, "?") for s in hdr]
        n_in, n_out = g.count("IN"), g.count("OUT")
        nf = 0
        for line in fh:
            nf += 1
            v = line.rstrip("\n").split("\t")[1:]
            pin = pout = 0
            for s, gg, x in zip(hdr, g, v):
                if x != "absent":
                    per_strain[s] += 1
                    if gg == "IN":
                        pin += 1
                    elif gg == "OUT":
                        pout += 1
            key = ("IN_only" if pin and not pout else "OUT_only" if pout and not pin
                   else "both" if pin and pout else "none")
            occ[key] += 1
            if pin == n_in and pout == n_out:
                occ["present_in_all_strains"] += 1
            if pin >= 0.95 * n_in and pout >= 0.95 * n_out:
                occ["present_ge95pct_both_groups"] += 1
    print(f"matrix_file\t{mfile.name}\tfamilies\t{nf}\tIN_strains\t{n_in}\tOUT_strains\t{n_out}")
    for k in ("IN_only", "OUT_only", "both", "none", "present_in_all_strains", "present_ge95pct_both_groups"):
        print(f"families_{k}\t{occ[k]}")
    for gg in ("IN", "OUT"):
        vals = [per_strain[s] for s, x in zip(hdr, g) if x == gg]
        if vals:
            vals.sort()
            print(f"families_per_strain_{gg}\tmin {vals[0]}\tmedian {vals[len(vals)//2]}\tmax {vals[-1]}")

    si = d / "significant_islands.tsv"
    if si.exists():
        isl = rows(si)
        print(f"significant_islands\t{len(isl)}")
    pc = d / "pair_classification.tsv.zst"
    if pc.exists():
        p = subprocess.Popen(["zstd", "-dc", str(pc)], stdout=subprocess.PIPE, text=True)
        c = Counter()
        hdrp = p.stdout.readline().rstrip("\n").split("\t")
        ci = hdrp.index("classification")
        for line in p.stdout:
            c[line.split("\t")[ci]] += 1
        p.wait()
        for k, n in c.most_common():
            print(f"pairs_{k}\t{n}")
    ms = d / "module_summary.tsv"
    if ms.exists():
        m = rows(ms)
        sizes = sorted((int(r["size"]) for r in m), reverse=True)
        print(f"leiden_modules\t{len(m)}\tlargest\t{sizes[:5]}")
    for name in ("assembly_quality_correlations.tsv", "report/pangenome_openness.tsv"):
        f = d / name
        if f.exists():
            print(f"--- {name}")
            print(f.read_text().rstrip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
