#!/usr/bin/env python3
"""Feasibility measurements for an exemplar-anchored island locus view.

Measures, on one pangenome.nf run:
  M1  exemplar flank availability: genes on the island's contig left/right of it
  M2  flank composition: bins of the N nearest flank genes each side
  M3  per-strain column states for the top islands, under the proposed rule:
        syntenic  - strain has a copy on some contig X, and a copy of at least one
                    other column family (any other column) on X within K ranks
        elsewhere - present (annotated or rescued), no such neighbour
        absent
      and flank-anchor status per strain (left flank block and right flank
      block both syntenic with each other, i.e. the locus is intact around the
      island).
  M4  island redundancy: located islands whose member set is >= 50% contained in a
      larger located island with the same example-strain contig.

Usage (pixi python or python3.12):
  island_locus_view_feasibility.py --run_dir <.../output/pangenome> --flank 5 --k 10 --top 50
"""
from __future__ import annotations

import argparse
import collections
import csv
import statistics
import subprocess
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".zst"):
        return subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True).stdout
    return open(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", required=True, type=Path)
    ap.add_argument("--flank", type=int, default=5)
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--min_strains", type=int, default=2)
    ap.add_argument("--sort", choices=["size", "strains"], default="size",
                    help="rank qualifying islands by island_size (the page's order) or n_strains")
    ap.add_argument("--empty_frac", type=float, default=1.0,
                    help="island block counts as an empty site when >= this fraction of its "
                    "columns are absent (1.0 = all absent)")
    a = ap.parse_args()
    R = a.run_dir

    bins = {}
    with open(R / "frequency_table.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bins[r["family"]] = r["bin"]

    islands = []
    with open(R / "report_tables" / "islands_with_domains.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["locus_id"] in ("-", ""):
                continue
            r["members"] = [f for f in r["member_families"].split(",") if f]
            islands.append(r)
    qual = [r for r in islands if int(r["n_strains"]) >= a.min_strains]
    if a.sort == "size":
        qual.sort(key=lambda r: -int(r["island_size"]))
    else:
        qual.sort(key=lambda r: (-int(r["n_strains"]), -int(r["island_size"])))
    top = qual[:a.top]
    print(f"located islands {len(islands)}; >= {a.min_strains} strains {len(qual)}; top {len(top)}")

    # M4: redundancy among located islands (same example contig not required for
    # member overlap; report both).
    by_size = sorted(islands, key=lambda r: -len(r["members"]))
    sets = [set(r["members"]) for r in by_size]
    fam_index = collections.defaultdict(list)
    for i, s in enumerate(sets):
        for f in s:
            fam_index[f].append(i)
    contained = 0
    for i, s in enumerate(sets):
        cand = collections.Counter(j for f in s for j in fam_index[f] if j < i and len(sets[j]) > len(s))
        if any(c / len(s) >= 0.5 for c in cand.values()):
            contained += 1
    print(f"M4 located islands >= 50% contained in a larger located island: {contained} / {len(sets)} "
          f"({contained / len(sets) * 100:.1f}%)")

    # Pass 1: exemplar strains' full gene order on their locus contigs.
    ex_contig = {(r["example_strain"], r["locus_contig"]) for r in top}
    order = collections.defaultdict(list)  # (strain, contig) -> [(rank, family)]
    with open_text(R / "family_positions.tsv.zst") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["Short"], r["contig"])
            if key in ex_contig:
                order[key].append((int(r["rank"]), r["family"]))
    for v in order.values():
        v.sort()

    # M1/M2 on the top islands; define columns = left flank + island + right flank
    left_avail, right_avail, flank_bins = [], [], collections.Counter()
    columns = {}
    for r in top:
        seq = order[(r["example_strain"], r["locus_contig"])]
        members = set(r["members"])
        idx = [i for i, (_, f) in enumerate(seq) if f in members]
        if not idx:
            continue
        lo, hi = min(idx), max(idx)
        left = [f for _, f in seq[max(0, lo - a.flank):lo]]
        right = [f for _, f in seq[hi + 1:hi + 1 + a.flank]]
        left_avail.append(lo)
        right_avail.append(len(seq) - 1 - hi)
        for f in left + right:
            flank_bins[bins.get(f, "?")] += 1
        columns[r["locus_id"]] = (left, [f for _, f in seq[lo:hi + 1]], right)
    n = len(left_avail)
    both = sum(1 for l, rr in zip(left_avail, right_avail) if l >= a.flank and rr >= a.flank)
    none_side = sum(1 for l, rr in zip(left_avail, right_avail) if l == 0 or rr == 0)
    print(f"M1 top islands with >= {a.flank} genes on BOTH sides on the exemplar contig: {both}/{n}; "
          f"touching a contig end on at least one side: {none_side}/{n}; "
          f"median genes left {statistics.median(left_avail)}, right {statistics.median(right_avail)}")
    tot = sum(flank_bins.values())
    print("M2 flank gene bins: " + ", ".join(f"{b} {c} ({c / tot * 100:.0f}%)" for b, c in flank_bins.most_common()))

    # Pass 2: all strains' copies of every column family
    need = {f for l, isl, rr in columns.values() for f in l + isl + rr}
    pos = collections.defaultdict(list)  # (strain, family) -> [(contig, rank)]
    with open_text(R / "family_positions.tsv.zst") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family"] in need:
                pos[(r["Short"], r["family"])].append((r["contig"], int(r["rank"])))
    strains = sorted({s for s, _ in pos} | {r["example_strain"] for r in top})
    with open(R / "presence_matrix.rescued.tsv") as fh:
        strains = fh.readline().rstrip("\n").split("\t")[1:]

    states = collections.Counter()
    anchored = []  # per island: strains with intact flanks, and among them island-block states
    indel_like = []
    for loc, (left, isl, right) in columns.items():
        cols = left + isl + right
        n_anch = 0
        n_anch_island_all_absent = 0
        n_anch_island_all_syn = 0
        for s in strains:
            st = []
            for i, f in enumerate(cols):
                mine = pos.get((s, f), [])
                if not mine:
                    st.append("absent")
                    continue
                syn = False
                for j, g in enumerate(cols):
                    if j == i:
                        continue
                    for (c1, r1) in mine:
                        for (c2, r2) in pos.get((s, g), []):
                            if c1 == c2 and abs(r1 - r2) <= a.k:
                                syn = True
                                break
                        if syn:
                            break
                    if syn:
                        break
                st.append("syntenic" if syn else "elsewhere")
            for x in st[len(left):len(left) + len(isl)]:
                states[x] += 1
            # flanks intact: some left-flank gene and some right-flank gene on the same
            # contig within (island length + 2K) ranks of each other
            ok = False
            if left and right:
                span = len(isl) + 2 * a.k
                for f in left:
                    for g in right:
                        for (c1, r1) in pos.get((s, f), []):
                            for (c2, r2) in pos.get((s, g), []):
                                if c1 == c2 and abs(r1 - r2) <= span:
                                    ok = True
                                    break
                            if ok:
                                break
                        if ok:
                            break
                    if ok:
                        break
            if ok:
                n_anch += 1
                block = st[len(left):len(left) + len(isl)]
                if block and sum(x == "absent" for x in block) / len(block) >= a.empty_frac:
                    n_anch_island_all_absent += 1
                if all(x == "syntenic" for x in block):
                    n_anch_island_all_syn += 1
        anchored.append(n_anch)
        indel_like.append((n_anch_island_all_absent, n_anch_island_all_syn))
    tot = sum(states.values())
    print(f"M3 island-column cell states over {len(columns)} islands x {len(strains)} strains: " +
          ", ".join(f"{k} {v} ({v / tot * 100:.1f}%)" for k, v in states.most_common()))
    print(f"M3 strains with intact flanks per island: median {statistics.median(anchored)}, "
          f"min {min(anchored)}, max {max(anchored)} (of {len(strains)})")
    ab = [x for x, _ in indel_like]
    sy = [y for _, y in indel_like]
    print(f"M3 among flank-anchored strains, island block >= {a.empty_frac:.0%} absent (empty site): median {statistics.median(ab)}; "
          f"ALL syntenic (full island in place): median {statistics.median(sy)}")
    print(f"M3 islands with >= 10 anchored strains having an empty site AND >= 2 with the full island: "
          f"{sum(1 for x, y in indel_like if x >= 10 and y >= 2)}/{len(indel_like)}")


if __name__ == "__main__":
    main()
