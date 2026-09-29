#!/usr/bin/env python3
"""Collect outputs/runs/n<N>_s<seed>/ into outputs/sweep_summary.tsv (one row
per run) and outputs/sweep_by_n.tsv (median and range over seeds per N)."""
from __future__ import annotations

import csv
import json
import re
import statistics
from pathlib import Path

AN = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/"
          "studies/fungi/Afumigatus_pangenome/analysis/island_genome_count")
METRICS = ["n_cooccurring_pairs", "n_islands", "n_islands_placed_af293",
           "n_islands_in_af293_starship", "af293_high_conf_starships_hit",
           "island_gene_frac_in_starship", "noncore_gene_frac_in_starship",
           "presence_starships_scored", "presence_starships_above_null95",
           "runtime_s_island_steps"]


def main() -> None:
    rows = []
    for d in sorted((AN / "outputs/runs").glob("n*_s*")):
        m = re.fullmatch(r"n(\d+)_s(\d+)", d.name)
        summ = d / "starship.summary.json"
        if not m or not summ.exists():
            continue
        row = {"n": int(m[1]), "seed": int(m[2]), **json.loads(summ.read_text())}
        row["n_cooccurring_pairs"] = int((d / "n_cooccurring_pairs.txt").read_text())
        t = dict(line.split("\t") for line in (d / "timings.tsv").read_text().splitlines())
        row["runtime_s_island_steps"] = sum(
            int(t[k]) for k in ("frequency_bins", "cooccurrence", "pair_classification", "build_islands"))
        rows.append(row)
    rows.sort(key=lambda r: (r["n"], r["seed"]))
    cols = ["n", "seed"] + METRICS
    with (AN / "outputs/sweep_summary.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    by_n = {}
    for r in rows:
        by_n.setdefault(r["n"], []).append(r)
    with (AN / "outputs/sweep_by_n.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["n", "n_seeds"] + [f"{m}_{s}" for m in METRICS for s in ("median", "min", "max")])
        for n, rs in sorted(by_n.items()):
            out = [n, len(rs)]
            for m in METRICS:
                vals = [r[m] for r in rs]
                out += [round(statistics.median(vals), 4), round(min(vals), 4), round(max(vals), 4)]
            w.writerow(out)
    print(f"{len(rows)} runs summarized")


if __name__ == "__main__":
    main()
