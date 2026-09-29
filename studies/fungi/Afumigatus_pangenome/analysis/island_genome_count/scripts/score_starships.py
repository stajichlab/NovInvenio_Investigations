#!/usr/bin/env python3
"""Score one run's significant islands against known A. fumigatus Starships.

Two checks (truth: Gluck-Thaler et al. 2025, mBio, doi:10.1128/mbio.01092-25):

1. Coordinate check, Af293 only (Table S7, 8 elements). An island is "placed"
   in Af293 when at least one member family has an annotated Af293 gene
   (tier1_cluster member "Asfu_Af293|<protein_id>", coordinates from
   gene_positions). Rescued genome_only hits have no gene coordinates and are
   not placed. An island "hits" a Starship when a placed gene's midpoint lies
   inside the element. Baseline: the share of Af293 genes in non-core families
   (this run's frequency_table bin != core) that lie inside any element.
   Islands are built from non-core genes, so this is the matched background.

2. Presence-pattern check, all kept ingroup strains that have a Table S21
   genotype (ground_truth_starships_by_short.tsv). A strain carries an island
   when at least ISLAND_CARRY_FRAC of its member families are present
   (present or genome_only, the pipeline's own is_present rule). For each
   Starship, best Jaccard over all islands. Null: the same statistic after
   shuffling the Starship's presence vector across strains (N_PERM times);
   report the 95th percentile. Starships with fewer than MIN_CLASS carriers
   or non-carriers in the subset are not scored.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

import numpy as np

AF293 = "Asfu_Af293"
ISLAND_CARRY_FRAC = 0.5   # majority of member families present
N_PERM = 200
MIN_CLASS = 3             # carriers and non-carriers needed to score a Starship
PERM_SEED = 20260929


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def zstd_lines(path: Path):
    proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
    yield from proc.stdout
    if proc.wait() != 0:
        raise SystemExit(f"zstd failed on {path}")


def load_af293_genes(gene_positions: Path) -> dict[str, tuple[str, int]]:
    """protein_id -> (contig, midpoint) for Af293 genes."""
    genes = {}
    lines = zstd_lines(gene_positions)
    header = next(lines).rstrip("\n").split("\t")
    if header[:5] != ["Short", "protein_id", "contig", "start", "end"]:
        raise SystemExit(f"unexpected gene_positions header {header}")
    for line in lines:
        short, pid, contig, start, end = line.rstrip("\n").split("\t")[:5]
        if short == AF293:
            genes[pid] = (contig, (int(start) + int(end)) // 2)
    return genes


def load_af293_members(cluster_tsv: Path) -> dict[str, list[str]]:
    """family rep -> Af293 protein_ids (prefix stripped)."""
    prefix = AF293 + "|"
    fam = {}
    with cluster_tsv.open() as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            if member.startswith(prefix):
                fam.setdefault(rep, []).append(member[len(prefix):])
    return fam


def in_starship(contig: str, mid: int, ships: list[dict]) -> str | None:
    for s in ships:
        if s["contig"] == contig and s["start"] <= mid <= s["end"]:
            return s["starship_id"]
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--islands", required=True, type=Path)
    ap.add_argument("--frequency_table", required=True, type=Path)
    ap.add_argument("--presence_matrix", required=True, type=Path)
    ap.add_argument("--cluster_tsv", required=True, type=Path)
    ap.add_argument("--gene_positions", required=True, type=Path)
    ap.add_argument("--starships", required=True, type=Path, help="af293_starships_table_s7.tsv")
    ap.add_argument("--presence_truth", required=True, type=Path,
                    help="ground_truth_starships_by_short.tsv")
    ap.add_argument("--subset_strains", required=True, type=Path)
    ap.add_argument("--out_prefix", required=True, type=Path)
    args = ap.parse_args()

    ships = read_tsv(args.starships)
    for s in ships:
        s["start"], s["end"] = int(s["start"]), int(s["end"])
    islands = read_tsv(args.islands)
    for isl in islands:
        isl["families"] = isl["member_families"].split(",")

    # ---- 1. Af293 coordinate check ----
    genes = load_af293_genes(args.gene_positions)
    members = load_af293_members(args.cluster_tsv)
    bins = {r["family"]: r["bin"] for r in read_tsv(args.frequency_table)}

    def gene_hits(fams):
        out = []
        for f in fams:
            for pid in members.get(f, []):
                if pid in genes:
                    contig, mid = genes[pid]
                    out.append(in_starship(contig, mid, ships))
        return out

    noncore = [h for f, b in bins.items() if b not in ("core", "-")
               for h in gene_hits([f])]
    base_frac = (sum(h is not None for h in noncore) / len(noncore)) if noncore else float("nan")

    per_ship = {s["starship_id"]: {"islands": 0} for s in ships}
    n_placed = n_in_ship = 0
    placed_genes = placed_in_ship = 0
    for isl in islands:
        hits = gene_hits(isl["families"])
        if not hits:
            continue
        n_placed += 1
        placed_genes += len(hits)
        placed_in_ship += sum(h is not None for h in hits)
        hit_ships = {h for h in hits if h is not None}
        if hit_ships:
            n_in_ship += 1
        for h in hit_ships:
            per_ship[h]["islands"] += 1
    island_gene_frac = placed_in_ship / placed_genes if placed_genes else float("nan")

    # ---- 2. presence-pattern check ----
    subset = [s for s in args.subset_strains.read_text().split() if s.startswith("Asfu_")]
    truth = {r["nameID"]: json.loads(r["presence_by_short_json"]) for r in read_tsv(args.presence_truth)}
    gt_strains = sorted(s for s in subset if all(s in p for p in truth.values()))
    fam_needed = {f for isl in islands for f in isl["families"]}
    pres: dict[str, np.ndarray] = {}
    with args.presence_matrix.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = [header.index(s) for s in gt_strains]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in fam_needed:
                pres[parts[0]] = np.array([parts[c] != "absent" for c in cols])
    if islands and gt_strains:
        carry = np.array([
            np.mean([pres[f] for f in isl["families"] if f in pres], axis=0) >= ISLAND_CARRY_FRAC
            for isl in islands
        ])
    else:
        carry = np.zeros((0, len(gt_strains)), dtype=bool)

    def best_jaccard(vec: np.ndarray) -> float:
        if carry.shape[0] == 0:
            return 0.0
        c = carry.astype(np.int32)
        v = vec.astype(np.int32)
        inter = c @ v   # int, not bool: bool @ bool returns a logical OR-AND
        union = c.sum(1) + v.sum() - inter
        with np.errstate(invalid="ignore", divide="ignore"):
            jac = np.where(union > 0, inter / union, 0.0)
        return float(jac.max())

    rng = np.random.default_rng(PERM_SEED)
    rows = []
    for s in ships:
        rows.append({"check": "af293_coord", "starship": s["starship_id"],
                     "name": f'{s["name"]} {s["haplotype"]}', "confidence": s["confidence"],
                     "value": per_ship[s["starship_id"]]["islands"], "null95": "", "n_carriers": ""})
    n_scored = n_recovered = 0
    for name, p in truth.items():
        vec = np.array([bool(p[s]) for s in gt_strains])
        k = int(vec.sum())
        if k < MIN_CLASS or len(vec) - k < MIN_CLASS:
            rows.append({"check": "presence_jaccard", "starship": name, "name": name,
                         "confidence": "high", "value": "", "null95": "", "n_carriers": k})
            continue
        obs = best_jaccard(vec)
        null = [best_jaccard(rng.permutation(vec)) for _ in range(N_PERM)]
        q95 = float(np.quantile(null, 0.95))
        n_scored += 1
        n_recovered += obs > q95
        rows.append({"check": "presence_jaccard", "starship": name, "name": name,
                     "confidence": "high", "value": round(obs, 4), "null95": round(q95, 4),
                     "n_carriers": k})

    with open(f"{args.out_prefix}.starships.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    high = [s["starship_id"] for s in ships if s["confidence"] == "high"]
    summary = {
        "n_islands": len(islands),
        "n_islands_placed_af293": n_placed,
        "n_islands_in_af293_starship": n_in_ship,
        "af293_high_conf_starships": len(high),
        "af293_high_conf_starships_hit": sum(per_ship[h]["islands"] > 0 for h in high),
        "island_gene_frac_in_starship": island_gene_frac,
        "noncore_gene_frac_in_starship": base_frac,
        "presence_gt_strains": len(gt_strains),
        "presence_starships_scored": n_scored,
        "presence_starships_above_null95": int(n_recovered),
    }
    Path(f"{args.out_prefix}.summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
