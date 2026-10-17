#!/usr/bin/env python3
"""Coccidioides: antiSMASH SM regions vs DUF3435 (captain) genes.

Everything is in the antiSMASH annotation's own gene space (per-strain funannotate
runs under 2025_All_Cocci/Assembly/annotation; see parse_antismash.py). No Starship
element calls exist for Coccidioides, so 'overlap with a called Starship region' cannot
be computed; the only Starship-related signal is the captain-type tyrosine recombinase
domain DUF3435 (PF11917.14), full-sequence E <= 1e-5. A captain gene is not a called
Starship.

Per run strain set (--config): regions by type, regions that contain >= 1 captain
gene (captain locus tag among the region's CDS), and a pooled hypergeometric test
(universe: all genes of the included strains' annotations; K = captain genes;
n = genes in SM regions; k = captain genes in SM regions).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import hypergeom


def rd(path: Path):
    if str(path).endswith(".zst"):
        out = subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout.decode()
        return out.splitlines(True)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt").readlines()
    return open(path).readlines()


def tsv(path: Path):
    return list(csv.DictReader(rd(path), delimiter="\t"))


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def type_label(products: str) -> str:
    ps = sorted({p for p in products.split(",") if p})
    return ps[0] if len(ps) == 1 else "hybrid:" + "+".join(ps)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--antismash-dir", required=True, type=Path, help="inputs/antismash from parse_antismash.py")
    ap.add_argument("--captain-tblout", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path, help="run samplesheet; restricts strains")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {r["Short"]: r for r in csv.DictReader(open(a.config))}
    regions = [r for r in tsv(a.antismash_dir / "sm_regions.tsv.gz") if r["Short"] in cfg]
    rgenes = [r for r in tsv(a.antismash_dir / "sm_region_genes.tsv.gz") if r["Short"] in cfg]
    counts = {r["Short"]: int(r["n_proteins"]) for r in tsv(a.antismash_dir / "gene_counts.tsv") if r["Short"] in cfg}
    cap = set()
    for line in rd(a.captain_tblout):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if float(p[4]) <= a.evalue:
            short, pid = p[0].split("|", 1)
            if short in cfg:
                cap.add((short, re.sub(r"-T\d+$", "", pid)))
    reg_genes = defaultdict(set)
    for g in rgenes:
        reg_genes[g["region_id"]].add((g["Short"], g["locus_tag"]))
    all_region_genes = set().union(*reg_genes.values()) if reg_genes else set()

    out_rows = []
    for r in regions:
        caps = reg_genes[r["region_id"]] & cap
        out_rows.append({**r, "sm_type": type_label(r["products"]), "species": cfg[r["Short"]]["Species"],
                         "n_captain_genes": len(caps), "captain_genes": ",".join(sorted(t for _, t in caps)),
                         "starship_overlap": "NA (no Starship calls)"})
    with gzip.open(a.out_dir / "sm_clusters.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    by = defaultdict(Counter)
    for r in out_rows:
        for key in (("all", r["sm_type"]), (r["species"], r["sm_type"])):
            by[key]["n_clusters"] += 1
            by[key]["with_captain"] += r["n_captain_genes"] > 0
    with open(a.out_dir / "sm_by_type.tsv", "w") as fh:
        fh.write("species\tsm_type\tn_clusters\tn_with_captain_gene\tpct_with_captain_gene\n")
        for (sp, t), c in sorted(by.items(), key=lambda kv: (kv[0][0] != "all", kv[0][0], -kv[1]["n_clusters"])):
            fh.write(f"{sp}\t{t}\t{c['n_clusters']}\t{c['with_captain']}\t{100 * c['with_captain'] / c['n_clusters']:.2f}\n")

    N = sum(counts.values())
    K = len(cap)
    n = len(all_region_genes)
    k = len(all_region_genes & cap)
    p = hypergeom.sf(k - 1, N, K, n)
    strains_ann = sorted(counts)
    m = [
        ("strains_in_run", len(cfg), "rows of the run samplesheet"),
        ("strains_with_antismash", len(strains_ann), "exact name match to an annotation dir with antiSMASH output"),
        ("strains_without_antismash", ",".join(sorted(set(cfg) - set(counts))), ""),
        ("sm_clusters", len(out_rows), "antiSMASH 8.0.1 regions"),
        ("sm_clusters_with_captain_gene", sum(r["n_captain_genes"] > 0 for r in out_rows), "region CDS include a DUF3435 hit"),
        ("sm_clusters_overlapping_starship_ge1gene", "NA", "no Starship element calls for Coccidioides"),
        ("sm_clusters_overlapping_starship_ge50pct", "NA", "no Starship element calls for Coccidioides"),
        ("genes_total", N, "proteins in the included strains' annotations"),
        ("captain_genes", K, f"DUF3435 hits, E <= {a.evalue}"),
        ("genes_in_sm_clusters", n, ""),
        ("captain_genes_in_sm_clusters", k, ""),
        ("expected_captain_genes_in_sm_clusters", round(K * n / N, 2) if N else "", "K*n/N"),
        ("p_hypergeom_greater", f"{p:.3g}", "pooled over strains"),
        ("p_hypergeom_less", f"{hypergeom.cdf(k, N, K, n):.3g}", "pooled over strains (depletion)"),
    ]
    with open(a.out_dir / "sm_starship_overlap.tsv", "w") as fh:
        fh.write("metric\tvalue\tdefinition\n")
        for row in m:
            fh.write("\t".join(map(str, row)) + "\n")
    with open(a.out_dir / "PROVENANCE_sm.yaml", "w") as fh:
        fh.write(f"derived_by: {Path(__file__).resolve()}\ndate: {dt.date.today().isoformat()}\n"
                 f"command: {' '.join(sys.argv)}\nsource_annotation: /bigdata/stajichlab/shared/projects/Coccidioides/"
                 "PopGenomics/2025_All_Cocci/Assembly/annotation/<strain>/{antismash_local,annotate_misc,annotate_results}"
                 " (antiSMASH 8.0.1, funannotate 1.8.17; not the 20260112 freeze used by the pangenome runs)\ninputs:\n")
        for nm in ("sm_regions.tsv.gz", "sm_region_genes.tsv.gz", "gene_counts.tsv"):
            fh.write(f"  {nm}:\n    path: {(a.antismash_dir / nm).resolve()}\n    sha256: {sha(a.antismash_dir / nm)}\n")
        for nm, pth in (("captain_tblout", a.captain_tblout), ("config", a.config)):
            fh.write(f"  {nm}:\n    path: {Path(pth).resolve()}\n    sha256: {sha(pth)}\n")


if __name__ == "__main__":
    main()
