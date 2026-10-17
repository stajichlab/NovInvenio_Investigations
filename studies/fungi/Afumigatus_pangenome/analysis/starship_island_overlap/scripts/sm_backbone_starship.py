#!/usr/bin/env python3
"""A. fumigatus: secondary-metabolite (SM) backbone genes vs Af293 Starships and islands.

No antiSMASH result exists for the proteins used by the pangenome run (NCBI/UniProt
annotation). The SM proxy is the existing hmmsearch of two backbone Pfam models
(results/sm_backbone/SM_backbone.hmm = PF00109 ketoacyl-synt [PKS KS] + PF00668
Condensation [NRPS C]) against all 2,788,402 study proteins
(results/full_293run/all_strains.fa). Hits with full-sequence E <= 1e-5 count.
So 'SM' here means PKS or NRPS backbone genes only; terpene, RiPP, DMATS and other
classes are not covered, and there are no cluster boundaries.

Starship coordinates exist for Af293 only (Table S7), so the gene-level overlap is
for Af293 genes. The family-level view uses all strains.

Outputs (--out-dir): sm_starship_overlap.tsv (metric table), sm_by_type.tsv,
sm_backbone_af293_genes.tsv.gz, sm_backbone_families.tsv.gz, PROVENANCE.yaml.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import sys
from collections import defaultdict
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "bin"))
import starship_island_overlap as sio  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

TYPE = {"ketoacyl-synt": "PKS (KS)", "Condensation": "NRPS (C)"}


def tblout_hits(path: Path, evalue: float) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with sio.open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if float(p[4]) <= evalue:
                out[p[0]].add(TYPE.get(p[2], p[2]))
    return out


def hyper(N, K, n, k):
    exp = K * n / N
    return {"N": N, "K": K, "n": n, "k": k, "expected": round(exp, 2),
            "fold": round(k / exp, 2) if exp else "", "p_greater": f"{hypergeom.sf(k - 1, N, K, n):.3g}"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--starship-annotation", required=True, type=Path)
    ap.add_argument("--sm-tblout", required=True, type=Path)
    ap.add_argument("--captain-tblout", required=True, type=Path)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--extra-gene-positions", type=Path, help="coordinates for unplaced Af293 members")
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rd = a.run_dir
    freq = sio.load_frequency(rd / "frequency_table.tsv")
    islands = sio.load_islands(rd / "significant_islands.tsv")
    island_fams = set().union(*(i["family_set"] for i in islands)) & set(freq)
    els = sio.load_elements(a.starship_annotation)
    strains = {e["Short"] for e in els}
    all_m2r, strain_m2r = sio.load_members(rd / "cluster" / "tier1_cluster.tsv", strains)
    genes = sio.load_strain_genes(rd / "gene_positions.tsv.zst", strains)
    if a.extra_gene_positions:
        have = {(g[0], g[1]) for g in genes}
        genes += [g for g in sio.load_strain_genes(a.extra_gene_positions, strains) if (g[0], g[1]) not in have]
    sm = tblout_hits(a.sm_tblout, a.evalue)
    cap = tblout_hits(a.captain_tblout, a.evalue)
    efams, egenes = sio.starship_families(genes, strain_m2r, els)
    gene_el = {(g[0], g[1]): eid for eid, gs in egenes.items() for g in gs}
    star_fams = set().union(*efams.values()) & set(freq)
    cap_fams = {all_m2r[m] for m in cap if m in all_m2r}

    # ---- Af293 gene level
    rows, metrics = [], []
    n_all = len(genes)
    in_el = {k for k in gene_el}
    sm_genes = [(g, sorted(sm[f"{g[0]}|{g[1]}"])) for g in genes if f"{g[0]}|{g[1]}" in sm]
    cap_genes = {(g[0], g[1]) for g in genes if f"{g[0]}|{g[1]}" in cap}
    for g, types in sm_genes:
        fam = strain_m2r.get(f"{g[0]}|{g[1]}")
        rows.append({"Short": g[0], "protein_id": g[1], "contig": g[2], "start": g[3], "end": g[4],
                     "sm_type": "+".join(types), "starship_element": gene_el.get((g[0], g[1]), ""),
                     "family": fam or "", "bin": freq.get(fam, {}).get("bin", ""),
                     "family_in_island": "yes" if fam in island_fams else "no",
                     "family_has_captain_hit": "yes" if fam in cap_fams else "no"})
    sio.write_tsv(a.out_dir / "sm_backbone_af293_genes.tsv.gz", rows)
    k = sum(1 for g, _ in sm_genes if (g[0], g[1]) in in_el)

    def M(metric, value, note=""):
        metrics.append({"metric": metric, "value": value, "definition": note})

    M("af293_genes", n_all, "Af293 rows of gene_positions" + (" + --extra-gene-positions" if a.extra_gene_positions else ""))
    M("af293_genes_in_starships", len(in_el), "gene midpoint inside a Table S7 element (8 elements)")
    M("af293_captain_genes", len(cap_genes), "DUF3435 hit, E <= 1e-5")
    M("af293_captain_genes_in_starships", len(cap_genes & in_el), "")
    M("af293_sm_backbone_genes", len(sm_genes), "PKS KS or NRPS C hit, E <= 1e-5")
    M("af293_sm_backbone_genes_in_starships", k, "")
    h = hyper(n_all, len(in_el), len(sm_genes), k)
    M("hypergeom_af293_genes", h, "universe: all Af293 genes in gene_positions")
    nc = [g for g in genes if freq.get(strain_m2r.get(f"{g[0]}|{g[1]}"), {}).get("bin", "core") not in sio.CORE_BINS | {"-"}]
    nc_keys = {(g[0], g[1]) for g in nc}
    h2 = hyper(len(nc), len(in_el & nc_keys), sum(1 for g, _ in sm_genes if (g[0], g[1]) in nc_keys),
               sum(1 for g, _ in sm_genes if (g[0], g[1]) in nc_keys and (g[0], g[1]) in in_el))
    M("hypergeom_af293_noncore_genes", h2, "universe: Af293 genes of non-core families")
    ships_with_sm = sorted({gene_el[(g[0], g[1])] for g, _ in sm_genes if (g[0], g[1]) in in_el})
    M("starships_with_sm_backbone_gene", f"{len(ships_with_sm)} of {len(els)}", ",".join(ships_with_sm))

    # ---- family level (all strains)
    sm_fams_by_type: dict[str, set[str]] = defaultdict(set)
    for m, types in sm.items():
        fam = all_m2r.get(m)
        if fam in freq:
            for t in types:
                sm_fams_by_type[t].add(fam)
    sm_fams = set().union(*sm_fams_by_type.values())
    fam_rows = []
    for f in sorted(sm_fams):
        fam_rows.append({"family": f, "sm_types": "+".join(sorted(t for t, fs in sm_fams_by_type.items() if f in fs)),
                         "bin": freq[f]["bin"], "strain_count": freq[f]["strain_count"],
                         "starship_family": "yes" if f in star_fams else "no",
                         "captain_family": "yes" if f in cap_fams else "no",
                         "in_island": "yes" if f in island_fams else "no"})
    sio.write_tsv(a.out_dir / "sm_backbone_families.tsv.gz", fam_rows)
    M("sm_backbone_families", len(sm_fams), "families with >= 1 member hit (any strain)")
    M("sm_backbone_families_starship", len(sm_fams & star_fams), "also a Starship family (Af293 gene in an element)")
    M("sm_backbone_families_captain", len(sm_fams & cap_fams), "also a DUF3435 family")
    M("sm_backbone_families_in_island", len(sm_fams & island_fams), "")
    nonc_f = {f for f, r in freq.items() if r["bin"] not in sio.CORE_BINS | {"-"}}
    M("sm_backbone_families_noncore", len(sm_fams & nonc_f), "")
    M("hypergeom_family_island_noncore", hyper(len(nonc_f), len(sm_fams & nonc_f), len(island_fams & nonc_f),
                                               len(sm_fams & nonc_f & island_fams)),
      "universe: non-core families; K = SM backbone families; n = island families")
    sio.write_tsv(a.out_dir / "sm_starship_overlap.tsv", metrics)

    by_type = []
    for t in sorted(sm_fams_by_type):
        ag = [r for r in rows if t in r["sm_type"]]
        fs = sm_fams_by_type[t]
        by_type.append({"sm_type": t, "af293_genes": len(ag),
                        "af293_genes_in_starship": sum(1 for r in ag if r["starship_element"]),
                        "families_all_strains": len(fs), "families_starship": len(fs & star_fams),
                        "families_in_island": len(fs & island_fams), "families_noncore": len(fs & nonc_f)})
    sio.write_tsv(a.out_dir / "sm_by_type.tsv", by_type)
    with open(a.out_dir / "PROVENANCE_sm.yaml", "w") as fh:
        fh.write(f"derived_by: {Path(__file__).resolve()}\ndate: {dt.date.today().isoformat()}\n"
                 f"command: {' '.join(sys.argv)}\ninputs:\n")
        for k_, p in (("run_dir_frequency_table", rd / "frequency_table.tsv"),
                      ("significant_islands", rd / "significant_islands.tsv"),
                      ("starship_annotation", a.starship_annotation), ("sm_tblout", a.sm_tblout),
                      ("captain_tblout", a.captain_tblout)):
            fh.write(f"  {k_}:\n    path: {Path(p).resolve()}\n    sha256: {sio.sha256(Path(p))}\n")


if __name__ == "__main__":
    main()
