#!/usr/bin/env python3
"""A. fumigatus Af293: antiSMASH SM regions vs Table S7 Starship elements and pangenome islands.

Cluster-level extension of ../../scripts/sm_backbone_starship.py (which had only
PKS KS / NRPS C backbone genes and no cluster boundaries).

Inputs
  --antismash-json      antiSMASH 8 JSON (plain or .gz) run on the study's Af293 genome
                        with the study GFF3 gene models (CDS names = GFF3 locus_tag).
  --gff3                the same GFF3 (gene features: locus_tag -> contig/start/end).
  --pep                 the study's Af293 protein FASTA (UniProt; GN= locus tag fallback).
  --genome              the genome FASTA given to antiSMASH (chromosome lengths).
  --run-dir             pangenome.nf output dir (frequency_table, significant_islands,
                        cluster/tier1_cluster.tsv, gene_positions.tsv.zst).
  --extra-gene-positions  ../inputs/af293_unplaced_swissprot_coords.tsv (has locus_tag).
  --starship-annotation ../inputs/af293_starships_table_s7.elements.tsv
  --captain-tblout      DUF3435 hmmsearch tblout (targets 'Short|protein_id').
  --gain-loss           family_gain_loss.tsv.gz (column 'direction').
  --proxy-genes         sm_backbone_af293_genes.tsv.gz of the proxy run (comparison only).

Definitions
  SM cluster            one antiSMASH region. type = region products joined by '+'.
  cluster genes         CDS features inside the region (antiSMASH CDS coordinates).
  core gene             CDS with >= 1 definition domain for a protocluster rule
                        (hmm_detection rule_results; antiSMASH 'biosynthetic' kind).
  gene inside element   CDS midpoint (start+end)//2 inside a Table S7 element (1-based
                        inclusive), the rule of bin/starship_island_overlap.py.
  overlap rule A        >= 1 cluster gene inside an element.
  overlap rule B        >= 50% of the cluster genes inside one element set.
  captain gene          Af293 gene whose pangenome member has a DUF3435 hit, E <= 1e-5.
  locus tag -> member   (1) exact (contig, start, end) match of the GFF3 gene to a
                        gene_positions / --extra-gene-positions row (extra rows carry the
                        locus tag), (2) else GN=<locus tag> in the protein FASTA header.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "bin"))
import starship_island_overlap as sio  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

SHORT = "Asfu_Af293"


def loc_span(loc: str) -> tuple[int, int]:
    """antiSMASH location string -> 1-based inclusive (start, end)."""
    nums = [int(v) for pair in re.findall(r"(\d+):(\d+)", loc) for v in pair]
    return min(nums) + 1, max(nums)


def fasta_lengths(path: Path) -> dict[str, int]:
    out, cur = {}, None
    with sio.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                out[cur] = 0
            elif cur:
                out[cur] += len(line.strip())
    return out


def hyper(N, K, n, k):
    exp = K * n / N if N else float("nan")
    return {"N": N, "K": K, "n": n, "k": k, "expected": round(exp, 2),
            "fold": round(k / exp, 2) if exp else "", "p_greater": f"{hypergeom.sf(k - 1, N, K, n):.3g}"}


def load_json(path: Path):
    with sio.open_text(path) as fh:
        return json.load(fh)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for n in ("antismash-json", "gff3", "pep", "genome", "run-dir", "extra-gene-positions",
              "starship-annotation", "captain-tblout", "gain-loss", "proxy-genes", "out-dir"):
        ap.add_argument(f"--{n}", required=True, type=Path)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261017)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rd = a.run_dir

    # ---------------- pangenome tables
    freq = sio.load_frequency(rd / "frequency_table.tsv")
    islands = sio.load_islands(rd / "significant_islands.tsv")
    island_fams = set().union(*(i["family_set"] for i in islands)) & set(freq)
    gl = sio.load_gain_loss(a.gain_loss)
    all_m2r, strain_m2r = sio.load_members(rd / "cluster" / "tier1_cluster.tsv", {SHORT})
    els = sio.load_elements(a.starship_annotation)
    nonc_f = {f for f, r in freq.items() if r["bin"] not in sio.CORE_BINS | {"-"}}

    # ---------------- locus tag -> member
    gff_gene = {}
    with open(a.gff3) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 9 and p[2] == "gene":
                m = re.search(r"locus_tag=([^;]+)", p[8])
                if m:
                    gff_gene[m.group(1)] = (p[0], int(p[3]), int(p[4]))
    by_coord = {}
    for g in sio.load_strain_genes(rd / "gene_positions.tsv.zst", {SHORT}):
        by_coord.setdefault((g[2], g[3], g[4]), g[1])
    extra_lt = {}
    with open(a.extra_gene_positions) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Short"] == SHORT and r.get("locus_tag"):
                extra_lt[r["locus_tag"]] = r["protein_id"]
                by_coord.setdefault((r["contig"], int(r["start"]), int(r["end"])), r["protein_id"])
    gn = {}
    with sio.open_text(a.pep) as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.search(r" GN=(\S+)", line)
                if m:
                    gn.setdefault(m.group(1), line[1:].split()[0])

    def member_of(lt: str) -> tuple[str, str]:
        if lt in gff_gene and gff_gene[lt] in by_coord:
            return by_coord[gff_gene[lt]], "coord"
        if lt in extra_lt:
            return extra_lt[lt], "extra_locus_tag"
        if lt in gn:
            return gn[lt], "fasta_GN"
        return "", "none"

    # ---------------- captain hits (members)
    cap_members = set()
    with sio.open_text(a.captain_tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if float(p[4]) <= a.evalue and p[0].startswith(SHORT + "|"):
                cap_members.add(p[0][len(SHORT) + 1:])

    # ---------------- antiSMASH JSON
    js = load_json(a.antismash_json)
    cds = {}          # locus -> dict
    regions = []
    for rec in js["records"]:
        contig = rec["id"]
        for f in rec["features"]:
            if f["type"] == "CDS":
                q = f["qualifiers"]
                name = (q.get("locus_tag") or q.get("gene") or q.get("protein_id") or [""])[0]
                s, e = loc_span(f["location"])
                cds[name] = {"contig": contig, "start": s, "end": e, "mid": (s + e) // 2}
        core = set()
        hm = rec["modules"].get("antismash.detection.hmm_detection", {}).get("rule_results", {})
        for _pc, hits in hm.get("cds_by_protocluster", []):
            for h in hits:
                if any(h.get("definition_domains", {}).values()):
                    core.add(h["cds_name"])
        kcb = {}
        cbm = rec["modules"].get("antismash.modules.clusterblast", {})
        for res in (cbm.get("knowncluster") or {}).get("results", []):
            ranking = res.get("ranking") or []
            if ranking:
                ref, score = ranking[0]
                kcb[res.get("region_number")] = (ref.get("accession", ""), ref.get("description", ""),
                                                 score.get("similarity", ""))
        rnum = 0
        for f in rec["features"]:
            if f["type"] != "region":
                continue
            rnum += 1
            q = f["qualifiers"]
            s, e = loc_span(f["location"])
            num = int(q.get("region_number", [rnum])[0])
            regions.append({"contig": contig, "start": s, "end": e, "region_number": num,
                            "products": q.get("product", []), "contig_edge": q.get("contig_edge", [""])[0],
                            "core_names": core, "kcb": kcb.get(num)})
    for r in regions:
        r["genes"] = sorted((n for n, c in cds.items() if c["contig"] == r["contig"] and r["start"] <= c["start"] and c["end"] <= r["end"]),
                            key=lambda n: cds[n]["start"])
        r["core"] = [g for g in r["genes"] if g in r["core_names"]]

    # ---------------- elements
    els_by_c = defaultdict(list)
    for e in els:
        els_by_c[e["contig"]].append(e)

    def el_of(contig, mid, placed=None):
        for e in (placed or els_by_c).get(contig, ()):
            if e["start"] <= mid <= e["end"]:
                return e["element_id"]
        return None

    for n, c in cds.items():
        c["element"] = el_of(c["contig"], c["mid"]) or ""
        c["member"], c["map"] = member_of(n)
        c["family"] = strain_m2r.get(f"{SHORT}|{c['member']}", "") if c["member"] else ""
        c["captain"] = c["member"] in cap_members

    # ---------------- per cluster
    rows = []
    for r in sorted(regions, key=lambda r: (r["contig"], r["start"])):
        g = r["genes"]
        ins = [x for x in g if cds[x]["element"]]
        frac = len(ins) / len(g) if g else 0.0
        r.update(n_in=len(ins), frac=frac, ruleA=len(ins) >= 1, ruleB=frac >= 0.5,
                 elements=sorted({cds[x]["element"] for x in ins}),
                 captain=[x for x in g if cds[x]["captain"]],
                 type="+".join(r["products"]),
                 cid=f"{r['contig']}.region{r['region_number']:03d}")
        kc = r["kcb"] or ("", "", "")
        rows.append({"cluster_id": r["cid"], "contig": r["contig"], "start": r["start"], "end": r["end"],
                     "length_bp": r["end"] - r["start"] + 1, "type": r["type"], "contig_edge": r["contig_edge"],
                     "n_genes": len(g), "n_core_genes": len(r["core"]),
                     "core_genes": ",".join(r["core"]),
                     "core_members": ",".join(cds[x]["member"] for x in r["core"]),
                     "core_families": ",".join(cds[x]["family"] for x in r["core"]),
                     "n_genes_in_starship": len(ins), "frac_genes_in_starship": round(frac, 4),
                     "starship_elements": ",".join(r["elements"]),
                     "overlap_ge1gene": "yes" if r["ruleA"] else "no",
                     "overlap_ge50pct": "yes" if r["ruleB"] else "no",
                     "captain_genes": ",".join(r["captain"]),
                     "knowncluster_top_hit": kc[0], "knowncluster_description": kc[1],
                     "knowncluster_similarity": kc[2],
                     "genes": ",".join(g), "members": ",".join(cds[x]["member"] for x in g)})
    sio.write_tsv(a.out_dir / "sm_clusters.tsv.gz", rows)

    # ---------------- per type
    by_type = []
    types = sorted({r["type"] for r in regions}, key=lambda t: (-sum(1 for r in regions if r["type"] == t), t))
    for t in types + ["ALL"]:
        rs = [r for r in regions if t == "ALL" or r["type"] == t]
        by_type.append({"type": t, "n_clusters": len(rs), "n_overlap_ge1gene": sum(r["ruleA"] for r in rs),
                        "n_overlap_ge50pct": sum(r["ruleB"] for r in rs),
                        "n_with_captain": sum(1 for r in rs if r["captain"])})
    sio.write_tsv(a.out_dir / "sm_by_type.tsv", by_type)

    # ---------------- per element
    per_el = []
    for e in els:
        eg = [n for n, c in cds.items() if c["element"] == e["element_id"]]
        rs = [r for r in regions if e["element_id"] in r["elements"]]
        per_el.append({"element_id": e["element_id"], "label": e["label"], "confidence": e["confidence"],
                       "contig": e["contig"], "start": e["start"], "end": e["end"],
                       "n_genes_in_element": len(eg),
                       "n_sm_clusters_ge1gene": len(rs), "n_sm_clusters_ge50pct": sum(r["ruleB"] for r in rs),
                       "sm_clusters": ";".join(f"{r['cid']}({r['type']};{r['n_in']}/{len(r['genes'])} genes)" for r in rs),
                       "n_sm_cluster_genes_in_element": sum(1 for n in eg if any(n in r["genes"] for r in regions)),
                       "captain_genes_in_element": ",".join(n for n in eg if cds[n]["captain"])})
    sio.write_tsv(a.out_dir / "sm_starship_overlap.tsv", per_el)

    # ---------------- permutation (shift each element on its own chromosome)
    clen = fasta_lengths(a.genome)

    def cluster_stats(placed):
        A = B = gin = 0
        for r in regions:
            k = sum(1 for x in r["genes"] if el_of(cds[x]["contig"], cds[x]["mid"], placed))
            A += k >= 1
            B += bool(r["genes"]) and k / len(r["genes"]) >= 0.5
            gin += k
        return A, B, gin

    obs = cluster_stats(els_by_c)
    rng = random.Random(a.seed)
    null = []
    for _ in range(a.n_perm):
        placed = defaultdict(list)
        for e in els:
            L = e["end"] - e["start"]
            s = rng.randint(1, max(1, clen[e["contig"]] - L))
            placed[e["contig"]].append({"element_id": e["element_id"], "start": s, "end": s + L})
        null.append(cluster_stats(placed))
    perm_rows = []
    for i, name in enumerate(("clusters_overlap_ge1gene", "clusters_overlap_ge50pct", "sm_cluster_genes_in_elements")):
        vals = sorted(x[i] for x in null)
        perm_rows.append({"statistic": name, "observed": obs[i], "n_perm": a.n_perm, "seed": a.seed,
                          "null_mean": round(sum(vals) / len(vals), 3), "null_q95": vals[int(0.95 * (len(vals) - 1))],
                          "null_max": vals[-1],
                          "p_perm_greater": round((sum(1 for v in vals if v >= obs[i]) + 1) / (a.n_perm + 1), 5),
                          "p_perm_less": round((sum(1 for v in vals if v <= obs[i]) + 1) / (a.n_perm + 1), 5)})
    sio.write_tsv(a.out_dir / "sm_permutation.tsv", perm_rows)

    # ---------------- gene-level hypergeometric
    N = len(cds)
    in_el = {n for n, c in cds.items() if c["element"]}
    sm_g = set().union(*(set(r["genes"]) for r in regions)) if regions else set()
    core_g = set().union(*(set(r["core"]) for r in regions)) if regions else set()
    h_all = hyper(N, len(in_el), len(sm_g), len(sm_g & in_el))
    h_core = hyper(N, len(in_el), len(core_g), len(core_g & in_el))

    # ---------------- family view (core genes -> families)
    fam_rows = []
    core_fams = defaultdict(set)
    gene_fams = set()
    for r in regions:
        for x in r["genes"]:
            if cds[x]["family"]:
                gene_fams.add(cds[x]["family"])
        for x in r["core"]:
            if cds[x]["family"]:
                core_fams[cds[x]["family"]].add(r["type"])
    for f in sorted(core_fams):
        fam_rows.append({"family": f, "sm_types": "+".join(sorted(core_fams[f])), "bin": freq.get(f, {}).get("bin", "NA"),
                         "strain_count": freq.get(f, {}).get("strain_count", ""), "direction": gl.get(f, "NA"),
                         "in_island": "yes" if f in island_fams else "no"})
    sio.write_tsv(a.out_dir / "sm_core_families.tsv.gz", fam_rows)
    brk = []
    for label, fs in (("core_gene_families", set(core_fams)), ("all_cluster_gene_families", gene_fams)):
        cnt = Counter((freq.get(f, {}).get("bin", "NA"), gl.get(f, "NA"), f in island_fams) for f in fs)
        for (b, d, i), n in sorted(cnt.items()):
            brk.append({"gene_set": label, "bin": b, "direction": d, "in_island": "yes" if i else "no", "n_families": n})
    sio.write_tsv(a.out_dir / "sm_family_class_direction.tsv", brk)
    cf_nc = set(core_fams) & nonc_f
    h_fam = hyper(len(nonc_f), len(cf_nc), len(island_fams & nonc_f), len(cf_nc & island_fams))
    gf_nc = gene_fams & nonc_f
    h_fam_all = hyper(len(nonc_f), len(gf_nc), len(island_fams & nonc_f), len(gf_nc & island_fams))

    # ---------------- proxy comparison
    with sio.open_text(a.proxy_genes) as fh:
        proxy = list(csv.DictReader(fh, delimiter="\t"))
    mem2lt = {c["member"]: n for n, c in cds.items() if c["member"]}
    proxy_in_region = [p for p in proxy if mem2lt.get(p["protein_id"]) in sm_g]
    proxy_core = [p for p in proxy if mem2lt.get(p["protein_id"]) in core_g]

    # ---------------- summary
    src = "full_v072_antismash"
    S = []

    def M(metric, value, source):
        S.append({"metric": metric, "value": value, "source": source})

    mapc = Counter(c["map"] for c in cds.values())
    M("antismash_version", js.get("version", ""), "Af293.json.gz")
    M("af293_cds_in_antismash", N, "CDS features in antiSMASH JSON (from GFF3)")
    M("af293_cds_mapped_to_member", sum(v for k, v in mapc.items() if k != "none"),
      "; ".join(f"{k}:{v}" for k, v in sorted(mapc.items())))
    M("af293_cds_with_family", sum(1 for c in cds.values() if c["family"]), "member in tier1_cluster.tsv")
    M("sm_clusters", len(regions), "antiSMASH regions (sm_clusters.tsv.gz)")
    M("sm_cluster_genes", len(sm_g), "CDS inside regions")
    M("sm_core_genes", len(core_g), "CDS with a rule definition domain")
    M("sm_clusters_contig_edge", sum(1 for r in regions if r["contig_edge"] == "True"), "region qualifier contig_edge")
    M("starship_elements", len(els), "Table S7 (inputs/af293_starships_table_s7.elements.tsv)")
    M("af293_cds_in_starships", len(in_el), "CDS midpoint inside an element")
    M("sm_clusters_overlap_ge1gene", sum(r["ruleA"] for r in regions), "rule A (sm_clusters.tsv.gz)")
    M("sm_clusters_overlap_ge50pct", sum(r["ruleB"] for r in regions), "rule B (sm_clusters.tsv.gz)")
    M("sm_clusters_with_captain", sum(1 for r in regions if r["captain"]), "DUF3435 hit E<=1e-5 in a cluster gene")
    M("af293_captain_cds", sum(1 for c in cds.values() if c["captain"]), "Af293 CDS whose member has a DUF3435 hit")
    M("elements_with_sm_cluster_ge1gene", sum(1 for p in per_el if p["n_sm_clusters_ge1gene"]), "sm_starship_overlap.tsv")
    M("elements_with_sm_cluster_ge50pct", sum(1 for p in per_el if p["n_sm_clusters_ge50pct"]), "sm_starship_overlap.tsv")
    for pr in perm_rows:
        M(f"perm_{pr['statistic']}_observed", pr["observed"], "sm_permutation.tsv")
        M(f"perm_{pr['statistic']}_null_mean", pr["null_mean"], "sm_permutation.tsv")
        M(f"perm_{pr['statistic']}_p_greater", pr["p_perm_greater"], "sm_permutation.tsv")
        M(f"perm_{pr['statistic']}_p_less", pr["p_perm_less"], "sm_permutation.tsv")
    for nm, h in (("hypergeom_sm_cluster_genes_in_starships", h_all), ("hypergeom_sm_core_genes_in_starships", h_core),
                  ("hypergeom_core_families_in_islands_noncore", h_fam),
                  ("hypergeom_cluster_gene_families_in_islands_noncore", h_fam_all)):
        for k, v in h.items():
            M(f"{nm}_{k}", v, "universe: all Af293 CDS" if "starships" in nm else "universe: non-core families; K = SM families (Af293 genes); n = island families")
    M("sm_core_families", len(core_fams), "families of core genes (sm_core_families.tsv.gz)")
    M("sm_core_families_noncore", len(cf_nc), "")
    M("sm_core_families_noncore_in_island", len(cf_nc & island_fams), "")
    M("proxy_backbone_genes", len(proxy), "full_v072/with_unplaced_af293_coords/sm_backbone_af293_genes.tsv.gz")
    M("proxy_backbone_genes_in_antismash_region", len(proxy_in_region), "proxy gene member -> locus tag in a region")
    M("proxy_backbone_genes_antismash_core", len(proxy_core), "proxy gene is an antiSMASH core gene")
    M("proxy_backbone_genes_in_starships", sum(1 for p in proxy if p["starship_element"]), "proxy result")
    sio.write_tsv(a.out_dir / "summary.tsv", S)

    with open(a.out_dir / "PROVENANCE_analysis.yaml", "w") as fh:
        fh.write(f"derived_by: {Path(__file__).resolve()}\ndate: {dt.date.today().isoformat()}\n"
                 f"command: {' '.join(sys.argv)}\nparameters: {{evalue: {a.evalue}, n_perm: {a.n_perm}, seed: {a.seed}}}\ninputs:\n")
        for k in ("antismash_json", "gff3", "pep", "genome", "extra_gene_positions", "starship_annotation",
                  "captain_tblout", "gain_loss", "proxy_genes"):
            p = Path(getattr(a, k))
            fh.write(f"  {k}:\n    path: {p.resolve()}\n    sha256: {sio.sha256(p)}\n")
        for n in ("frequency_table.tsv", "significant_islands.tsv", "cluster/tier1_cluster.tsv", "gene_positions.tsv.zst"):
            fh.write(f"  run_dir/{n}:\n    path: {(rd / n).resolve()}\n    sha256: {sio.sha256(rd / n)}\n")
    print(f"wrote {a.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
