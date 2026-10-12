#!/usr/bin/env python3
"""Prototype: GO / Pfam / Pfam-clan / InterPro over-representation for the novelty
candidates of one nf_NovInvenio pairwise run.

Two scopes:
  A  whole_set   all candidates of all IN species vs all IN proteins
  B  <Short>     one IN species' candidates vs that species' proteins

Two counting levels:
  protein  one unit per protein
  family   candidate proteins collapsed to one unit per mmseqs family
           (novelties.<Short>.tsv family_id; empty family_id = singleton).
           Scope A: one unit per family across all IN species.
           Scope B: one unit per family within that species.
           The unit's terms are the union of its annotated members' terms.
           Non-candidate background proteins stay one unit per protein, because
           the pipeline does not cluster the background. This makes the family
           level conservative on the candidate side only.

Inputs (all read-only, from the run's results dir):
  presence_matrix.uniprot.tsv  every IN protein the pipeline tested, with the
                               merged uniprot_go_ids / uniprot_pfam_ids /
                               uniprot_pfam_names / uniprot_interpro_ids columns
  candidates.txt               "<Short>::<protein_id>" per line
  novelties.<Short>.tsv        family_id per candidate

Statistics: one-sided hypergeometric P(X >= k). Universe per test = units with >=1
term in that vocabulary (unannotated units are dropped, as in bin/domain_enrichment.py).
Terms with fewer than --min-bg background units are not tested. BH-FDR within each
test (scope x level x vocabulary) -> fdr_bh. For scope B, a second BH across all
species tests of the same level x vocabulary -> fdr_bh_all_species.

GO: true-path propagation over is_a edges of go-basic.obo (goatools GODag); the
three namespace roots are dropped. IEA evidence is included by default.
Pfam clans: from Pfam-A.hmm.dat "#=GF CL" lines (optional --pfam-dat).
InterPro names/types: optional --interpro-entries (InterProScan entries.json).
Names are labels only. They do not change any statistic.

Outputs (--outdir): coverage.tsv, tests.tsv (one row per test, with the number of
tested terms), enrichment.tsv.gz (every tested term).
"""
import argparse
import csv
import glob
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import hypergeom

# Hardcoded per CLAUDE.md (study-specific script location; no __file__-relative paths).
NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from stats_util import bh_fdr  # noqa: E402

csv.field_size_limit(10**9)
GO_ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
VOCABS = ["go", "pfam", "pfam_clan", "interpro"]


def split_ids(s):
    return [x for x in s.split("|") if x] if s else []


def parse_go(s, exclude_iea):
    out = set()
    for entry in split_ids(s):
        go_id, _, ev = entry.rpartition(":")
        if not go_id.startswith("GO:"):
            go_id, ev = entry, ""
        if exclude_iea and ev == "IEA":
            continue
        out.add(go_id)
    return out


def load_pfam_clans(path):
    clan, acc = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                acc = line.split()[2].split(".")[0]
            elif line.startswith("#=GF CL") and acc:
                clan[acc] = line.split()[2]
            elif line.startswith("//"):
                acc = None
    return clan


def load_obo(path):
    from goatools.obo_parser import GODag
    dag = GODag(str(path), prt=None)
    return dag


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--obo", type=Path, default=NII_ROOT / "data/go/go-basic.obo")
    ap.add_argument("--pfam-dat", type=Path, default=None, help="Pfam-A.hmm.dat for clan mapping (optional)")
    ap.add_argument("--interpro-entries", type=Path, default=None, help="InterProScan entries.json for IPR names (optional)")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--min-bg", type=int, default=5, help="Minimum background units carrying a term to test it (default 5)")
    ap.add_argument("--exclude-iea", action="store_true")
    ap.add_argument("--coverage-only", action="store_true", help="Write coverage.tsv and stop")
    ap.add_argument("--report-ids-fdr", type=float, default=0.05, help="List candidate IDs for terms with fdr_bh below this")
    args = ap.parse_args()
    rd = args.results_dir
    args.outdir.mkdir(parents=True, exist_ok=True)

    candidates = {l.strip() for l in open(rd / "candidates.txt") if l.strip()}

    # family_id per candidate
    family = {}
    for f in glob.glob(str(rd / "novelties.*.tsv")):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = f"{row['source_proteome']}::{row['protein_id']}"
                family[key] = row.get("family_id") or f"singleton:{key}"

    # GO DAG
    dag = None
    if not args.coverage_only:
        print(f"loading {args.obo}", file=sys.stderr)
        dag = load_obo(args.obo)
    anc_cache = {}

    def propagate(terms):
        out = set()
        for t in terms:
            if t not in anc_cache:
                rec = dag.get(t) if dag is not None else None
                if rec is None:
                    anc_cache[t] = None  # obsolete or missing in this obo release
                else:
                    anc_cache[t] = {rec.id} | rec.get_all_parents()
            a = anc_cache[t]
            if a:
                out |= a
        return out - GO_ROOTS

    clans = load_pfam_clans(args.pfam_dat) if args.pfam_dat else {}
    names = defaultdict(dict)  # vocab -> id -> name
    pfam_c = args.pfam_dat.parent / "Pfam-C" if args.pfam_dat else None
    if pfam_c and pfam_c.exists():  # clan labels: "#=GF ID" precedes "#=GF AC" in Pfam-C
        cid = None
        with open(pfam_c, errors="replace") as fh:
            for line in fh:
                if line.startswith("#=GF ID"):
                    cid = line.split(maxsplit=2)[2].strip()
                elif line.startswith("#=GF AC") and cid:
                    names["pfam_clan"][line.split()[2].split(".")[0]] = cid
                    cid = None

    # per protein annotations
    prot = {}  # key -> dict(species, cand, vocab->set)
    with open(rd / "presence_matrix.uniprot.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sp = row["source_proteome"]
            key = f"{sp}::{row['protein_id']}"
            pf = split_ids(row["uniprot_pfam_ids"])
            pfn = split_ids(row["uniprot_pfam_names"])
            if len(pf) == len(pfn):
                for a, n in zip(pf, pfn):
                    names["pfam"][a] = n
            go_raw = parse_go(row["uniprot_go_ids"], args.exclude_iea)
            ann = {
                "go_raw": go_raw,
                "pfam": set(pf),
                "pfam_clan": {clans[p] for p in pf if p in clans},
                "interpro": set(split_ids(row["uniprot_interpro_ids"])),
            }
            prot[key] = {"sp": sp, "cand": key in candidates, "ann": ann}

    missing = candidates - prot.keys()
    if missing:
        sys.exit(f"ERROR: {len(missing)} candidates not in presence_matrix.uniprot.tsv, e.g. {sorted(missing)[:3]}")
    species = sorted({p["sp"] for p in prot.values()})

    # coverage table
    cov_rows = []
    for sp in species + ["ALL_IN"]:
        sel = [p for p in prot.values() if sp == "ALL_IN" or p["sp"] == sp]
        for subset, rows in (("all", sel), ("candidates", [p for p in sel if p["cand"]])):
            n = len(rows)
            r = {"species": sp, "subset": subset, "n_proteins": n}
            for v, f in (("go", "go_raw"), ("pfam", "pfam"), ("pfam_clan", "pfam_clan"), ("interpro", "interpro")):
                k = sum(1 for p in rows if p["ann"][f])
                r[f"n_{v}"] = k
                r[f"frac_{v}"] = round(k / n, 4) if n else 0
            k = sum(1 for p in rows if p["ann"]["go_raw"] or p["ann"]["pfam"] or p["ann"]["interpro"])
            r["n_any"] = k
            r["frac_any"] = round(k / n, 4) if n else 0
            cov_rows.append(r)
    with open(args.outdir / "coverage.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cov_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(cov_rows)
    print(f"wrote {args.outdir/'coverage.tsv'}", file=sys.stderr)
    if args.coverage_only:
        return 0

    for p in prot.values():
        p["ann"]["go"] = propagate(p["ann"]["go_raw"])

    if dag is not None:
        for t, rec in dag.items():
            names["go"][t] = f"{rec.namespace}:{rec.name}"
    if args.interpro_entries:
        ent = json.load(open(args.interpro_entries))
        for k, v in ent.items():
            if k.startswith("IPR"):
                names["interpro"][k] = f"{v.get('type','')}:{v.get('description') or v.get('name','')}"

    def build_units(sp_filter, level):
        """Return list of (unit_id, is_candidate, member_keys, terms_by_vocab)."""
        units = []
        fam_units = defaultdict(list)
        for key, p in prot.items():
            if sp_filter and p["sp"] != sp_filter:
                continue
            if p["cand"] and level == "family":
                fam_units[family.get(key, f"singleton:{key}")].append(key)
            else:
                units.append((key, p["cand"], [key]))
        for fid, members in fam_units.items():
            units.append((f"family:{fid}", True, members))
        out = []
        for uid, is_c, members in units:
            terms = {v: set().union(*(prot[m]["ann"][v] for m in members)) for v in VOCABS}
            out.append((uid, is_c, members, terms))
        return out

    scopes = [("whole_set", None)] + [(sp, sp) for sp in species]
    all_rows, test_rows = [], []
    for scope, spf in scopes:
        for level in ("protein", "family"):
            units = build_units(spf, level)
            for v in VOCABS:
                if v == "pfam_clan" and not clans:
                    continue
                uni = [u for u in units if u[3][v]]
                M = len(uni)
                cands = [u for u in uni if u[1]]
                N = len(cands)
                n_c_total = sum(1 for u in units if u[1])
                bg = defaultdict(int)
                kc = defaultdict(int)
                cid = defaultdict(list)
                for u in uni:
                    for t in u[3][v]:
                        bg[t] += 1
                        if u[1]:
                            kc[t] += 1
                            cid[t].extend(u[2])
                rows = []
                if N > 0:
                    for t, n in bg.items():
                        if n < args.min_bg:
                            continue
                        k = kc.get(t, 0)
                        p = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
                        exp = n * N / M
                        rows.append({
                            "scope": scope, "level": level, "vocab": v, "term": t,
                            "name": names[v].get(t, ""),
                            "k_cand_with_term": k, "N_cand_annotated": N,
                            "n_bg_with_term": n, "M_bg_annotated": M,
                            "expected": round(exp, 3),
                            "fold": round(k / exp, 3) if exp else "",
                            "pvalue": p,
                        })
                    for r, q in zip(rows, bh_fdr([r["pvalue"] for r in rows])):
                        r["fdr_bh"] = q
                        r["candidate_ids"] = ";".join(sorted(cid[r["term"]])) if q < args.report_ids_fdr else ""
                n_sig = sum(1 for r in rows if r["fdr_bh"] < 0.05)
                test_rows.append({
                    "scope": scope, "level": level, "vocab": v,
                    "n_cand_units": n_c_total, "N_cand_annotated": N,
                    "frac_cand_annotated": round(N / n_c_total, 4) if n_c_total else 0,
                    "M_bg_annotated": M, "n_bg_units": len(units),
                    "n_terms_tested": len(rows), "n_fdr_lt_0.05": n_sig,
                    "status": "ok" if N > 0 else "skipped_no_annotated_candidates",
                })
                all_rows.extend(rows)
                print(f"{scope}\t{level}\t{v}\tN={N}/{n_c_total}\tM={M}\ttested={len(rows)}\tsig={n_sig}", file=sys.stderr)

    # BH across all species-level (scope B) tests of one level x vocab
    grp = defaultdict(list)
    for i, r in enumerate(all_rows):
        if r["scope"] != "whole_set":
            grp[(r["level"], r["vocab"])].append(i)
    for idx in grp.values():
        for i, q in zip(idx, bh_fdr([all_rows[i]["pvalue"] for i in idx])):
            all_rows[i]["fdr_bh_all_species"] = q
    for r in all_rows:
        r.setdefault("fdr_bh_all_species", "")

    all_rows.sort(key=lambda r: (r["scope"] != "whole_set", r["scope"], r["level"], r["vocab"], r["pvalue"]))
    fields = ["scope", "level", "vocab", "term", "name", "k_cand_with_term", "N_cand_annotated",
              "n_bg_with_term", "M_bg_annotated", "expected", "fold", "pvalue", "fdr_bh",
              "fdr_bh_all_species", "candidate_ids"]
    with gzip.open(args.outdir / "enrichment.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)
    with open(args.outdir / "tests.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(test_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(test_rows)
    print(f"wrote {args.outdir}/enrichment.tsv.gz ({len(all_rows)} rows) and tests.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
