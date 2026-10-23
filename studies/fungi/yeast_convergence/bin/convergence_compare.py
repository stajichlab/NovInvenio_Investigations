#!/usr/bin/env python3
"""Cross-contrast convergence test for the yeast_convergence studies.

Reads every contrast's enrichment tables (enrich_novelties_v2.py output, written
by run_enrichment_all.sh) and asks which functional terms recur across the
independent yeast origins more often than chance, and whether the same terms
also recur in the filamentous pseudo-origin controls.

Inputs, per contrast <id> in study_index.tsv and direction in {novelty, loss}:
    <folder>/<id>/analysis/enrichment/outputs/<direction>/proteome_clustered/
        enrichment_domains.tsv.gz   (pfam, pfam_clan, interpro)
        enrichment_go.tsv.gz        (go)
A contrast or direction with no table is reported as missing and left out.

A term is "significant" in a contrast when, at the chosen scope (default
whole_set) and level (default family): fdr_bh <= --fdr, k_cand_with_term >=
--min-k and fold > 1. A term is "tested" in a contrast when it has a row there.

Recurrence test (per direction x vocabulary, over the ORIGIN contrasts):
    r_obs(t)  = number of origin contrasts where t is significant.
    Null: contrast i calls a tested term significant with probability
    p_i = |significant_i| / |tested_i|, independently of other contrasts.
    So r(t) ~ Poisson-binomial over the origin contrasts that tested t, and
    p_recur(t) = P(r >= r_obs). BH-adjusted over every term tested in >= 2
    origin contrasts.
    This null asks "more recurrence than independent draws of the same size
    would give"; it does not model terms that are enriched everywhere because
    of family size. The pseudo-origin column is the check for that.

Classification of terms with r_obs >= --min-origins and q_recur <= --alpha:
    yeast_convergent      significant in no pseudo-origin control
    streamlining_generic  also significant in >= 1 pseudo-origin control
n_pseudo_deep_sig counts the controls whose split_rank is class or deeper
(the depth of most yeast origins). A yeast_convergent term that the deep
controls tested but did not call is the strongest result; a
streamlining_generic term called only by shallow controls is worth a look,
since shallow splits differ from yeast origins in depth as well as biology.

Small ingroups: the novelty tables of contrasts with fewer than
--min-in-novelty ingroup species (default 2) are not read; their loss tables are.

Annotation tiers:
    primary    pfam, pfam_clan: every contrast.
    secondary  go, interpro: only contrasts whose species.csv lists every
               species with Protein_Source=uniprot (all others are excluded
               from that vocabulary and listed in annotation_coverage.tsv).
               --include-mixed-annotation turns the exclusion off.
    Classes are assigned from primary vocabularies; secondary rows carry the
    same columns but are labelled tier=secondary and are a cross-check.

Outputs (in --out-dir):
    term_recurrence.tsv   one row per (direction, vocab, term) tested in >= 1 origin
    pairwise_jaccard.tsv  Jaccard of significant sets for every contrast pair,
                          over terms tested in both, with the role pair
    reversal_overlap.tsv  terms LOST in an origin and gained (novelty) in its
                          reversal, and the opposite
    focus_terms.tsv       status of focus_terms.tsv families in every contrast
    annotation_coverage.tsv  per contrast: species, UniProt species, whether
                          GO/InterPro were used
    summary.md            counts, coverage, top terms
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import sys
from collections import defaultdict
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER = NII_ROOT / "studies/fungi/yeast_convergence"
DIRECTIONS = ("novelty", "loss")
TABLE_FOR_VOCAB = {"pfam": "enrichment_domains.tsv.gz", "pfam_clan": "enrichment_domains.tsv.gz",
                   "interpro": "enrichment_domains.tsv.gz", "go": "enrichment_go.tsv.gz"}
# Every proteome gets Pfam (UniProt records, or the HMMER step for NCBI/JGI
# proteomes), so Pfam and Pfam clan are comparable across all contrasts. GO and
# InterPro exist only for UniProt proteomes (CLAUDE.md "Annotating non-UniProt
# proteomes"), so a contrast with any non-UniProt species under-calls them.
PRIMARY_VOCABS = ("pfam", "pfam_clan")
# split_rank values (study_index.tsv) counted as deep controls.
DEEP_RANKS = {"class", "subphylum", "phylum"}
SECONDARY_VOCABS = ("go", "interpro")


def read_index(folder: Path) -> list[dict]:
    with open(folder / "study_index.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_focus(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = [line for line in fh if not line.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def load_calls(table: Path, vocab: str, scope: str, level: str, fdr: float, min_k: int):
    """Return (tested: {term: (name, k, fold, fdr)}, significant: set(term))."""
    tested, sig = {}, set()
    with gzip.open(table, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["vocab"] != vocab or row["scope"] != scope or row["level"] != level:
                continue
            try:
                k = int(row["k_cand_with_term"])
                fold = float(row["fold"])
                q = float(row["fdr_bh"]) if row["fdr_bh"] not in ("", "NA") else 1.0
            except ValueError:
                continue
            tested[row["term"]] = (row["name"], k, fold, q)
            if q <= fdr and k >= min_k and fold > 1:
                sig.add(row["term"])
    return tested, sig


def poisson_binomial_sf(ps: list[float], r: int) -> float:
    """P(X >= r) for X = sum of independent Bernoulli(p_i)."""
    if r <= 0:
        return 1.0
    dist = [1.0]
    for p in ps:
        nxt = [0.0] * (len(dist) + 1)
        for j, v in enumerate(dist):
            nxt[j] += v * (1 - p)
            nxt[j + 1] += v * p
        dist = nxt
    return min(1.0, sum(dist[r:]))


def bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [1.0] * m
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


def annotation_coverage(folder: Path, index: list[dict]) -> dict[str, dict]:
    """Per contrast: number of species and of UniProt-sourced species, from
    <id>/species.csv. A blank Protein_Source counts as not UniProt."""
    cov = {}
    for st in index:
        path = folder / st["id"] / "species.csv"
        n = n_up = 0
        if path.exists():
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh):
                    n += 1
                    n_up += row.get("Protein_Source") == "uniprot"
        cov[st["id"]] = {"id": st["id"], "role": st["role"], "n_species": n, "n_uniprot": n_up,
                         "all_uniprot": n > 0 and n == n_up}
    return cov


def gather(folder: Path, index: list[dict], vocabs, scope, level, fdr, min_k,
           coverage: dict[str, dict] | None = None, include_mixed: bool = False,
           min_in_novelty: int = 1):
    """calls[(id, direction, vocab)] = (tested, sig); missing = list of (id, direction).
    Secondary vocabularies are skipped for contrasts that are not all-UniProt
    unless include_mixed. The novelty direction is skipped for contrasts with
    fewer than min_in_novelty ingroup species: a one-genome "novelty" list is
    mostly that genome's own orphans and annotation artefacts."""
    calls, missing = {}, []
    for st in index:
        all_up = include_mixed or coverage is None or coverage[st["id"]]["all_uniprot"]
        for d in DIRECTIONS:
            if d == "novelty" and int(st.get("n_in") or 0) < min_in_novelty:
                continue
            base = folder / st["id"] / "analysis/enrichment/outputs" / d / "proteome_clustered"
            got_any = False
            for v in vocabs:
                if v in SECONDARY_VOCABS and not all_up:
                    continue
                table = base / TABLE_FOR_VOCAB[v]
                if table.exists():
                    calls[(st["id"], d, v)] = load_calls(table, v, scope, level, fdr, min_k)
                    got_any = True
            if not got_any:
                missing.append((st["id"], d))
    return calls, missing


def recurrence(index, calls, vocabs, alpha, min_origins):
    origins = [st["id"] for st in index if st["role"] == "origin"]
    pseudos = [st["id"] for st in index if st["role"] == "pseudo"]
    deep = {st["id"] for st in index if st["role"] == "pseudo" and st.get("split_rank") in DEEP_RANKS}
    reversals = [st["id"] for st in index if st["role"] == "reversal"]
    out = []
    for d, v in itertools.product(DIRECTIONS, vocabs):
        have = [o for o in origins if (o, d, v) in calls]
        if not have:
            continue
        rate = {o: (len(calls[(o, d, v)][1]) / len(calls[(o, d, v)][0])) if calls[(o, d, v)][0] else 0.0
                for o in have}
        terms = sorted(set().union(*(calls[(o, d, v)][0] for o in have)))
        rows = []
        for t in terms:
            tested_in = [o for o in have if t in calls[(o, d, v)][0]]
            sig_in = [o for o in tested_in if t in calls[(o, d, v)][1]]
            name = next(calls[(o, d, v)][0][t][0] for o in tested_in)
            p = poisson_binomial_sf([rate[o] for o in tested_in], len(sig_in))
            ps_tested = [s for s in pseudos if (s, d, v) in calls and t in calls[(s, d, v)][0]]
            ps_sig = [s for s in ps_tested if t in calls[(s, d, v)][1]]
            rv_sig = [s for s in reversals if (s, d, v) in calls and t in calls[(s, d, v)][1]]
            rows.append({
                "tier": "primary" if v in PRIMARY_VOCABS else "secondary",
                "direction": d, "vocab": v, "term": t, "name": name,
                "n_origin_tested": len(tested_in), "n_origin_sig": len(sig_in),
                "origins_sig": ";".join(sig_in), "p_recur": p, "q_recur": None,
                "n_pseudo_tested": len(ps_tested), "n_pseudo_sig": len(ps_sig),
                "n_pseudo_deep_tested": len([s for s in ps_tested if s in deep]),
                "n_pseudo_deep_sig": len([s for s in ps_sig if s in deep]),
                "pseudo_sig": ";".join(ps_sig), "reversal_sig": ";".join(rv_sig),
                "class": "",
            })
        testable = [i for i, r in enumerate(rows) if r["n_origin_tested"] >= 2]
        qs = bh([rows[i]["p_recur"] for i in testable])
        for i, q in zip(testable, qs):
            rows[i]["q_recur"] = q
        for r in rows:
            if r["q_recur"] is None:
                r["q_recur"] = 1.0
            if r["n_origin_sig"] >= min_origins:
                if r["q_recur"] <= alpha:
                    r["class"] = "yeast_convergent" if r["n_pseudo_sig"] == 0 else "streamlining_generic"
                else:
                    r["class"] = "recurrent_not_significant"
        out.extend(rows)
    return out


def pairwise_jaccard(index, calls, vocabs):
    role = {st["id"]: st["role"] for st in index}
    ids = [st["id"] for st in index]
    out = []
    for d, v in itertools.product(DIRECTIONS, vocabs):
        for a, b in itertools.combinations(ids, 2):
            if (a, d, v) not in calls or (b, d, v) not in calls:
                continue
            ta, sa = calls[(a, d, v)]
            tb, sb = calls[(b, d, v)]
            shared = set(ta) & set(tb)
            sa, sb = sa & shared, sb & shared
            union = sa | sb
            out.append({
                "direction": d, "vocab": v, "a": a, "b": b,
                "role_pair": "-".join(sorted((role[a], role[b]))),
                "n_shared_tested": len(shared), "n_sig_a": len(sa), "n_sig_b": len(sb),
                "n_both": len(sa & sb), "jaccard": (len(sa & sb) / len(union)) if union else 0.0,
            })
    return out


def reversal_overlap(index, calls, vocabs):
    out = []
    for st in index:
        if st["role"] != "reversal" or not st.get("compare_with"):
            continue
        rev, ori = st["id"], st["compare_with"]
        for v in vocabs:
            for ori_dir, rev_dir, label in (("loss", "novelty", "lost_in_yeast_origin__gained_in_reversal"),
                                            ("novelty", "loss", "gained_in_yeast_origin__lost_in_reversal")):
                if (ori, ori_dir, v) not in calls or (rev, rev_dir, v) not in calls:
                    continue
                t_o, s_o = calls[(ori, ori_dir, v)]
                t_r, s_r = calls[(rev, rev_dir, v)]
                for t in sorted(s_o & s_r):
                    out.append({"reversal": rev, "origin": ori, "vocab": v, "pattern": label,
                                "term": t, "name": t_o[t][0],
                                "origin_k": t_o[t][1], "origin_fdr": t_o[t][3],
                                "reversal_k": t_r[t][1], "reversal_fdr": t_r[t][3]})
    return out


def focus_status(index, calls, focus):
    out = []
    for f in focus:
        for st in index:
            for d in DIRECTIONS:
                key = (st["id"], d, "pfam")
                if key not in calls:
                    status, k, fold, q = "no_table", "", "", ""
                else:
                    tested, sig = calls[key]
                    if f["term"] not in tested:
                        status, k, fold, q = "not_tested", "", "", ""
                    else:
                        _, k, fold, q = tested[f["term"]]
                        status = "significant" if f["term"] in sig else "tested_ns"
                out.append({"term": f["term"], "name": f["name"],
                            "expected_in_yeast": f["expected_in_yeast"], "id": st["id"],
                            "role": st["role"], "direction": d, "status": status,
                            "k": k, "fold": fold, "fdr": q})
    return out


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{x:.4g}" if isinstance(x, float) else x) for k, x in r.items()})


def summary_md(args, index, missing, rec, jac, rev, focus_rows, coverage=None) -> str:
    lines = ["# yeast_convergence: cross-contrast comparison", "",
             f"Scope `{args.scope}`, level `{args.level}`, significant = FDR <= {args.fdr}, "
             f"k >= {args.min_k}, fold > 1. Recurrence: Poisson-binomial null, BH q <= {args.alpha}, "
             f"in >= {args.min_origins} origins.", ""]
    n_by_role = defaultdict(int)
    for st in index:
        n_by_role[st["role"]] += 1
    lines.append("Contrasts: " + ", ".join(f"{n} {r}" for r, n in sorted(n_by_role.items())))
    single = [st["id"] for st in index if int(st.get("n_in") or 0) < args.min_in_novelty]
    if single:
        lines += ["", f"Novelty lists not used (fewer than {args.min_in_novelty} ingroup species; "
                  f"losses still used): {', '.join(single)}."]
    if coverage is not None:
        mixed = [c["id"] for c in coverage.values() if not c["all_uniprot"]]
        if args.include_mixed_annotation:
            lines += ["", "GO/InterPro used for every contrast (--include-mixed-annotation); "
                      f"{len(mixed)} contrasts have non-UniProt species, so read those vocabularies with care."]
        else:
            lines += ["", f"Primary call: {', '.join(PRIMARY_VOCABS)} over all contrasts. "
                      f"GO/InterPro (secondary) only over the {len(coverage) - len(mixed)} all-UniProt contrasts; "
                      f"excluded: {', '.join(mixed) if mixed else 'none'}."]
    if missing:
        lines += ["", "Missing tables (left out): " + ", ".join(f"{i}/{d}" for i, d in missing)]
    for d in DIRECTIONS:
        lines += ["", f"## {d}", ""]
        for cls in ("yeast_convergent", "streamlining_generic"):
            hits = sorted((r for r in rec if r["direction"] == d and r["class"] == cls
                           and r["tier"] == "primary"),
                          key=lambda r: (r["q_recur"], -r["n_origin_sig"]))
            n_sec = sum(1 for r in rec if r["direction"] == d and r["class"] == cls and r["tier"] == "secondary")
            lines.append(f"**{cls}** (primary): {len(hits)} terms; secondary GO/InterPro cross-check: {n_sec}")
            for r in hits[:15]:
                lines.append(f"- {r['vocab']} {r['term']} {r['name']}: {r['n_origin_sig']}/"
                             f"{r['n_origin_tested']} origins, q={r['q_recur']:.3g}, "
                             f"controls {r['n_pseudo_sig']}/{r['n_pseudo_tested']} "
                             f"(deep {r['n_pseudo_deep_sig']}/{r['n_pseudo_deep_tested']})")
            lines.append("")
        means = defaultdict(list)
        for j in jac:
            if j["direction"] == d and j["vocab"] == "pfam":
                means[j["role_pair"]].append(j["jaccard"])
        if means:
            lines.append("Mean Pfam Jaccard of significant sets by role pair: " + ", ".join(
                f"{k} {sum(v) / len(v):.3f} (n={len(v)})" for k, v in sorted(means.items())))
    if rev:
        lines += ["", "## Reversals", ""]
        for r in rev:
            lines.append(f"- {r['reversal']} vs {r['origin']} ({r['pattern']}): {r['vocab']} {r['term']} {r['name']}")
    sig_focus = [f for f in focus_rows if f["status"] == "significant"]
    lines += ["", "## Focus families (Nagy et al. 2014)", "",
              f"{len(sig_focus)} significant calls across contrasts; full grid in focus_terms.tsv."]
    for f in sig_focus:
        lines.append(f"- {f['term']} {f['name']} ({f['expected_in_yeast']} expected): "
                     f"{f['id']} [{f['role']}] {f['direction']}, k={f['k']}, fdr={f['fdr']:.3g}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folder", type=Path, default=FOLDER)
    ap.add_argument("--out-dir", type=Path, default=None, help="default: <folder>/analysis/convergence")
    ap.add_argument("--vocab", nargs="+", default=["pfam", "pfam_clan", "interpro", "go"],
                    choices=list(TABLE_FOR_VOCAB))
    ap.add_argument("--scope", default="whole_set")
    ap.add_argument("--level", default="family", choices=["family", "protein"])
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-k", type=int, default=3)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-origins", type=int, default=3)
    ap.add_argument("--min-in-novelty", type=int, default=2,
                    help="ignore the novelty tables of contrasts with fewer ingroup species (default 2)")
    ap.add_argument("--include-mixed-annotation", action="store_true",
                    help="use GO/InterPro tables even for contrasts with non-UniProt species")
    args = ap.parse_args(argv)

    out_dir = args.out_dir or args.folder / "analysis/convergence"
    out_dir.mkdir(parents=True, exist_ok=True)
    index = read_index(args.folder)
    focus = read_focus(args.folder / "focus_terms.tsv")
    coverage = annotation_coverage(args.folder, index)
    calls, missing = gather(args.folder, index, args.vocab, args.scope, args.level, args.fdr, args.min_k,
                            coverage=coverage, include_mixed=args.include_mixed_annotation,
                            min_in_novelty=args.min_in_novelty)
    if not calls:
        print("ERROR: no enrichment tables found; run run_enrichment_all.sh first", file=sys.stderr)
        return 1

    rec = recurrence(index, calls, args.vocab, args.alpha, args.min_origins)
    jac = pairwise_jaccard(index, calls, args.vocab)
    rev = reversal_overlap(index, calls, args.vocab)
    foc = focus_status(index, calls, focus)

    write_tsv(out_dir / "term_recurrence.tsv", rec,
              ["tier", "direction", "vocab", "term", "name", "n_origin_tested", "n_origin_sig", "origins_sig",
               "p_recur", "q_recur", "n_pseudo_tested", "n_pseudo_sig", "n_pseudo_deep_tested",
               "n_pseudo_deep_sig", "pseudo_sig", "reversal_sig", "class"])
    write_tsv(out_dir / "pairwise_jaccard.tsv", jac,
              ["direction", "vocab", "a", "b", "role_pair", "n_shared_tested", "n_sig_a", "n_sig_b",
               "n_both", "jaccard"])
    write_tsv(out_dir / "reversal_overlap.tsv", rev,
              ["reversal", "origin", "vocab", "pattern", "term", "name", "origin_k", "origin_fdr",
               "reversal_k", "reversal_fdr"])
    write_tsv(out_dir / "focus_terms.tsv", foc,
              ["term", "name", "expected_in_yeast", "id", "role", "direction", "status", "k", "fold", "fdr"])
    write_tsv(out_dir / "annotation_coverage.tsv",
              [dict(c, secondary_vocabs_used=args.include_mixed_annotation or c["all_uniprot"])
               for c in coverage.values()],
              ["id", "role", "n_species", "n_uniprot", "all_uniprot", "secondary_vocabs_used"])
    (out_dir / "summary.md").write_text(summary_md(args, index, missing, rec, jac, rev, foc, coverage))
    print(f"wrote {out_dir}/: term_recurrence.tsv pairwise_jaccard.tsv reversal_overlap.tsv "
          f"focus_terms.tsv annotation_coverage.tsv summary.md ({len(missing)} missing tables)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
