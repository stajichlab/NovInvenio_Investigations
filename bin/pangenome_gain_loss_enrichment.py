#!/usr/bin/env python3
"""Over-representation of Pfam / Pfam-clan terms among GAINED and LOST pangenome families.

Input 1, --gain-loss-table: family_gain_loss.tsv.gz of bin/pangenome_gain_loss.py (one row
per family: family, class, rule_freq, direction, outgroup_count, ...).
Input 2, --family-annotation: one row per family in the column layout of
bin/hmm_domtbl_to_annotations.py / bin/extract_dat_annotations.py (key column `accession`
or `family`; `pfam_ids` "|"-separated; optional `go_ids` "GO:nnnnnnn[:EVIDENCE]|...",
`interpro_ids`). For a pangenome run the key is the mmseqs cluster representative id,
which is also the family id; a family "has" a term when its representative has it.

Scope: the families whose `class` is in --classes (default shell,soft_core).
Test sets: families in scope whose --direction-column (default rule_freq) equals each
value of --test-sets (default gain,loss). Each test set is tested against the scope.

Statistics (same as studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py):
  per vocabulary v, universe M = families in scope with >= 1 term of v;
  N = test-set families in that universe; K = universe families with the term;
  k = test-set families with the term.
  p = P(X >= k), X ~ Hypergeom(M, K, N) = scipy.stats.hypergeom.sf(k - 1, M, K, N)
  (p = 1 when k = 0). A term is tested only if K >= --min-universe (default 5).
  BH FDR (lib/stats_util.bh_fdr) within one (test set, vocabulary).
  Vocabularies: pfam; pfam_clan (Pfam -> clan from Pfam-A.hmm.dat, --pfam-dat);
  go and interpro only when the annotation table has values for them (GO terms are used
  as given, without propagation to parent terms).

Outputs (--out-dir):
  enrichment_domains.tsv.gz  pfam + pfam_clan tests (columns below), sorted by
                             test_set, vocab, p, term
  enrichment_go.tsv.gz       go + interpro tests (only when present)
  tests.tsv                  per (test set, vocab): M, N, terms tested, terms FDR < 0.05
                             (all, and with k >= 3)
  coverage.tsv               per direction value in scope: families, with >= 1 term per vocab
  PROVENANCE.yaml
Columns of the enrichment tables: label, scope, test_set, direction (= test_set), vocab,
term, name, k, n_set (= N), K, M, expected (= K*N/M), fold (= k/expected), p, q,
n_families_with_term_in_set (= k), family_ids (test-set families with the term, only
when q < --report-ids-fdr).

Usage:
  bin/pangenome_gain_loss_enrichment.py \\
      --gain-loss-table studies/fungi/<s>/analysis/gain_loss/<run>/family_gain_loss.tsv.gz \\
      --family-annotation studies/fungi/<s>/analysis/gain_loss_enrichment/<run>/family_pfam.tsv.gz \\
      --out-dir studies/fungi/<s>/analysis/gain_loss_enrichment/<run>/shell_soft_core \\
      --classes shell,soft_core --label <run>
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import platform
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))

DEFAULT_PFAM_DAT = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat"
DOMAIN_VOCABS = ("pfam", "pfam_clan")
OTHER_VOCABS = ("go", "interpro")
COLUMNS = ["label", "scope", "test_set", "direction", "vocab", "term", "name", "k", "n_set",
           "K", "M", "expected", "fold", "p", "q", "n_families_with_term_in_set", "family_ids"]


# ---------------------------------------------------------------- I/O
def open_text(path: Path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    if p.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dc", p], stdout=subprocess.PIPE)
        return io.TextIOWrapper(proc.stdout)
    return open(p)


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    text = "\t".join(header) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows)
    if str(path).endswith(".gz"):
        with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            gz.write(text.encode())
    else:
        path.write_text(text)


def read_gain_loss(path: Path, direction_column: str) -> list[dict]:
    with open_text(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if rows and (direction_column not in rows[0] or "class" not in rows[0] or "family" not in rows[0]):
        raise SystemExit(f"ERROR: {path} needs columns family, class, {direction_column}")
    return rows


def split_terms(value: str, vocab: str) -> set[str]:
    out = set()
    for t in (value or "").split("|"):
        t = t.strip()
        if not t:
            continue
        if vocab == "pfam":
            t = t.split(".")[0]
        elif vocab == "go":
            t = ":".join(t.split(":")[:2])   # GO:0030729:IBA -> GO:0030729
        out.add(t)
    return out


def read_annotation(path: Path) -> tuple[dict[str, dict[str, set[str]]], dict[str, str]]:
    """family -> {vocab: terms}; pfam accession -> short name (from pfam_names)."""
    ann: dict[str, dict[str, set[str]]] = {}
    pfam_short: dict[str, str] = {}
    with open_text(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        key = "family" if "family" in rd.fieldnames else "accession"
        if key not in rd.fieldnames or "pfam_ids" not in rd.fieldnames:
            raise SystemExit(f"ERROR: {path} needs a family/accession column and pfam_ids")
        for r in rd:
            terms = {"pfam": split_terms(r.get("pfam_ids", ""), "pfam"),
                     "go": split_terms(r.get("go_ids", ""), "go"),
                     "interpro": split_terms(r.get("interpro_ids", ""), "interpro")}
            fam = r[key]
            if fam in ann:
                raise SystemExit(f"ERROR: family {fam} appears twice in {path}")
            ann[fam] = terms
            ids = [x.split(".")[0] for x in (r.get("pfam_ids") or "").split("|") if x]
            names = (r.get("pfam_names") or "").split("|")
            if len(ids) == len(names):
                for a, n in zip(ids, names):
                    pfam_short.setdefault(a, n)
    return ann, pfam_short


def load_pfam_dat(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Return (Pfam accession -> clan accession, Pfam accession -> 'ID: DE').
    Same parser as enrich_novelties_v2.py load_pfam_dat()."""
    clan, names, acc, pid, de = {}, {}, None, None, ""
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                pid = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF AC"):
                acc = line.split()[2].split(".")[0]
            elif line.startswith("#=GF DE"):
                de = line[10:].strip()
            elif line.startswith("#=GF CL") and acc:
                clan[acc] = line.split()[2]
            elif line.startswith("//"):
                if acc:
                    names[acc] = f"{pid}: {de}" if de else (pid or "")
                acc, pid, de = None, None, ""
    return clan, names


def load_pfam_clan_names(path: Path) -> dict[str, str]:
    out, cid = {}, None
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                cid = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF AC") and cid:
                out[line.split()[2].split(".")[0]] = cid
                cid = None
    return out


# ---------------------------------------------------------------- statistics
def enrich(universe: dict[str, set[str]], test_set: set[str], min_universe: int,
           report_ids_fdr: float = 0.05) -> list[dict]:
    """universe: family -> terms (every family has >= 1 term). test_set: families (subset
    of the universe keys are used; others ignored). Returns rows sorted by (p, term)."""
    from scipy.stats import hypergeom
    from stats_util import bh_fdr
    M = len(universe)
    in_set = [f for f in universe if f in test_set]
    N = len(in_set)
    K = Counter(t for terms in universe.values() for t in terms)
    k = Counter(t for f in in_set for t in universe[f])
    fams = defaultdict(list)
    for f in in_set:
        for t in universe[f]:
            fams[t].append(f)
    rows = []
    for term in sorted(K):
        if K[term] < min_universe:
            continue
        kk = k.get(term, 0)
        p = float(hypergeom.sf(kk - 1, M, K[term], N)) if kk > 0 else 1.0
        exp = K[term] * N / M if M else 0.0
        rows.append({"term": term, "k": kk, "n_set": N, "K": K[term], "M": M,
                     "expected": exp, "fold": (kk / exp) if exp else float("nan"), "p": p,
                     "_fams": sorted(fams.get(term, []))})
    for r, q in zip(rows, bh_fdr([r["p"] for r in rows])):
        r["q"] = q
        r["family_ids"] = ";".join(r["_fams"]) if q < report_ids_fdr else ""
        del r["_fams"]
    rows.sort(key=lambda r: (r["p"], r["term"]))
    return rows


def fmt(x) -> str:
    if isinstance(x, float):
        return "NA" if x != x else f"{x:.10g}"
    return str(x)


# ---------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gain-loss-table", required=True, type=Path)
    ap.add_argument("--family-annotation", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--classes", default="shell,soft_core",
                    help="comma list of family classes that form the scope (default shell,soft_core)")
    ap.add_argument("--direction-column", default="rule_freq",
                    help="column of the gain/loss table with the test-set label (default rule_freq)")
    ap.add_argument("--test-sets", default="gain,loss", help="comma list of direction values to test")
    ap.add_argument("--min-universe", type=int, default=5,
                    help="test a term only if >= this many universe families carry it (default 5)")
    ap.add_argument("--report-ids-fdr", type=float, default=0.05)
    ap.add_argument("--pfam-dat", type=Path, default=Path(DEFAULT_PFAM_DAT),
                    help="Pfam-A.hmm.dat (clan mapping and names); Pfam-C next to it gives clan names")
    ap.add_argument("--no-clans", action="store_true", help="skip the pfam_clan vocabulary")
    ap.add_argument("--label", default="")
    args = ap.parse_args(argv)

    classes = [c.strip() for c in args.classes.split(",") if c.strip()]
    test_sets = [c.strip() for c in args.test_sets.split(",") if c.strip()]
    scope = "+".join(classes)
    label = args.label or args.gain_loss_table.resolve().parent.name

    gl = read_gain_loss(args.gain_loss_table, args.direction_column)
    ann, pfam_short = read_annotation(args.family_annotation)
    clans, pfam_names, clan_names = {}, {}, {}
    if not args.no_clans and args.pfam_dat and args.pfam_dat.exists():
        clans, pfam_names = load_pfam_dat(args.pfam_dat)
        pc = args.pfam_dat.parent / "Pfam-C"
        if pc.exists():
            clan_names = load_pfam_clan_names(pc)
    elif not args.no_clans:
        raise SystemExit(f"ERROR: --pfam-dat {args.pfam_dat} not found (use --no-clans to skip clans)")

    in_scope = [r for r in gl if r["class"] in classes]
    missing = [r["family"] for r in in_scope if r["family"] not in ann]
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} scope families have no row in {args.family_annotation} "
                         f"(e.g. {missing[:3]}); every family must be listed (empty terms allowed)")
    dir_of = {r["family"]: r[args.direction_column] for r in in_scope}

    def terms_of(fam: str, v: str) -> set[str]:
        if v == "pfam_clan":
            return {clans[p] for p in ann[fam]["pfam"] if p in clans}
        return ann[fam][v]

    vocabs = [v for v in DOMAIN_VOCABS if not (v == "pfam_clan" and not clans)]
    vocabs += [v for v in OTHER_VOCABS if any(ann[f][v] for f in dir_of)]
    name_of = {"pfam": lambda t: pfam_names.get(t) or pfam_short.get(t, ""),
               "pfam_clan": lambda t: clan_names.get(t, ""),
               "go": lambda t: "", "interpro": lambda t: ""}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    # coverage per direction value in scope
    cov_rows = []
    dvals = sorted(set(dir_of.values()))
    for d in dvals + ["ALL"]:
        fams = [f for f, x in dir_of.items() if d == "ALL" or x == d]
        row = [label, scope, d, len(fams)]
        for v in vocabs:
            n = sum(1 for f in fams if terms_of(f, v))
            row += [n, f"{n / len(fams):.4f}" if fams else "NA"]
        cov_rows.append(row)
    write_tsv(args.out_dir / "coverage.tsv",
              ["label", "scope", "direction", "n_families"]
              + [x for v in vocabs for x in (f"n_with_{v}", f"frac_with_{v}")], cov_rows)

    out_rows = {"domains": [], "go": []}
    test_rows = []
    for ts in test_sets:
        members = {f for f, x in dir_of.items() if x == ts}
        for v in vocabs:
            uni = {f: terms_of(f, v) for f in sorted(dir_of)}
            uni = {f: t for f, t in uni.items() if t}
            rows = enrich(uni, members, args.min_universe, args.report_ids_fdr)
            n_set = sum(1 for f in uni if f in members)
            test_rows.append([label, scope, ts, v, len(uni), n_set, len(members), len(rows),
                              sum(1 for r in rows if r["q"] < 0.05),
                              sum(1 for r in rows if r["q"] < 0.05 and r["k"] >= 3),
                              sum(1 for r in rows if r["q"] < 0.05 and r["fold"] > 1)])
            tgt = out_rows["domains" if v in DOMAIN_VOCABS else "go"]
            for r in rows:
                tgt.append([label, scope, ts, ts, v, r["term"], name_of[v](r["term"]), r["k"],
                            r["n_set"], r["K"], r["M"], fmt(r["expected"]), fmt(r["fold"]),
                            fmt(r["p"]), fmt(r["q"]), r["k"], r["family_ids"]])
    write_tsv(args.out_dir / "enrichment_domains.tsv.gz", COLUMNS, out_rows["domains"])
    if out_rows["go"]:
        write_tsv(args.out_dir / "enrichment_go.tsv.gz", COLUMNS, out_rows["go"])
    write_tsv(args.out_dir / "tests.tsv",
              ["label", "scope", "test_set", "vocab", "M_universe_annotated", "N_set_annotated",
               "n_set_all", "n_terms_tested", "n_fdr05", "n_fdr05_k3", "n_fdr05_fold_gt1"], test_rows)

    import yaml
    from provenance import sha256_of
    inputs = [args.gain_loss_table, args.family_annotation]
    if clans:
        inputs.append(args.pfam_dat)
    prov = {
        "script": "bin/pangenome_gain_loss_enrichment.py",
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "command": " ".join(shlex.quote(a) for a in [sys.argv[0]] + (argv if argv is not None else sys.argv[1:])),
        "date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": platform.node(),
        "label": label,
        "parameters": {"classes": classes, "direction_column": args.direction_column,
                       "test_sets": test_sets, "min_universe": args.min_universe,
                       "report_ids_fdr": args.report_ids_fdr, "vocabularies": vocabs},
        "definitions": {
            "scope": f"families with class in {classes} ({len(dir_of)} families)",
            "universe": "per vocabulary: scope families with >= 1 term of that vocabulary",
            "test_set": f"scope families with {args.direction_column} == <test_set>, inside the universe",
            "statistic": "one-sided hypergeometric P(X >= k) = scipy.stats.hypergeom.sf(k-1, M, K, N); "
                         "BH FDR (lib/stats_util.bh_fdr) within each (test_set, vocab)",
            "annotation_unit": "family representative protein (mmseqs cluster representative); "
                               "other family members are not annotated",
            "go_note": "GO terms used as given, no propagation to parent terms",
        },
        "inputs": [{"path": str(p.resolve()), "sha256": sha256_of(p), "bytes": p.stat().st_size}
                   for p in inputs],
        "outputs": sorted(x.name for x in args.out_dir.iterdir() if x.name != "PROVENANCE.yaml"),
        "license": "derived from the lab's own pipeline output and Pfam (CC0 1.0)",
    }
    with open(args.out_dir / "PROVENANCE.yaml", "w") as fh:
        yaml.safe_dump(prov, fh, sort_keys=False)
    for r in test_rows:
        print("\t".join(str(x) for x in r), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
