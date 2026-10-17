#!/usr/bin/env python3
"""Family-level gain/loss polarity for ALL families of one nf_NovInvenio pangenome run.

The pipeline (nf_NovInvenio bin/pangenome_cooccurrence.py) polarises a family only
when it is `family_a` of an FDR-significant co-occurring pair. This script applies
the SAME outgroup-frequency rule to every family in the run's frequency_table.

Rule (copied from nf_NovInvenio bin/pangenome_cooccurrence.py, identical at
commits 05fcf3a, a2836ad and d03f5a2):
  frac = outgroup strains carrying the family / outgroup strains in the matrix
  polarize_by_frequency():  frac >= loss_min -> loss ; frac <= gain_max -> gain ;
                            otherwise (or no outgroup strain) -> ambiguous
  polarize_direction() (strict): all outgroup strains -> loss ; none -> gain ;
                                 some -> ambiguous
  "present" = matrix cell `present` or `genome_only` (PresenceMatrix.is_present).
  The outgroup denominator is EVERY samplesheet OUT strain that is a matrix column,
  representative or not (pipeline main(), `outgroup_in_matrix`). The matrix is
  presence_matrix.rescued.tsv (workflows/pangenome_profile.nf, COOCCURRENCE input).

Direction column written here:
  conserved  class == core (rule label still reported in rule_freq)
  absent     no ingroup strain carries the family (class outgroup_only / absent)
  gain | loss | ambiguous   rule_freq, for every other class

Inputs (found under --run-dir; .gz / .zst variants are read transparently):
  frequency_table.tsv, presence_matrix.rescued.tsv, samplesheet.with_clades.csv,
  strain_inventory.tsv (optional), clade_assignments.tsv (optional),
  cooccurring_pairs.tsv.zst (optional; used for validation only).

Outputs (--out-dir):
  family_gain_loss.tsv.gz            one row per family
  summary_direction_by_class.tsv     counts, every class x direction (nothing dropped)
  summary_scopes.tsv                 counts per scope (selected classes, + nonrep_only, ...)
  summary_size_distribution.tsv      families per ingroup-representative-strain-count bin
  outgroup_fraction_distribution.tsv families per outgroup present count / fraction bin
  clade_counts.<source>.tsv          per-clade counts (or a "not available" line)
  validation.tsv                     pair family_a labels: pipeline vs this script
  PROVENANCE.yaml

Usage:
  bin/pangenome_gain_loss.py --run-dir studies/fungi/<study>/results/<run> \\
      --out-dir studies/fungi/<study>/analysis/gain_loss/<run>
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import platform
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))

DEFAULT_LOSS_MIN_FRAC = 0.9   # nextflow.config pangenome_polarity_loss_min_frac
DEFAULT_GAIN_MAX_FRAC = 0.1   # nextflow.config pangenome_polarity_gain_max_frac
DEFAULT_CLASSES = ("soft_core", "shell", "cloud", "singleton")
ALL_CLASSES = ("core", "soft_core", "shell", "cloud", "singleton", "nonrep_only",
               "outgroup_only", "ingroup_only", "absent")
DIRECTIONS = ("gain", "loss", "ambiguous", "conserved", "absent")
PRESENT_STATES = ("present", "genome_only")
# Ingroup-representative strain-count bins (label, low, high inclusive).
SIZE_BINS = (("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3-5", 3, 5), ("6-10", 6, 10),
             ("11-20", 11, 20), ("21-50", 21, 50), ("51-100", 51, 100),
             (">100", 101, 10**9))


# ---------------------------------------------------------------- rule (pipeline copy)
def polarize_direction(out_count: int, out_total: int) -> str:
    if out_total == 0 or 0 < out_count < out_total:
        return "ambiguous"
    if out_count == out_total:
        return "loss"
    return "gain"


def outgroup_fraction(out_count: int, out_total: int) -> float | None:
    if out_total == 0:
        return None
    return out_count / out_total


def polarize_by_frequency(out_count: int, out_total: int,
                          loss_min: float = DEFAULT_LOSS_MIN_FRAC,
                          gain_max: float = DEFAULT_GAIN_MAX_FRAC) -> str:
    frac = outgroup_fraction(out_count, out_total)
    if frac is None:
        return "ambiguous"
    if frac >= loss_min:
        return "loss"
    if frac <= gain_max:
        return "gain"
    return "ambiguous"


def validate_thresholds(loss_min: float, gain_max: float) -> None:
    if not (0.0 <= gain_max < loss_min <= 1.0):
        raise ValueError(f"need 0 <= gain_max < loss_min <= 1 (got {gain_max}, {loss_min})")


def direction_for(klass: str, rule_label: str) -> str:
    if klass == "core":
        return "conserved"
    if klass in ("outgroup_only", "absent"):
        return "absent"
    return rule_label


def size_bin(n: int) -> str:
    for label, lo, hi in SIZE_BINS:
        if lo <= n <= hi:
            return label
    raise ValueError(n)


# ---------------------------------------------------------------- I/O helpers
def open_text(path: Path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    if p.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dc", p], stdout=subprocess.PIPE)
        return io.TextIOWrapper(proc.stdout)
    return open(p)


def find_input(pdir: Path, name: str, required: bool = True) -> Path | None:
    for cand in (name, name + ".gz", name + ".zst"):
        if (pdir / cand).exists():
            return pdir / cand
    if name.endswith(".zst") and (pdir / name[:-4]).exists():
        return pdir / name[:-4]
    if required:
        raise SystemExit(f"ERROR: {name} not found in {pdir}")
    return None


def resolve_pangenome_dir(run_dir: Path) -> Path:
    if find_input(run_dir, "frequency_table.tsv", required=False):
        return run_dir
    hits = sorted({p.parent for p in run_dir.glob("*/pangenome/frequency_table.tsv*")}
                  | {p.parent for p in run_dir.glob("pangenome/frequency_table.tsv*")})
    if len(hits) != 1:
        raise SystemExit(f"ERROR: expected one pangenome/frequency_table.tsv under {run_dir}, "
                         f"found {len(hits)}: {hits}")
    return hits[0]


def write_tsv(path: Path, header: list[str], rows) -> None:
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def fmt_frac(x: float | None) -> str:
    return "NA" if x is None else f"{x:.4f}"


# ---------------------------------------------------------------- readers
def read_samplesheet(path: Path) -> list[dict]:
    with open_text(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["stratum_label"] = (r.get("Stratum") or "").strip() or (r.get("TaxonGroup") or "").strip()
    return rows


def read_representatives(path: Path | None) -> set[str] | None:
    if path is None or path.stat().st_size == 0:
        return None
    with open_text(path) as fh:
        return {r["Short"] for r in csv.DictReader(fh, delimiter="\t")
                if r["is_representative"].strip() == "1"}


def read_frequency_table(path: Path) -> dict[str, dict]:
    with open_text(path) as fh:
        return {r["family"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_mash_clades(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with open_text(path) as fh:
        return {r["Short"]: r["mash_clade"] for r in csv.DictReader(fh, delimiter="\t")}


def scan_matrix(path: Path, in_reps: set[str], in_nonreps: set[str], out_all: set[str],
                out_reps: set[str], clade_maps: dict[str, dict[str, str]]):
    """Stream the matrix once. Per family: counts + per-clade rep carriers."""
    res = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        strains = header[1:]
        idx = {s: i + 1 for i, s in enumerate(strains)}
        i_inrep = [idx[s] for s in strains if s in in_reps]
        i_innon = [idx[s] for s in strains if s in in_nonreps]
        i_out = [idx[s] for s in strains if s in out_all]
        i_outrep = [idx[s] for s in strains if s in out_reps]
        clade_cols = {src: [(idx[s], m.get(s, "unknown")) for s in strains if s in in_reps]
                      for src, m in clade_maps.items()}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pres = [p in PRESENT_STATES for p in parts]
            clades = {}
            for src, cols in clade_cols.items():
                c = Counter(lab for i, lab in cols if pres[i])
                clades[src] = dict(c)
            res[parts[0]] = {
                "in_rep": sum(pres[i] for i in i_inrep),
                "in_nonrep": sum(pres[i] for i in i_innon),
                "out_all": sum(pres[i] for i in i_out),
                "out_rep": sum(pres[i] for i in i_outrep),
                "clades": clades,
            }
    n = {"matrix_strains": len(strains), "in_rep": len(i_inrep), "in_nonrep": len(i_innon),
         "out_all": len(i_out), "out_rep": len(i_outrep)}
    return res, n


def read_pair_family_labels(path: Path) -> tuple[dict, int]:
    """Unique family_a -> (direction_a, direction_a_freq, asymmetry_a); conflicts counted."""
    labels: dict[str, tuple] = {}
    conflicts = 0
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ia, js, jf, jx = (header.index(c) for c in
                          ("family_a", "direction_a", "direction_a_freq", "asymmetry_a"))
        for line in fh:
            p = line.rstrip("\n").split("\t")
            val = (p[js], p[jf], p[jx])
            old = labels.setdefault(p[ia], val)
            if old != val:
                conflicts += 1
    return labels, conflicts


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path,
                    help="results dir of one pangenome run (or its pangenome/ subdir)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--loss-min-frac", type=float, default=DEFAULT_LOSS_MIN_FRAC)
    ap.add_argument("--gain-max-frac", type=float, default=DEFAULT_GAIN_MAX_FRAC)
    ap.add_argument("--classes", default=",".join(DEFAULT_CLASSES),
                    help="comma list of ingroup classes for the 'selected' scope "
                         f"(default {','.join(DEFAULT_CLASSES)})")
    ap.add_argument("--outgroup-strains", choices=("all", "representatives"), default="all",
                    help="outgroup denominator: 'all' = pipeline rule (every OUT strain in "
                         "the matrix); 'representatives' = OUT representatives only")
    ap.add_argument("--ingroup-label", default="IN")
    ap.add_argument("--outgroup-label", default="OUT")
    ap.add_argument("--matrix", type=Path, help="override presence matrix path")
    ap.add_argument("--samplesheet", type=Path, help="override samplesheet path")
    ap.add_argument("--label", help="study/run label written to the summaries")
    args = ap.parse_args(argv)
    validate_thresholds(args.loss_min_frac, args.gain_max_frac)
    selected = tuple(c.strip() for c in args.classes.split(",") if c.strip())
    bad = [c for c in selected if c not in ALL_CLASSES]
    if bad:
        ap.error(f"unknown classes {bad}; known: {ALL_CLASSES}")

    pdir = resolve_pangenome_dir(args.run_dir)
    f_freq = find_input(pdir, "frequency_table.tsv")
    f_matrix = args.matrix or find_input(pdir, "presence_matrix.rescued.tsv")
    f_sheet = args.samplesheet or find_input(pdir, "samplesheet.with_clades.csv")
    f_inv = find_input(pdir, "strain_inventory.tsv", required=False)
    f_clades = find_input(pdir, "clade_assignments.tsv", required=False)
    f_pairs = find_input(pdir, "cooccurring_pairs.tsv.zst", required=False)
    label = args.label or f"{args.run_dir.resolve().name}"

    sheet = read_samplesheet(f_sheet)
    reps = read_representatives(f_inv)
    in_all = {r["Short"] for r in sheet if r["GROUP"] == args.ingroup_label}
    out_all = {r["Short"] for r in sheet if r["GROUP"] == args.outgroup_label}
    in_reps = in_all if reps is None else in_all & reps
    in_nonreps = in_all - in_reps
    out_reps = out_all if reps is None else out_all & reps
    clade_maps = {"stratum": {r["Short"]: r["stratum_label"] or "unknown" for r in sheet}}
    mash = read_mash_clades(f_clades)
    if mash:
        clade_maps["mash_clade"] = mash

    freq = read_frequency_table(f_freq)
    mat, n = scan_matrix(f_matrix, in_reps, in_nonreps, out_all, out_reps, clade_maps)
    out_key, out_total = (("out_all", n["out_all"]) if args.outgroup_strains == "all"
                          else ("out_rep", n["out_rep"]))
    rule_text = (f"outgroup_frequency: frac=carriers/{out_total} {args.outgroup_label} strains "
                 f"({args.outgroup_strains}); loss if frac>={args.loss_min_frac}, "
                 f"gain if frac<={args.gain_max_frac}, else ambiguous")

    missing_in_matrix = sorted(set(freq) - set(mat))
    count_mismatch = sum(1 for f, r in freq.items() if f in mat
                         and int(r["strain_count"]) != mat[f]["in_rep"])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fam_rows = []
    for fam in sorted(freq):
        r = freq[fam]
        m = mat.get(fam)
        if m is None:
            continue
        klass = r["bin"]
        oc = m[out_key]
        rf = polarize_by_frequency(oc, out_total, args.loss_min_frac, args.gain_max_frac)
        rs = polarize_direction(oc, out_total)
        d = direction_for(klass, rf)
        used = ("class=core -> conserved" if d == "conserved" else
                "no ingroup carrier -> absent" if d == "absent" else "outgroup_frequency")
        fam_rows.append({
            "family": fam, "class": klass, "ingroup_rep_count": m["in_rep"],
            "ingroup_rep_freq": fmt_frac(m["in_rep"] / n["in_rep"] if n["in_rep"] else None),
            "ingroup_nonrep_count": m["in_nonrep"], "outgroup_count": oc,
            "outgroup_total": out_total, "outgroup_freq": fmt_frac(outgroup_fraction(oc, out_total)),
            "outgroup_rep_freq_table": r.get("frequency_out", "-"),
            "rule_freq": rf, "rule_strict": rs, "direction": d, "rule_used": used,
            "_clades": m["clades"],
        })

    cols = ["family", "class", "ingroup_rep_count", "ingroup_rep_freq", "ingroup_nonrep_count",
            "outgroup_count", "outgroup_total", "outgroup_freq", "outgroup_rep_freq_table",
            "rule_freq", "rule_strict", "direction", "rule_used"]
    raw = args.out_dir / "family_gain_loss.tsv.gz"
    with open(raw, "wb") as bfh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=bfh, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as fh:
                fh.write("\t".join(cols) + "\n")
                for row in fam_rows:
                    fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # direction x class, every class
    dc = Counter((r["class"], r["direction"]) for r in fam_rows)
    rc = Counter((r["class"], r["rule_freq"]) for r in fam_rows)
    classes_seen = [c for c in ALL_CLASSES if any(k[0] == c for k in dc)] + \
        sorted({k[0] for k in dc} - set(ALL_CLASSES))
    write_tsv(args.out_dir / "summary_direction_by_class.tsv",
              ["label", "class", "n_families"] + [f"direction_{d}" for d in DIRECTIONS]
              + [f"rule_freq_{x}" for x in ("gain", "loss", "ambiguous")],
              [[label, c, sum(v for k, v in dc.items() if k[0] == c)]
               + [dc.get((c, d), 0) for d in DIRECTIONS]
               + [rc.get((c, x), 0) for x in ("gain", "loss", "ambiguous")]
               for c in classes_seen])

    # scopes
    pair_labels, pair_conflicts = (read_pair_family_labels(f_pairs) if f_pairs else ({}, 0))
    by_fam = {r["family"]: r for r in fam_rows}
    scopes = [
        ("selected:" + "+".join(selected), lambda r: r["class"] in selected),
        ("selected+nonrep_only", lambda r: r["class"] in selected or r["class"] == "nonrep_only"),
        ("nonrep_only", lambda r: r["class"] == "nonrep_only"),
        ("outgroup_only(rule_freq)", lambda r: r["class"] == "outgroup_only"),
        ("ingroup_only(rule_freq)", lambda r: r["class"] == "ingroup_only"),
        ("core(rule_freq)", lambda r: r["class"] == "core"),
        ("pair_family_a", lambda r: r["family"] in pair_labels),
    ]
    scope_rows = []
    for name, pred in scopes:
        sub = [r for r in fam_rows if pred(r)]
        cf = Counter(r["rule_freq"] for r in sub)
        cs = Counter(r["rule_strict"] for r in sub)
        scope_rows.append([label, name, len(sub)] + [cf.get(x, 0) for x in ("gain", "loss", "ambiguous")]
                          + [cs.get(x, 0) for x in ("gain", "loss", "ambiguous")])
    write_tsv(args.out_dir / "summary_scopes.tsv",
              ["label", "scope", "n_families", "freq_gain", "freq_loss", "freq_ambiguous",
               "strict_gain", "strict_loss", "strict_ambiguous"], scope_rows)

    # size distribution (selected + nonrep_only, by class, so either scope can be rebuilt)
    sz = Counter((r["class"], size_bin(r["ingroup_rep_count"]), r["direction"]) for r in fam_rows
                 if r["class"] in selected or r["class"] == "nonrep_only")
    size_rows = []
    for c in [c for c in ALL_CLASSES if c in selected or c == "nonrep_only"]:
        for b, _, _ in SIZE_BINS:
            vals = [sz.get((c, b, d), 0) for d in ("gain", "loss", "ambiguous")]
            if sum(vals):
                size_rows.append([label, c, b] + vals)
    write_tsv(args.out_dir / "summary_size_distribution.tsv",
              ["label", "class", "ingroup_rep_strains_bin", "gain", "loss", "ambiguous"], size_rows)

    # outgroup fraction distribution (selected classes)
    sel = [r for r in fam_rows if r["class"] in selected]
    if out_total <= 10:
        key = lambda r: f"{r['outgroup_count']}/{out_total}"  # noqa: E731
        order = [f"{k}/{out_total}" for k in range(out_total + 1)]
    else:
        def key(r):
            f = r["outgroup_count"] / out_total
            return ("0" if f == 0 else "(0,gain_max]" if f <= args.gain_max_frac else
                    "(gain_max,loss_min)" if f < args.loss_min_frac else
                    "[loss_min,1)" if f < 1 else "1")
        order = ["0", "(0,gain_max]", "(gain_max,loss_min)", "[loss_min,1)", "1"]
    of = Counter((key(r), r["rule_freq"]) for r in sel)
    write_tsv(args.out_dir / "outgroup_fraction_distribution.tsv",
              ["label", "scope", "outgroup_carriers", "n_families", "rule_freq"],
              [[label, "selected", k, of[(k, d)], d] for k in order
               for d in ("gain", "loss", "ambiguous") if of.get((k, d))])

    # clade counts
    clade_files = []
    for src in clade_maps:
        rows = []
        pres = defaultdict(Counter)
        priv = defaultdict(Counter)
        for r in sel:
            cl = r["_clades"].get(src, {})
            for c in cl:
                pres[c][r["direction"]] += 1
            if len(cl) == 1:
                priv[next(iter(cl))][r["direction"]] += 1
        n_reps_clade = Counter(clade_maps[src].get(s, "unknown") for s in in_reps)
        for c in sorted(n_reps_clade):
            rows.append([label, src, c, n_reps_clade[c]]
                        + [pres[c].get(d, 0) for d in ("gain", "loss", "ambiguous")]
                        + [priv[c].get(d, 0) for d in ("gain", "loss", "ambiguous")])
        fn = args.out_dir / f"clade_counts.{src}.tsv"
        write_tsv(fn, ["label", "clade_source", "clade", "n_ingroup_reps",
                       "present_gain", "present_loss", "present_ambiguous",
                       "private_gain", "private_loss", "private_ambiguous"], rows)
        clade_files.append(fn.name)

    # validation against the pipeline pair labels
    val_rows = []
    if pair_labels:
        miss = [f for f in pair_labels if f not in by_fam]
        for name, pi, mine in (("direction_a_freq", 1, "rule_freq"), ("direction_a", 0, "rule_strict")):
            pc = Counter(v[pi] for v in pair_labels.values())
            mc = Counter(by_fam[f][mine] for f in pair_labels if f in by_fam)
            nmis = sum(1 for f, v in pair_labels.items() if f in by_fam and by_fam[f][mine] != v[pi])
            for x in ("gain", "loss", "ambiguous"):
                val_rows.append([label, name, x, pc.get(x, 0), mc.get(x, 0)])
            val_rows.append([label, name, "family_mismatches", 0, nmis])
        nfrac = sum(1 for f, v in pair_labels.items() if f in by_fam
                    and v[2] != by_fam[f]["outgroup_freq"])
        val_rows += [[label, "asymmetry_a", "4dp_mismatches", 0, nfrac],
                     [label, "pairs", "unique_family_a", len(pair_labels), len(pair_labels) - len(miss)],
                     [label, "pairs", "family_a_label_conflicts", pair_conflicts, ""]]
    else:
        val_rows.append([label, "pairs", "not_available", "", ""])
    write_tsv(args.out_dir / "validation.tsv",
              ["label", "column", "value", "pipeline", "this_script"], val_rows)

    # provenance
    import yaml
    from provenance import sha256_of
    inputs = [p for p in (f_freq, f_matrix, f_sheet, f_inv, f_clades, f_pairs) if p]
    prov = {
        "script": "bin/pangenome_gain_loss.py",
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "command": " ".join(shlex.quote(a) for a in [sys.argv[0]] + (argv or sys.argv[1:])),
        "date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": platform.node(),
        "derived_by": "bin/pangenome_gain_loss.py (outgroup-frequency rule copied from "
                      "nf_NovInvenio bin/pangenome_cooccurrence.py polarize_by_frequency)",
        "parameters": {"loss_min_frac": args.loss_min_frac, "gain_max_frac": args.gain_max_frac,
                       "selected_classes": list(selected), "outgroup_strains": args.outgroup_strains,
                       "ingroup_label": args.ingroup_label, "outgroup_label": args.outgroup_label},
        "rule": rule_text,
        "strain_counts": n,
        "checks": {"families_in_frequency_table": len(freq),
                   "families_missing_from_matrix": len(missing_in_matrix),
                   "ingroup_rep_count_mismatch_vs_frequency_table": count_mismatch,
                   "pair_family_a_label_conflicts": pair_conflicts},
        "clade_tables": clade_files or "not available",
        "inputs": [{"path": str(p.resolve()), "sha256": sha256_of(p),
                    "bytes": p.stat().st_size} for p in inputs],
        "outputs": sorted(x.name for x in args.out_dir.iterdir() if x.name != "PROVENANCE.yaml"),
        "license": "derived from the lab's own pipeline output; no external redistribution terms",
    }
    with open(args.out_dir / "PROVENANCE.yaml", "w") as fh:
        yaml.safe_dump(prov, fh, sort_keys=False)
    print(f"{label}: {len(fam_rows)} families; strains {n}; rule: {rule_text}", file=sys.stderr)
    if count_mismatch or missing_in_matrix:
        print(f"WARNING: {count_mismatch} ingroup rep-count mismatches vs frequency_table; "
              f"{len(missing_in_matrix)} families missing from matrix", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
