#!/usr/bin/env python3
"""Collect the issue #132 tier-1 sweep metrics into one TSV row per grid point.

Each grid point is one pangenome.nf run directory
(``<outdir>/<project>/pangenome/``) produced by ``run_tier1_sweep.sh``.

Per point this reports:

- min_id, cov (parsed from the project name ``..._id<min_id>_cov<cov>``);
- families and core/soft_core/shell/cloud/singleton counts
  (``frequency_table.tsv``, column ``bin``);
- significant islands (data rows of ``significant_islands.tsv``);
- trans pairs (``pair_classification.tsv[.zst]``, ``classification == trans``),
  plus the trans_unconfirmed count for context;
- Leiden module count and largest module size (``module_summary.tsv``);
- AMI, ARI and NMI of protein-level module labels against the baseline point.
  Family IDs are NOT stable across runs (a family is named by its mmseqs
  representative), so labels are never joined on family ID. Instead every
  protein is mapped protein -> family (``cluster/tier1_cluster.tsv``,
  ``rep<TAB>member``) -> module (``family_modules.tsv``). A protein whose
  family is in no module gets the label ``none``. Scores are given twice:
  over all proteins (``none`` is one label), and over only the proteins that
  carry a real module in BOTH runs;
- rho_partial_length for accessory_present vs n50 and vs n_contigs
  (``assembly_quality_correlations.tsv``).

AMI uses the arithmetic-mean normalisation (scikit-learn's default), so the
numbers are comparable with ``sklearn.metrics.adjusted_mutual_info_score``.
The implementation is numpy/scipy only because scikit-learn is not in the
NovInvenio pixi environment; ``test_collect_tier1_sweep.py`` checks it against
reference values.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.special import gammaln

BINS = ("core", "soft_core", "shell", "cloud", "singleton")
NONE_LABEL = "none"


# ---------------------------------------------------------------- file I/O
def open_text(path: Path):
    """Open plain, .gz or .zst text transparently."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if path.suffix == ".zst":
        proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
        return proc.stdout
    return open(path)


def first_existing(*paths: Path) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError(f"none of: {', '.join(str(p) for p in paths)}")


def parse_point(project: str) -> tuple[float, float]:
    m = re.search(r"_id([0-9.]+)_cov([0-9.]+)$", project)
    if not m:
        raise ValueError(f"cannot parse min_id/cov from project name {project!r}")
    return float(m.group(1)), float(m.group(2))


# ---------------------------------------------------------------- label mapping
def read_protein_to_family(cluster_tsv: Path) -> dict[str, str]:
    """protein -> family (mmseqs representative) from a rep<TAB>member TSV."""
    out: dict[str, str] = {}
    with open_text(cluster_tsv) as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            if member in out and out[member] != rep:
                raise ValueError(f"protein {member} in two families in {cluster_tsv}")
            out[member] = rep
    return out


def read_family_to_module(family_modules: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open_text(family_modules) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["family"]] = str(row["module_id"])
    return out


def protein_module_labels(protein_to_family: dict[str, str],
                          family_to_module: dict[str, str]) -> dict[str, str]:
    """protein -> module label; proteins whose family has no module get 'none'."""
    return {p: family_to_module.get(f, NONE_LABEL) for p, f in protein_to_family.items()}


def aligned_label_arrays(labels_a: dict[str, str], labels_b: dict[str, str],
                         drop_none: bool) -> tuple[np.ndarray, np.ndarray]:
    """Integer-code the labels of the proteins present in both runs.

    Proteins missing from one run's cluster TSV are an input error (every run
    clusters the same proteome set), so they raise rather than being dropped.
    """
    if labels_a.keys() != labels_b.keys():
        only_a = len(labels_a.keys() - labels_b.keys())
        only_b = len(labels_b.keys() - labels_a.keys())
        raise ValueError(f"protein sets differ between runs ({only_a} only in A, {only_b} only in B)")
    codes_a: dict[str, int] = {}
    codes_b: dict[str, int] = {}
    xa, xb = [], []
    for p, la in labels_a.items():
        lb = labels_b[p]
        if drop_none and (la == NONE_LABEL or lb == NONE_LABEL):
            continue
        xa.append(codes_a.setdefault(la, len(codes_a)))
        xb.append(codes_b.setdefault(lb, len(codes_b)))
    return np.asarray(xa, dtype=np.int64), np.asarray(xb, dtype=np.int64)


# ---------------------------------------------------------------- partition scores
def contingency(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    na, nb = int(a.max()) + 1, int(b.max()) + 1
    return np.bincount(a * nb + b, minlength=na * nb).reshape(na, nb).astype(np.int64)


def _entropy(counts: np.ndarray) -> float:
    counts = counts[counts > 0].astype(float)
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum())


def _mutual_info(c: np.ndarray) -> float:
    n = c.sum()
    a = c.sum(axis=1)
    b = c.sum(axis=0)
    nz = c > 0
    nij = c[nz].astype(float)
    ai = np.broadcast_to(a[:, None], c.shape)[nz].astype(float)
    bj = np.broadcast_to(b[None, :], c.shape)[nz].astype(float)
    return float((nij / n * (np.log(nij * n) - np.log(ai * bj))).sum())


def _expected_mutual_info(a: np.ndarray, b: np.ndarray, n: int) -> float:
    """E[MI] under the permutation (hypergeometric) model, Vinh et al. 2010."""
    emi = 0.0
    lg_n = gammaln(n + 1)
    for ai in a:
        for bj in b:
            lo = max(1, ai + bj - n)
            hi = min(ai, bj)
            if hi < lo:
                continue
            nij = np.arange(lo, hi + 1, dtype=float)
            term1 = nij / n * (np.log(n * nij) - np.log(float(ai) * float(bj)))
            log_p = (gammaln(ai + 1) + gammaln(bj + 1) + gammaln(n - ai + 1) + gammaln(n - bj + 1)
                     - lg_n - gammaln(nij + 1) - gammaln(ai - nij + 1) - gammaln(bj - nij + 1)
                     - gammaln(n - ai - bj + nij + 1))
            emi += float((term1 * np.exp(log_p)).sum())
    return emi


def ami(a: np.ndarray, b: np.ndarray) -> float:
    """Adjusted mutual information, arithmetic normalisation (sklearn default)."""
    c = contingency(a, b)
    ra, rb = c.sum(axis=1), c.sum(axis=0)
    ra, rb = ra[ra > 0], rb[rb > 0]
    # sklearn's special cases: identical trivial partitions score 1.0
    if (len(ra) == 1 and len(rb) == 1) or (len(ra) == len(rb) == len(a)):
        return 1.0
    n = int(c.sum())
    mi = _mutual_info(c)
    emi = _expected_mutual_info(ra, rb, n)
    ha, hb = _entropy(ra), _entropy(rb)
    denom = (ha + hb) / 2.0 - emi
    if denom < 0:
        denom = min(denom, -np.finfo("float64").eps)
    else:
        denom = max(denom, np.finfo("float64").eps)
    return float((mi - emi) / denom)


def nmi(a: np.ndarray, b: np.ndarray) -> float:
    c = contingency(a, b)
    ha, hb = _entropy(c.sum(axis=1)), _entropy(c.sum(axis=0))
    if ha == 0 and hb == 0:
        return 1.0
    return float(_mutual_info(c) / ((ha + hb) / 2.0))


def ari(a: np.ndarray, b: np.ndarray) -> float:
    c = contingency(a, b).astype(float)
    n = c.sum()

    def comb2(x):
        return x * (x - 1) / 2.0

    sum_ij = comb2(c).sum()
    sum_a = comb2(c.sum(axis=1)).sum()
    sum_b = comb2(c.sum(axis=0)).sum()
    expected = sum_a * sum_b / comb2(n)
    max_idx = (sum_a + sum_b) / 2.0
    if max_idx == expected:
        return 1.0
    return float((sum_ij - expected) / (max_idx - expected))


# ---------------------------------------------------------------- per-run metrics
def count_bins(freq_table: Path) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with open(freq_table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts[row["bin"]] += 1
    unknown = set(counts) - set(BINS)
    if unknown:
        raise ValueError(f"unexpected frequency bins {unknown} in {freq_table}")
    out = {b: counts.get(b, 0) for b in BINS}
    out["families"] = sum(counts.values())
    return out


def count_data_rows(path: Path) -> int:
    with open_text(path) as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def count_pair_classes(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index("classification")
        for line in fh:
            # split only as far as needed
            counts[line.split("\t", idx + 1)[idx]] += 1
    return counts


def module_stats(module_summary: Path) -> tuple[int, int]:
    sizes = []
    with open(module_summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sizes.append(int(row["size"]))
    return len(sizes), (max(sizes) if sizes else 0)


def read_partial_rhos(path: Path) -> dict[str, float]:
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["comparison"]] = float(row["rho_partial_length"])
    return out


def run_labels(pdir: Path) -> dict[str, str]:
    cluster = first_existing(pdir / "cluster" / "tier1_cluster.tsv",
                             pdir / "cluster" / "tier1_cluster.tsv.zst")
    return protein_module_labels(read_protein_to_family(cluster),
                                 read_family_to_module(pdir / "family_modules.tsv"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sweep_dir", required=True, type=Path,
                    help="outdir holding one <project>/pangenome/ per grid point")
    ap.add_argument("--baseline_project", required=True,
                    help="project name of the AMI/ARI reference point (e.g. coccidioides_sweep_id0.9_cov0.8)")
    ap.add_argument("--projects", nargs="*", default=None,
                    help="project names to collect (default: every coccidioides_sweep_* under --sweep_dir)")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args(argv)

    projects = args.projects or sorted(p.name for p in args.sweep_dir.glob("coccidioides_sweep_*")
                                       if (p / "pangenome").is_dir())
    if args.baseline_project not in projects:
        raise SystemExit(f"baseline {args.baseline_project} not among projects {projects}")

    print(f"loading baseline labels from {args.baseline_project}", file=sys.stderr)
    base_labels = run_labels(args.sweep_dir / args.baseline_project / "pangenome")

    cols = ["project", "min_id", "cov", "families", *BINS, "significant_islands",
            "trans_pairs", "trans_unconfirmed_pairs", "leiden_modules", "largest_module",
            "proteins", "proteins_in_module", "ami_all", "ari_all", "nmi_all",
            "n_proteins_both_in_module", "ami_modules_only", "ari_modules_only", "nmi_modules_only",
            "rho_partial_accessory_n50", "rho_partial_accessory_contigs"]
    rows = []
    for proj in projects:
        pdir = args.sweep_dir / proj / "pangenome"
        print(f"collecting {proj}", file=sys.stderr)
        min_id, cov = parse_point(proj)
        row = {"project": proj, "min_id": min_id, "cov": cov}
        row.update(count_bins(pdir / "frequency_table.tsv"))
        row["significant_islands"] = count_data_rows(pdir / "significant_islands.tsv")
        classes = count_pair_classes(first_existing(pdir / "pair_classification.tsv.zst",
                                                    pdir / "pair_classification.tsv"))
        row["trans_pairs"] = classes.get("trans", 0)
        row["trans_unconfirmed_pairs"] = classes.get("trans_unconfirmed", 0)
        row["leiden_modules"], row["largest_module"] = module_stats(pdir / "module_summary.tsv")

        labels = base_labels if proj == args.baseline_project else run_labels(pdir)
        row["proteins"] = len(labels)
        row["proteins_in_module"] = sum(1 for v in labels.values() if v != NONE_LABEL)
        a, b = aligned_label_arrays(base_labels, labels, drop_none=False)
        row["ami_all"], row["ari_all"], row["nmi_all"] = ami(a, b), ari(a, b), nmi(a, b)
        a, b = aligned_label_arrays(base_labels, labels, drop_none=True)
        row["n_proteins_both_in_module"] = len(a)
        if len(a):
            row["ami_modules_only"], row["ari_modules_only"], row["nmi_modules_only"] = ami(a, b), ari(a, b), nmi(a, b)
        else:
            row["ami_modules_only"] = row["ari_modules_only"] = row["nmi_modules_only"] = float("nan")

        rhos = read_partial_rhos(pdir / "assembly_quality_correlations.tsv")
        row["rho_partial_accessory_n50"] = rhos["accessory_present_vs_n50"]
        row["rho_partial_accessory_contigs"] = rhos["accessory_present_vs_n_contigs"]
        rows.append(row)
        if proj != args.baseline_project:
            del labels

    rows.sort(key=lambda r: (r["min_id"], r["cov"]))
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) and k not in ("min_id", "cov") else v)
                        for k, v in r.items()})
    print(f"wrote {len(rows)} rows to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
