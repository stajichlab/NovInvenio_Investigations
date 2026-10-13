#!/usr/bin/env python3
"""Looser loss-direction candidate sets, recomputed offline from the pipeline search cache.

Loss direction of nf_NovInvenio (workflows/loss_search.nf + bin/build_presence_matrix.py
--query-group OUT): an OUT (outgroup) protein is a loss candidate when it is present in
>= outgroup_min_frac of the OUT proteomes (its own proteome counts) and present in
<= loss_ingroup_max_frac of the IN proteomes (0.0 = in none).

This script rebuilds the loss presence matrix from the cached diamond searches
(results/<study>/search_cache/<OUT>_vs_<T>.diamond.tsv.gz) with the pipeline hit
filters, then relaxes the IN side: a protein passes at K when it passes the OUT rule
and has a kept hit in at most K IN proteomes. K = 0 is the pipeline rule and must
reproduce loss_candidates.txt; the symmetric difference is reported.

Hit filters (same order and values as the pipeline at NovInvenio main d03f5a2):
  parse   PARSE_HITS drops rows with evalue > parse_evalue (0.01). Only these rows enter
          the paralog lookup below.
  1       evalue < evalue (1e-5), strict.
  2       paralog competition, scope 'target': drop the hit when the query's
          within-proteome paralog (self_hits/<OUT>.paralog_cutoffs.tsv, OUT proteomes
          only, as LOSS_SEARCH runs PARSE_SELF_HITS on the outgroup) hits the same target
          protein with a strictly lower evalue, unless the query's own evalue <=
          paralog_rescue_evalue (1e-20). paralog_rescue_delta is off (null) in the run.
  3       other-group coverage floor: off in the run (other_coverage_floor_qcov = null),
          so not applied here.
The paralog table is read from the pipeline's published self_hits/ files, not recomputed
from the <OUT>_vs_<OUT> self searches.

The script also cross-checks every cell of the rebuilt matrix against
loss_presence_matrix.tsv, and clusters each candidate set with mmseqs easy-cluster at the
pipeline MMSEQS_CLUSTER settings (--min-seq-id 0.3 -c 0.8 --cov-mode 0) so that the
enrichment script can use family units (--families).

Outputs in --out-dir:
  loss_loose.tsv.gz            one row per OUT protein that passes the OUT rule and has
                               a kept hit in <= max(K) IN proteomes: n_out, n_in, the IN
                               species with a hit, their lineages and best evalues,
                               and the smallest K at which it passes
  candidates_K<k>.txt          Short::protein_id per line (enrichment --candidates input)
  families_K<k>.tsv            mmseqs cluster TSV (rep, member) for that set
  counts.tsv                   proteins and mmseqs families per K
  reproduction.txt             matrix cell check and K=0 vs loss_candidates.txt
  PROVENANCE.yaml              inputs (sha256), parameters, command, tool versions

Run with the NII pixi env python (-I). Heavy reading (about 217 gz files); use sbatch.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
MMSEQS_ARGS = ["--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def filter_pair_hits(rows, paralog_of, evalue=1e-5, parse_evalue=0.01, rescue_evalue=1e-20):
    """Pipeline filters for one query-proteome x target-proteome hit file.

    rows: iterable of (query_id, target_id, evalue). Return {query_id: (best_evalue, target_id)}
    for queries with >= 1 kept hit. Filter 2 uses scope 'target': the paralog's best
    evalue against the same target protein, over rows with evalue <= parse_evalue.
    """
    rows = [(q, t, e) for q, t, e in rows if not e > parse_evalue]
    best = {}
    for q, t, e in rows:
        k = (q, t)
        if k not in best or e < best[k]:
            best[k] = e
    kept = {}
    for q, t, e in rows:
        if not e < evalue:
            continue
        p = paralog_of.get(q)
        if p is not None:
            pev = best.get((p, t))
            if pev is not None and pev < e and not (rescue_evalue and e <= rescue_evalue):
                continue
        if q not in kept or e < kept[q][0]:
            kept[q] = (e, t)
    return kept


def passes(n_out_present, n_out, n_in_present, min_frac, max_in):
    """Loose loss rule: OUT fraction (own proteome counted) >= min_frac and <= max_in IN hits."""
    return n_out_present / n_out >= min_frac and n_in_present <= max_in


def read_fasta_subset(path, keep):
    seqs, cur = {}, None
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                cur = pid if pid in keep else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def mmseqs_cluster(keys, seqs, mmseqs, threads):
    """keys: list of Short::protein_id. Return {member_key: rep_key}."""
    base = os.environ.get("SCRATCH") or os.environ.get("TMPDIR") or None
    tmp = Path(tempfile.mkdtemp(prefix="loss_loose_mmseqs_", dir=base))
    try:
        idx = {}
        with open(tmp / "in.fa", "w") as fh:
            for i, k in enumerate(keys):
                idx[f"s{i}"] = k  # plain ids: mmseqs rewrites tr|ACC|NAME headers
                fh.write(f">s{i}\n{seqs[k]}\n")
        cmd = [mmseqs, "easy-cluster", str(tmp / "in.fa"), str(tmp / "c"), str(tmp / "w"),
               "--threads", str(threads), "-v", "1"] + MMSEQS_ARGS
        subprocess.run(cmd, check=True)
        rep = {}
        with open(tmp / "c_cluster.tsv") as fh:
            for line in fh:
                r, m = line.rstrip("\n").split("\t")[:2]
                rep[idx[m]] = idx[r]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if set(rep) != set(keys):
        sys.exit(f"mmseqs output has {len(rep)} members, expected {len(keys)}")
    return rep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path, help="study config.csv (GROUP, Short, Protein, TaxonGroup)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--proteins-dir", type=Path, default=None, help="default <config dir>/data_dir/pep")
    ap.add_argument("--outgroup-min-frac", type=float, default=0.75)
    ap.add_argument("--ks", default="0,1,2,3", help="comma list of max IN species with a hit")
    ap.add_argument("--in-max-frac", type=float, default=0.15,
                    help="extra set: IN fraction with a hit <= this (reported as its integer K)")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--parse-evalue", type=float, default=0.01)
    ap.add_argument("--rescue-evalue", type=float, default=1e-20)
    ap.add_argument("--mmseqs", default=str(NII_ROOT / ".pixi/envs/default/bin/mmseqs"))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_ON_NODE", "2")))
    args = ap.parse_args()
    t0 = time.time()
    rd, out = args.results, args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    inputs = {}

    cfg = list(csv.DictReader(open(args.config)))
    in_sp = [r["Short"] for r in cfg if r["GROUP"] == "IN"]
    out_sp = [r["Short"] for r in cfg if r["GROUP"] == "OUT"]
    lineage = {r["Short"]: r.get("TaxonGroup", "") for r in cfg}
    pep = {r["Short"]: r["Protein"] for r in cfg}
    inputs[str(args.config)] = sha256(args.config)
    log(f"IN ({len(in_sp)}) OUT ({len(out_sp)})")

    paralog_of = {}
    for o in out_sp:
        f = rd / "self_hits" / f"{o}.paralog_cutoffs.tsv"
        inputs[str(f)] = sha256(f)
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                paralog_of[r["protein_ID"]] = r["paralog_protein_ID"]
    log(f"paralog pairs (OUT proteomes): {len(paralog_of)}")

    # presence[(Q, q)] = {T: (best kept evalue, target)}
    presence = defaultdict(dict)
    n_files = n_rows = 0
    for q_sp in out_sp:
        for t_sp in out_sp + in_sp:
            if t_sp == q_sp:
                continue
            f = rd / "search_cache" / f"{q_sp}_vs_{t_sp}.diamond.tsv.gz"
            if not f.exists():
                sys.exit(f"missing search cache file {f}")
            rows = []
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    c = line.split("\t", 3)
                    rows.append((c[0], c[1], float(c[2])))
            n_files += 1
            n_rows += len(rows)
            for qid, (e, t) in filter_pair_hits(rows, paralog_of, args.evalue, args.parse_evalue,
                                                args.rescue_evalue).items():
                presence[(q_sp, qid)][t_sp] = (e, t)
        log(f"{q_sp}: done ({time.time() - t0:.0f}s, {n_files} files, {n_rows} rows)")
    inputs["search_cache"] = f"{n_files} files {rd / 'search_cache'}/<OUT>_vs_<T>.diamond.tsv.gz (not hashed)"

    # ---- cross-check against the pipeline matrix ----
    mpath = rd / "loss_presence_matrix.tsv"
    inputs[str(mpath)] = sha256(mpath)
    ours = set(presence)
    pipe_rows, cell_diff, cell_n = set(), 0, 0
    diff_examples = []
    with open(mpath) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["source_proteome"], r["protein_id"])
            pipe_rows.add(key)
            mine = presence.get(key, {})
            for s in in_sp + out_sp:
                v = 1 if (s == key[0] or s in mine) else 0
                cell_n += 1
                if v != int(r[s]):
                    cell_diff += 1
                    if len(diff_examples) < 10:
                        diff_examples.append(f"{key[0]}::{key[1]} {s} ours={v} pipeline={r[s]}")
    rows_only_ours = ours - pipe_rows
    rows_only_pipe = pipe_rows - ours

    # ---- loose sets ----
    ks = sorted({int(x) for x in args.ks.split(",") if x != ""})
    k_frac = int(args.in_max_frac * len(in_sp) + 1e-9)  # n_in / n_IN <= f  <=>  n_in <= floor(f * n_IN)
    kmax = max(ks + [k_frac])
    n_out = len(out_sp)
    table = []
    for (q_sp, qid), hits in presence.items():
        n_o = 1 + sum(1 for s in out_sp if s in hits)
        in_hit = [s for s in in_sp if s in hits]
        if not passes(n_o, n_out, len(in_hit), args.outgroup_min_frac, kmax):
            continue
        in_hit.sort(key=lambda s: hits[s][0])
        table.append({
            "candidate": f"{q_sp}::{qid}", "protein_id": qid, "source_proteome": q_sp,
            "n_out_present": n_o, "n_out_total": n_out, "n_in_present": len(in_hit), "n_in_total": len(in_sp),
            "min_k": len(in_hit),
            "in_species_with_hit": ",".join(in_hit),
            "in_lineages": ",".join(sorted({lineage.get(s, "") for s in in_hit})),
            "in_best_evalue": ",".join(f"{hits[s][0]:.3g}" for s in in_hit),
            "in_best_target": ",".join(hits[s][1] for s in in_hit),
        })
    table.sort(key=lambda r: (r["min_k"], r["source_proteome"], r["protein_id"]))

    pipe_c = {l.strip() for l in open(rd / "loss_candidates.txt") if l.strip()}
    inputs[str(rd / "loss_candidates.txt")] = sha256(rd / "loss_candidates.txt")
    sets = {}
    for k in ks:
        sets[f"K{k}"] = (k, sorted(r["candidate"] for r in table if r["n_in_present"] <= k))
    sets[f"INfrac{args.in_max_frac:g}"] = (k_frac, sorted(r["candidate"] for r in table if r["n_in_present"] <= k_frac))
    k0 = set(sets["K0"][1]) if "K0" in sets else set()
    sym = k0 ^ pipe_c

    # ---- families ----
    allkeys = {r["candidate"] for r in table}
    seqs = {}
    pdir = args.proteins_dir or args.config.parent / "data_dir" / "pep"
    by_sp = defaultdict(set)
    for kk in allkeys:
        s, p = kk.split("::", 1)
        by_sp[s].add(p)
    for s, ids in by_sp.items():
        f = pdir / pep[s]
        inputs[str(f)] = sha256(f)
        got = read_fasta_subset(f, ids)
        if set(got) != ids:
            sys.exit(f"{len(ids - set(got))} ids not found in {f}")
        for p, sq in got.items():
            seqs[f"{s}::{p}"] = sq
    count_rows = []
    for name, (k, keys) in sets.items():
        rep = mmseqs_cluster(keys, seqs, args.mmseqs, args.threads) if keys else {}
        with open(out / f"candidates_{name}.txt", "w") as fh:
            fh.write("".join(f"{x}\n" for x in keys))
        with open(out / f"families_{name}.tsv", "w") as fh:
            for m in keys:
                # protein_id only, as the pipeline's clusters/loss_clusters_cluster.tsv
                fh.write(f"{rep[m].split('::', 1)[1]}\t{m.split('::', 1)[1]}\n")
        count_rows.append({"set": name, "max_in_species": k, "n_proteins": len(keys),
                           "n_mmseqs_families": len(set(rep.values())),
                           "n_new_vs_K0": len(set(keys) - k0),
                           "per_out_species": ";".join(
                               f"{o}:{sum(1 for x in keys if x.startswith(o + '::'))}" for o in out_sp)})
    with open(out / "counts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(count_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(count_rows)
    with gzip.open(out / "loss_loose.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(table)

    # pipeline families vs our K0 clustering (partition comparison)
    pipe_fam_path = rd / "clusters" / "loss_clusters_cluster.tsv"
    fam_cmp = "not checked"
    if pipe_fam_path.exists() and "K0" in sets:
        pipe_part = defaultdict(set)
        for line in open(pipe_fam_path):
            r, m = line.rstrip("\n").split("\t")[:2]
            pipe_part[r].add(m)
        our_part = defaultdict(set)
        for line in open(out / "families_K0.tsv"):
            r, m = line.rstrip("\n").split("\t")[:2]
            our_part[r].add(m)
        a = {frozenset(v) for v in pipe_part.values()}
        b = {frozenset(v) for v in our_part.values()}
        fam_cmp = (f"pipeline loss_clusters {len(a)} families; our K0 mmseqs {len(b)} families; "
                   f"identical partitions: {a == b} (shared families {len(a & b)})")

    lines = [
        f"results: {rd}",
        f"search cache files read: {n_files}; rows: {n_rows}",
        f"rebuilt matrix rows (OUT proteins with >= 1 kept hit): {len(ours)}; pipeline loss_presence_matrix.tsv rows: {len(pipe_rows)}",
        f"rows only in rebuilt: {len(rows_only_ours)}; rows only in pipeline: {len(rows_only_pipe)}",
        f"cells compared (pipeline rows x {len(in_sp) + len(out_sp)} species): {cell_n}; differing cells: {cell_diff}",
    ] + [f"  e.g. {x}" for x in diff_examples] + [
        f"rule: OUT present / {n_out} >= {args.outgroup_min_frac} (own proteome counts) and IN present <= K",
        f"K=0 rebuilt: {len(k0)}; loss_candidates.txt: {len(pipe_c)}; symmetric difference: {len(sym)}",
    ] + [f"  {x} {'(rebuilt only)' if x in k0 else '(pipeline only)'}" for x in sorted(sym)[:50]] + [
        f"IN fraction <= {args.in_max_frac:g} of {len(in_sp)} equals IN present <= {k_frac}",
        f"family check: {fam_cmp}",
        f"runtime_s: {time.time() - t0:.0f}",
    ]
    txt = "\n".join(lines) + "\n"
    (out / "reproduction.txt").write_text(txt)
    print(txt)
    print(open(out / "counts.tsv").read())

    prov = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": {"path": str(Path(sys.argv[0]).resolve()), "sha256": sha256(sys.argv[0])},
        "command": " ".join([sys.executable] + sys.argv),
        "class": "ephemeral analysis output (not tracked; regenerate with the command above)",
        "pipeline_reference": "nf_NovInvenio main d03f5a2: workflows/loss_search.nf, bin/build_presence_matrix.py",
        "parameters": {"outgroup_min_frac": args.outgroup_min_frac, "ks": ks, "in_max_frac": args.in_max_frac,
                       "in_max_frac_as_k": k_frac, "evalue": args.evalue, "parse_evalue": args.parse_evalue,
                       "paralog_competition_scope": "target", "paralog_rescue_evalue": args.rescue_evalue,
                       "paralog_rescue_delta": None, "other_coverage_floor_qcov": None,
                       "mmseqs_args": MMSEQS_ARGS},
        "reproduction": {"matrix_rows_rebuilt": len(ours), "matrix_rows_pipeline": len(pipe_rows),
                         "rows_only_rebuilt": len(rows_only_ours), "rows_only_pipeline": len(rows_only_pipe),
                         "cells_compared": cell_n, "cells_differing": cell_diff,
                         "k0_rebuilt": len(k0), "loss_candidates_txt": len(pipe_c),
                         "k0_symmetric_difference": len(sym), "family_check": fam_cmp},
        "counts": count_rows,
        "inputs_sha256": inputs,
        "tools": {"python": platform.python_version(),
                  "mmseqs": subprocess.run([args.mmseqs, "version"], capture_output=True, text=True).stdout.strip()},
        "host": platform.node(), "slurm_job": os.environ.get("SLURM_JOB_ID"),
    }
    try:
        import yaml
        with open(out / "PROVENANCE.yaml", "w") as fh:
            yaml.safe_dump(prov, fh, sort_keys=False, width=120)
    except ImportError:
        json.dump(prov, open(out / "PROVENANCE.yaml", "w"), indent=1)
    return 1 if sym else 0


if __name__ == "__main__":
    raise SystemExit(main())
