#!/usr/bin/env python3
"""Over-representation of Pfam / Pfam-clan / InterPro (primary) and GO (secondary)
terms among the novelty candidates of one nf_NovInvenio pairwise run.

Version 2 of studies/fungi/neocallimastigomycota_vs_early_fungi/analysis/enrichment/
enrich_novelties.py (that prototype is unchanged). Changes in v2:

  1. Universe choice (--universe):
       proteome     every protein in the IN proteome FASTAs (default)
       matrix       rows of presence_matrix.uniprot.tsv (IN proteins with >=1 kept hit;
                    the prototype's universe)
       in_conserved matrix rows that pass the pipeline's ingroup rule
                    (IN count / n_IN >= --ingroup-min-frac), whatever their OUT count.
                    Candidates are the subset with 0 OUT; this asks "which terms are
                    over-represented among IN-conserved proteins that have no OUT hit".
  2. Annotation source (--annotations):
       matrix       uniprot_* columns of presence_matrix.uniprot.tsv (only valid
                    with --universe matrix or in_conserved)
       <DIR>        per-proteome UniProt TSVs <DIR>/<Protein_Accession>.tsv.gz
                    (columns accession, go_ids, pfam_ids, pfam_names, interpro_ids,
                    description), mapped by --species-csv (Short -> Protein_Accession).
                    Species without a file get no annotation; coverage.tsv shows it.
  3. Background clustering (--cluster-background): non-candidate IN proteins of the
     whole proteome are clustered once with mmseqs easy-cluster at the pipeline's
     MMSEQS_CLUSTER settings (--min-seq-id 0.3 -c 0.8 --cov-mode 0, default
     sensitivity and cluster mode). The result is cached (--cache-dir), keyed by the
     sha256 of the background FASTA plus the mmseqs settings. Family level then uses
     clustered units on both sides:
       candidates  pipeline family_id from novelties.<Short>.tsv (mmseqs of the
                   candidates only, same settings); empty family_id = singleton
       background  mmseqs cluster of the non-candidate IN proteins
     Scope A: one unit per family / cluster across species.
     Scope B: one unit per family / cluster restricted to that species.
     A unit's terms are the union of its members' terms. Without the flag the
     background stays one unit per protein (prototype behaviour).
  4. Outputs split by role: enrichment_domains.tsv.gz (pfam, pfam_clan, interpro =
     primary) and enrichment_go.tsv.gz (secondary). GO rows carry go_evidence:
     "IEA_only" if every candidate-unit annotation behind the term (after true-path
     propagation) is IEA, else "has_non_IEA", plus n_cand_units_non_IEA and
     cand_evidence_codes (candidate units per evidence code). Note: IBA (phylogenetic
     inference) counts as non-IEA but is also not experimental evidence.
  5. Optional --controls TSV: presence counts of named domain sets among all IN
     proteins, matrix rows, candidates and near-misses (ingroup rule passes, present
     in 1 or 2 OUT proteomes per presence_matrix.uniprot.tsv). Counts only.
  6. Writes PROVENANCE.yaml in --out-dir (inputs with sha256, tool versions, mmseqs
     settings, go-basic.obo data-version, Pfam .dat path, command, runtimes).
  7. --direction {novelty,loss} (added 2026-10-12; default novelty, which gives the
     same tests.tsv as before the flag existed). The "target" group is the group the
     candidates come from:
       novelty  target = IN,  other = OUT; candidates.txt; presence_matrix.uniprot.tsv;
                families from novelties.<Short>.tsv family_id
       loss     target = OUT, other = IN;  loss_candidates.txt; loss_presence_matrix.tsv;
                families from clusters/loss_clusters_cluster.tsv (the pipeline's
                LOSS_CLUSTER mmseqs run on loss_candidates.fa, same settings)
     Everything that the text above says about "IN" applies to the target group:
     universe = target-group proteomes, background = non-candidate target proteins,
     scope B = each target species. The matrix rule check is the pipeline rule of
     bin/build_presence_matrix.py: target count / n_target >= --ingroup-min-frac
     (in the loss direction pass the pipeline's outgroup_min_frac) and other count /
     n_other <= --other-max-frac (pipeline loss_ingroup_max_frac; 0.0 = absent).
     --candidates / --presence-matrix / --families override the default files.
     tests.tsv has the same columns in both directions; the direction is recorded
     in PROVENANCE.yaml and in the output directory name.
  8. Mixed annotation sources (added 2026-10-14). An annotation file can also come from
     a Pfam-A HMM search (bin/hmm_domtbl_to_annotations.py, for NCBI/JGI proteomes); such
     a file has Pfam only, with empty GO and InterPro columns. Its species then adds
     units to the Pfam and Pfam-clan tests only (the universe of each test is the units
     with >= 1 term in that vocabulary). PROVENANCE.yaml records, per vocabulary, the
     target species with and without terms (vocabulary_coverage) and the source of each
     annotation file read from its .provenance.yaml sidecar (annotation_sources). The
     log prints a NOTE line for each vocabulary that some species lack. coverage.tsv is
     unchanged: its per-species n_<vocab> columns already give the counts.

Statistics (same as the prototype): one-sided hypergeometric P(X >= k)
(scipy.stats.hypergeom.sf(k-1, M, n, N)). Universe per test = units with >=1 term
in that vocabulary. Terms with < --min-bg universe units are not tested. BH-FDR
(lib/stats_util.bh_fdr) within each test = scope x level x vocabulary -> fdr_bh.
Each vocabulary is corrected on its own, so the large GO test count does not change
Pfam q-values. Scope B also gets fdr_bh_all_species: BH across all species tests of
the same level x vocabulary. goatools is used only to read go-basic.obo and get
is_a ancestors; no goatools test is run.

Run with the NII pixi env: .pixi/envs/default/bin/python3 enrich_novelties_v2.py ...
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from scipy.stats import hypergeom

# Hardcoded per CLAUDE.md (study-specific script location; no __file__-relative paths).
NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from stats_util import bh_fdr  # noqa: E402

csv.field_size_limit(10**9)
GO_ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
DOMAIN_VOCABS = ["pfam", "pfam_clan", "interpro"]
VOCABS = DOMAIN_VOCABS + ["go"]
MMSEQS_ARGS = ["--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0"]
T0 = time.time()
TIMES = {}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def split_ids(s):
    return [x for x in s.split("|") if x] if s else []


def parse_go(s, exclude_iea):
    """Return {go_id: set(evidence codes)}."""
    out = defaultdict(set)
    for entry in split_ids(s):
        go_id, _, ev = entry.rpartition(":")
        if not go_id.startswith("GO:"):
            go_id, ev = entry, ""
        if exclude_iea and ev == "IEA":
            continue
        out[go_id].add(ev)
    return out


def bare_acc(pid):
    parts = pid.split("|")
    if len(parts) >= 2 and parts[0] in ("tr", "sp"):
        return parts[1]
    return pid


def read_fasta(path):
    """Yield (id, seq). id = first header token."""
    pid, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(seq)
                pid, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if pid is not None:
        yield pid, "".join(seq)


def load_pfam_dat(path):
    """Return (family -> clan, family -> 'ID: DE')."""
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


def load_pfam_clan_names(path):
    out, cid = {}, None
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                cid = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF AC") and cid:
                out[line.split()[2].split(".")[0]] = cid
                cid = None
    return out


def rule_candidates(counts, n_target, n_other, min_frac, other_max_frac=0.0):
    """Pipeline candidate rule (nf_NovInvenio bin/build_presence_matrix.py, both directions).

    counts: {key: (target_count, other_count)} from a presence matrix (the protein's own
    proteome counts in target_count). Return (target_rule, candidates):
      target_rule = keys with target_count / n_target >= min_frac
      candidates  = target_rule keys with other_count / n_other <= other_max_frac
    """
    target_rule = {k for k, (tc, _) in counts.items() if tc / n_target >= min_frac}
    cands = {k for k in target_rule
             if (counts[k][1] / n_other if n_other else 0) <= other_max_frac}
    return target_rule, cands


def vocabulary_coverage(cov_rows, vocabs=VOCABS):
    """Per vocabulary: which target species have >= 1 annotated protein in the universe.

    cov_rows: coverage.tsv rows (species, subset, n_proteins, n_<vocab>). ALL_* rows are
    ignored. Returns {vocab: {"species_with_terms": [...], "species_without_terms": [...],
    "n_with": int, "n_species": int}}. A species without terms for a vocabulary adds
    no units to that vocabulary's tests (for example, a Pfam-only hmmsearch annotation
    has no GO or InterPro).
    """
    sp_rows = [r for r in cov_rows if r["subset"] == "universe" and not str(r["species"]).startswith("ALL_")
               and int(r["n_proteins"]) > 0]
    out = {}
    for v in vocabs:
        with_t = [r["species"] for r in sp_rows if int(r[f"n_{v}"]) > 0]
        without = [r["species"] for r in sp_rows if int(r[f"n_{v}"]) == 0]
        out[v] = {"species_with_terms": with_t, "species_without_terms": without,
                  "n_with": len(with_t), "n_species": len(sp_rows)}
    return out


def tool_version(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip().splitlines()[0]
    except Exception as e:  # noqa: BLE001
        return f"unknown ({e})"


def cluster_background(bg_seqs, cache_dir, mmseqs, threads):
    """bg_seqs: list of (key, seq) sorted by key. Return {key: cluster_rep_key}, meta."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    for k, s in bg_seqs:
        h.update(f">{k}\n{s}\n".encode())
    h.update(" ".join(MMSEQS_ARGS).encode())
    key = h.hexdigest()[:16]
    cached = cache_dir / f"bg_clusters.{key}.tsv.gz"
    meta_path = cache_dir / f"bg_clusters.{key}.json"
    if cached.exists() and meta_path.exists():
        log(f"background clusters: cache hit {cached}")
        rep = {}
        with gzip.open(cached, "rt") as fh:
            for line in fh:
                r, m = line.rstrip("\n").split("\t")
                rep[m] = r
        meta = json.load(open(meta_path))
        meta["cache_hit"] = True
        return rep, meta
    scratch = os.environ.get("SCRATCH") or os.environ.get("TMPDIR") or None
    tmp = Path(tempfile.mkdtemp(prefix="bgclust_", dir=scratch))
    idx = {}
    fa = tmp / "bg.fa"
    with open(fa, "w") as fh:
        for i, (k, s) in enumerate(bg_seqs):
            sid = f"p{i:07d}"  # plain ids: avoids mmseqs header rewriting of tr|ACC|NAME
            idx[sid] = k
            fh.write(f">{sid}\n{s}\n")
    cmd = [mmseqs, "easy-cluster", str(fa), str(tmp / "bg"), str(tmp / "mmtmp"),
           "--threads", str(threads)] + MMSEQS_ARGS
    log("running: " + " ".join(cmd))
    t = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"mmseqs failed:\n{res.stderr[-3000:]}")
    secs = round(time.time() - t, 1)
    rep = {}
    with open(tmp / "bg_cluster.tsv") as fh, gzip.open(cached, "wt") as out:
        for line in fh:
            r, m = line.rstrip("\n").split("\t")
            rep[idx[m]] = idx[r]
            out.write(f"{idx[r]}\t{idx[m]}\n")
    if len(rep) != len(bg_seqs):
        sys.exit(f"mmseqs cluster tsv has {len(rep)} members, expected {len(bg_seqs)}")
    meta = {"key": key, "n_sequences": len(bg_seqs), "n_clusters": len(set(rep.values())),
            "mmseqs": mmseqs, "mmseqs_version": tool_version([mmseqs, "version"]),
            "command": "mmseqs easy-cluster <bg.fa> <out> <tmp> --threads N " + " ".join(MMSEQS_ARGS),
            "runtime_s": secs, "threads": threads,
            "created": datetime.datetime.now().isoformat(timespec="seconds")}
    json.dump(meta, open(meta_path, "w"), indent=1)
    shutil.rmtree(tmp, ignore_errors=True)
    meta["cache_hit"] = False
    return rep, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, type=Path, help="pipeline results dir (results/<study>)")
    ap.add_argument("--config", required=True, type=Path, help="study config.csv (GROUP, Short, Protein)")
    ap.add_argument("--annotations", required=True,
                    help="'matrix' or a dir of <Protein_Accession>.tsv.gz UniProt annotation files")
    ap.add_argument("--species-csv", type=Path, default=None,
                    help="species.csv with Short, Protein_Accession (default: <config dir>/species.csv)")
    ap.add_argument("--proteins-dir", type=Path, default=None,
                    help="protein FASTAs named as in config Protein (default: <config dir>/data_dir/pep)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--study", required=True, help="study name (labels and provenance only)")
    ap.add_argument("--universe", choices=["proteome", "matrix", "in_conserved"], default="proteome")
    ap.add_argument("--direction", choices=["novelty", "loss"], default="novelty",
                    help="novelty: candidates from IN, absent from OUT (default). "
                         "loss: candidates from OUT, absent from IN (see docstring item 7)")
    ap.add_argument("--candidates", type=Path, default=None,
                    help="candidate list (Short::protein_id per line); default "
                         "<results>/candidates.txt (novelty) or <results>/loss_candidates.txt (loss)")
    ap.add_argument("--presence-matrix", type=Path, default=None,
                    help="presence matrix; default <results>/presence_matrix.uniprot.tsv (novelty) "
                         "or <results>/loss_presence_matrix.tsv (loss)")
    ap.add_argument("--families", type=Path, default=None,
                    help="loss direction only: mmseqs cluster TSV (rep, member) of the candidates; "
                         "default <results>/clusters/loss_clusters_cluster.tsv")
    ap.add_argument("--ingroup-min-frac", type=float, default=0.75,
                    help="presence threshold within the target group (pipeline ingroup_min_frac for "
                         "novelty, outgroup_min_frac for loss); in_conserved universe and near-miss rule")
    ap.add_argument("--other-max-frac", type=float, default=0.0,
                    help="max fraction of the other group with a hit for the rule check (pipeline "
                         "--other-max-frac; loss_ingroup_max_frac in the loss direction). Default 0.0")
    ap.add_argument("--cluster-background", action="store_true")
    ap.add_argument("--cache-dir", type=Path, default=None, help="default: <out-dir>/cache")
    ap.add_argument("--mmseqs", default=str(NII_ROOT / ".pixi/envs/default/bin/mmseqs"))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_ON_NODE", "4")))
    ap.add_argument("--obo", type=Path, default=NII_ROOT / "data/go/go-basic.obo")
    ap.add_argument("--pfam-dat", type=Path, default=None, help="Pfam-A.hmm.dat (clans, names); Pfam-C next to it gives clan names")
    ap.add_argument("--interpro-entries", type=Path, default=None, help="InterProScan entries.json (IPR names)")
    ap.add_argument("--controls", type=Path, default=None,
                    help="TSV: control, pfam_ids ('|'-separated), desc_regex (optional, case-insensitive)")
    ap.add_argument("--min-bg", type=int, default=5)
    ap.add_argument("--exclude-iea", action="store_true")
    ap.add_argument("--report-ids-fdr", type=float, default=0.05)
    args = ap.parse_args()

    rd = args.results
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    cfg_dir = args.config.parent
    species_csv = args.species_csv or cfg_dir / "species.csv"
    prot_dir = args.proteins_dir or cfg_dir / "data_dir" / "pep"
    cache_dir = args.cache_dir or out / "cache"
    use_matrix_ann = args.annotations == "matrix"
    if use_matrix_ann and args.universe == "proteome":
        sys.exit("--annotations matrix cannot annotate the whole proteome; give an annotation dir or --universe matrix")
    inputs = {}  # path -> sha256

    # ---------- config ----------
    cfg = list(csv.DictReader(open(args.config)))
    loss = args.direction == "loss"
    # From here on in_sp = TARGET group (candidates' source) and out_sp = OTHER group.
    # Novelty: IN / OUT as named. Loss: swapped (OUT is the target group).
    tgt_label, oth_label = ("OUT", "IN") if loss else ("IN", "OUT")
    in_sp = [r["Short"] for r in cfg if r["GROUP"] == tgt_label]
    out_sp = [r["Short"] for r in cfg if r["GROUP"] == oth_label]
    pep = {r["Short"]: r["Protein"] for r in cfg}
    inputs[str(args.config)] = sha256(args.config)
    log(f"direction={args.direction} target({tgt_label})={in_sp} other({oth_label})={out_sp}")

    # ---------- candidates, families ----------
    cand_path = args.candidates or rd / ("loss_candidates.txt" if loss else "candidates.txt")
    candidates = {l.strip() for l in open(cand_path) if l.strip()}
    inputs[str(cand_path)] = sha256(cand_path)
    family = {}
    fam_source = None
    if loss:
        fpath = args.families or rd / "clusters" / "loss_clusters_cluster.tsv"
        inputs[str(fpath)] = sha256(fpath)
        fam_source = str(fpath)
        key_of = defaultdict(list)
        for c in candidates:
            key_of[c.split("::", 1)[1]].append(c)
        with open(fpath) as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")[:2]
                for key in key_of.get(mem, []):
                    family[key] = rep
    else:
        fam_source = "novelties.<Short>.tsv family_id"
        for sp in in_sp:
            f = rd / f"novelties.{sp}.tsv"
            if not f.exists():
                continue
            inputs[str(f)] = sha256(f)
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    key = f"{row['source_proteome']}::{row['protein_id']}"
                    family[key] = row.get("family_id") or f"singleton:{key}"
    nofam = [c for c in candidates if c not in family]
    if nofam:
        log(f"WARNING: {len(nofam)} candidates have no family row; treated as singletons")

    # ---------- matrix ----------
    mpath = args.presence_matrix or rd / ("loss_presence_matrix.tsv" if loss else "presence_matrix.uniprot.tsv")
    inputs[str(mpath)] = sha256(mpath)
    matrix = {}  # key -> (target_count, other_count, row-annotation dict)
    with open(mpath) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if use_matrix_ann and "uniprot_pfam_ids" not in (rdr.fieldnames or []):
            sys.exit(f"--annotations matrix needs uniprot_* columns; {mpath} has none "
                     "(use the .uniprot.tsv matrix via --presence-matrix)")
        miss_cols = [s for s in in_sp + out_sp if s not in (rdr.fieldnames or [])]
        if miss_cols:
            sys.exit(f"config species not in {mpath}: {miss_cols}")
        for row in rdr:
            key = f"{row['source_proteome']}::{row['protein_id']}"
            ic = sum(int(row[s] or 0) for s in in_sp)
            oc = sum(int(row[s] or 0) for s in out_sp)
            ann = None
            if use_matrix_ann:
                ann = {"go": row["uniprot_go_ids"], "pfam": row["uniprot_pfam_ids"],
                       "pfam_names": row["uniprot_pfam_names"], "interpro": row["uniprot_interpro_ids"],
                       "desc": row.get("uniprot_description", "")}
            matrix[key] = (ic, oc, ann)
    in_rule, rule_cands = rule_candidates({k: (v[0], v[1]) for k, v in matrix.items()},
                                          len(in_sp), len(out_sp), args.ingroup_min_frac,
                                          args.other_max_frac)
    log(f"matrix rows {len(matrix)}; target rule {len(in_rule)}; target rule & other frac <= "
        f"{args.other_max_frac} {len(rule_cands)}; {cand_path.name} {len(candidates)}")
    # Key names kept from the novelty-only version; in the loss direction "in_rule" means the
    # target-group (OUT) rule and "0_out" means 0 hits in the other group (IN).
    cand_check = {"matrix_in_rule_and_0_out": len(rule_cands), "candidates_txt": len(candidates),
                  "symmetric_difference": len(rule_cands ^ candidates),
                  "candidates_file": str(cand_path), "presence_matrix": str(mpath),
                  "rule": (f"{tgt_label} count / {len(in_sp)} >= {args.ingroup_min_frac} and "
                           f"{oth_label} count / {len(out_sp)} <= {args.other_max_frac} "
                           "(own proteome counts toward the target group)"),
                  "only_in_matrix_rule": sorted(rule_cands - candidates)[:50],
                  "only_in_candidates_file": sorted(candidates - rule_cands)[:50]}

    # ---------- proteomes ----------
    seqs = {}
    sp_of = {}
    need_fasta = args.universe == "proteome" or args.cluster_background
    if need_fasta:
        for sp in in_sp:
            f = prot_dir / pep[sp]
            inputs[str(f)] = sha256(f)
            for pid, s in read_fasta(f):
                k = f"{sp}::{pid}"
                seqs[k] = s
                sp_of[k] = sp
    if args.universe == "proteome":
        universe = set(seqs)
    elif args.universe == "matrix":
        universe = set(matrix)
    else:
        universe = set(in_rule)
    for k in universe:
        sp_of.setdefault(k, k.split("::", 1)[0])
    miss = candidates - universe
    if miss:
        sys.exit(f"ERROR: {len(miss)} candidates not in universe '{args.universe}', e.g. {sorted(miss)[:3]}")
    log(f"universe '{args.universe}': {len(universe)} {tgt_label} (target group) proteins")

    # ---------- annotations ----------
    ann_raw = {}  # key -> dict
    ann_files = {}
    all_sp_ann = defaultdict(dict)  # species -> acc -> ann (all config species; for controls)
    if not use_matrix_ann:
        adir = Path(args.annotations)
        upid = {r["Short"]: r.get("Protein_Accession", "") for r in csv.DictReader(open(species_csv))}
        inputs[str(species_csv)] = sha256(species_csv)
        for sp in in_sp + out_sp:
            f = adir / f"{upid.get(sp, '')}.tsv.gz"
            if not upid.get(sp) or not f.exists():
                log(f"no annotation file for {sp} ({f})")
                continue
            ann_files[sp] = str(f)
            inputs[str(f)] = sha256(f)
            with gzip.open(f, "rt") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    all_sp_ann[sp][row["accession"]] = {
                        "go": row.get("go_ids", ""), "pfam": row.get("pfam_ids", ""),
                        "pfam_names": row.get("pfam_names", ""), "interpro": row.get("interpro_ids", ""),
                        "desc": row.get("description", "")}
        for k in universe:
            sp, pid = k.split("::", 1)
            a = all_sp_ann.get(sp, {}).get(bare_acc(pid))
            if a:
                ann_raw[k] = a
    else:
        for k in universe:
            ann_raw[k] = matrix[k][2]
    log(f"annotation rows matched: {len(ann_raw)} of {len(universe)}")

    # ---------- vocab helpers ----------
    clans, pfam_names = ({}, {})
    names = defaultdict(dict)
    if args.pfam_dat:
        clans, pfam_names = load_pfam_dat(args.pfam_dat)
        names["pfam"].update(pfam_names)
        pc = args.pfam_dat.parent / "Pfam-C"
        if pc.exists():
            names["pfam_clan"].update(load_pfam_clan_names(pc))
    from goatools.obo_parser import GODag
    dag = GODag(str(args.obo), prt=None)
    obo_version = ""
    with open(args.obo) as fh:
        for line in fh:
            if line.startswith("data-version:"):
                obo_version = line.split(":", 1)[1].strip()
                break
    for t, rec in dag.items():
        names["go"][t] = f"{rec.namespace}:{rec.name}"
    if args.interpro_entries:
        ent = json.load(open(args.interpro_entries))
        for k, v in ent.items():
            if k.startswith("IPR"):
                names["interpro"][k] = f"{v.get('type', '')}:{v.get('description') or v.get('name', '')}"
    anc_cache = {}

    def ancestors(t):
        if t not in anc_cache:
            rec = dag.get(t)
            anc_cache[t] = ({rec.id} | rec.get_all_parents()) - GO_ROOTS if rec is not None else set()
        return anc_cache[t]

    prot = {}
    for k in universe:
        a = ann_raw.get(k)
        pf = set(split_ids(a["pfam"])) if a else set()
        if a and not args.pfam_dat:
            pfn = split_ids(a["pfam_names"])
            if len(pfn) == len(split_ids(a["pfam"])):
                for x, n in zip(split_ids(a["pfam"]), pfn):
                    names["pfam"].setdefault(x, n)
        go_ev = defaultdict(set)
        if a:
            for t, evs in parse_go(a["go"], args.exclude_iea).items():
                for anc in ancestors(t):
                    go_ev[anc] |= evs
        prot[k] = {"sp": sp_of[k], "cand": k in candidates,
                   "pfam": pf, "pfam_clan": {clans[p] for p in pf if p in clans},
                   "interpro": set(split_ids(a["interpro"])) if a else set(),
                   "go": dict(go_ev), "go_raw_n": len(split_ids(a["go"])) if a else 0,
                   "matrix": k in matrix}
    TIMES["load_s"] = round(time.time() - T0, 1)

    # ---------- coverage ----------
    cov_rows = []
    all_label = f"ALL_{tgt_label}"
    for sp in in_sp + [all_label]:
        sel = [p for p in prot.values() if sp == all_label or p["sp"] == sp]
        for subset, rows in (("universe", sel), ("candidates", [p for p in sel if p["cand"]])):
            n = len(rows)
            r = {"species": sp, "subset": subset, "n_proteins": n,
                 "n_in_matrix": sum(1 for p in rows if p["matrix"])}
            for v in VOCABS:
                kk = sum(1 for p in rows if p[v])
                r[f"n_{v}"] = kk
                r[f"frac_{v}"] = round(kk / n, 4) if n else 0
            kk = sum(1 for p in rows if p["go"] or p["pfam"] or p["interpro"])
            r["n_any"] = kk
            r["frac_any"] = round(kk / n, 4) if n else 0
            cov_rows.append(r)
    with open(out / "coverage.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cov_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(cov_rows)
    # Per-vocabulary species coverage (also in PROVENANCE.yaml). The per-species counts
    # are the n_<vocab> columns of coverage.tsv; this names the species without terms.
    vocab_cov = vocabulary_coverage(cov_rows)
    for v, c in vocab_cov.items():
        if c["species_without_terms"]:
            log(f"NOTE: {v} terms in {c['n_with']} of {c['n_species']} {tgt_label} species; "
                f"none in: {', '.join(c['species_without_terms'])}")
    # Annotation source per species: the sidecar <file>.provenance.yaml of a non-UniProt
    # annotation file (bin/hmm_domtbl_to_annotations.py) names it; no sidecar = not recorded.
    ann_sources = {}
    for sp, f in ann_files.items():
        side = Path(f + ".provenance.yaml")
        src = "not recorded (no sidecar)"
        if side.exists():
            try:
                import yaml
                rec = yaml.safe_load(side.read_text()) or {}
                src = rec.get("annotation_source") or src
                inputs[str(side)] = sha256(side)
            except Exception as e:  # noqa: BLE001
                src = f"sidecar unreadable ({e})"
        ann_sources[sp] = src

    # ---------- background clusters ----------
    bg_rep, clust_meta = {}, None
    if args.cluster_background:
        bg = sorted((k, seqs[k]) for k in seqs if k not in candidates)
        t = time.time()
        bg_rep, clust_meta = cluster_background(bg, cache_dir, args.mmseqs, args.threads)
        TIMES["background_cluster_s"] = round(time.time() - t, 1)
        log(f"background clusters: {clust_meta['n_clusters']} clusters from {clust_meta['n_sequences']} proteins")

    # ---------- units ----------
    def build_units(sp_filter, level):
        groups = defaultdict(list)
        for key, p in prot.items():
            if sp_filter and p["sp"] != sp_filter:
                continue
            if level == "protein":
                groups[(p["cand"], key)].append(key)
            elif p["cand"]:
                groups[(True, "fam:" + family.get(key, f"singleton:{key}"))].append(key)
            elif args.cluster_background:
                groups[(False, "bgc:" + bg_rep[key])].append(key)
            else:
                groups[(False, key)].append(key)
        units = []
        for (is_c, uid), members in groups.items():
            terms = {v: set().union(*(prot[m][v] for m in members)) for v in DOMAIN_VOCABS}
            go = defaultdict(set)
            for m in members:
                for t, evs in prot[m]["go"].items():
                    go[t] |= evs
            terms["go"] = go
            units.append((uid, is_c, members, terms))
        return units

    scopes = [("whole_set", None)] + [(sp, sp) for sp in in_sp]
    rows_by_vocab = {"domains": [], "go": []}
    test_rows = []
    t = time.time()
    for scope, spf in scopes:
        for level in ("protein", "family"):
            units = build_units(spf, level)
            n_c_total = sum(1 for u in units if u[1])
            for v in VOCABS:
                if v == "pfam_clan" and not clans:
                    continue
                uni = [u for u in units if u[3][v]]
                M = len(uni)
                N = sum(1 for u in uni if u[1])
                bgc, kc, cid, non_iea = defaultdict(int), defaultdict(int), defaultdict(list), defaultdict(int)
                evc = defaultdict(lambda: defaultdict(int))  # GO term -> evidence code -> candidate units
                for u in uni:
                    for term in u[3][v]:
                        bgc[term] += 1
                        if u[1]:
                            kc[term] += 1
                            cid[term].extend(u[2])
                            if v == "go":
                                if u[3][v][term] - {"IEA"}:
                                    non_iea[term] += 1
                                for ev in u[3][v][term]:
                                    evc[term][ev or "NA"] += 1
                rows = []
                if N > 0:
                    for term, n in bgc.items():
                        if n < args.min_bg:
                            continue
                        k = kc.get(term, 0)
                        p = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
                        exp = n * N / M
                        r = {"scope": scope, "level": level, "vocab": v, "term": term,
                             "name": names[v].get(term, ""),
                             "k_cand_with_term": k, "N_cand_annotated": N,
                             "n_bg_with_term": n, "M_bg_annotated": M,
                             "expected": round(exp, 3), "fold": round(k / exp, 3) if exp else "",
                             "pvalue": p}
                        if v == "go":
                            r["n_cand_units_non_IEA"] = non_iea.get(term, 0)
                            r["go_evidence"] = ("none" if k == 0 else
                                                "IEA_only" if non_iea.get(term, 0) == 0 else "has_non_IEA")
                            r["cand_evidence_codes"] = ";".join(
                                f"{e}:{c}" for e, c in sorted(evc[term].items(), key=lambda x: (-x[1], x[0])))
                        rows.append(r)
                    for r, q in zip(rows, bh_fdr([r["pvalue"] for r in rows])):
                        r["fdr_bh"] = q
                        r["candidate_ids"] = ";".join(sorted(cid[r["term"]])) if q < args.report_ids_fdr else ""
                n_sig = sum(1 for r in rows if r["fdr_bh"] < 0.05)
                n_sig_k3 = sum(1 for r in rows if r["fdr_bh"] < 0.05 and r["k_cand_with_term"] >= 3)
                test_rows.append({
                    "study": args.study, "universe": args.universe,
                    "background_clustered": bool(args.cluster_background),
                    "scope": scope, "level": level, "vocab": v,
                    "role": "secondary" if v == "go" else "primary",
                    "n_cand_units": n_c_total, "N_cand_annotated": N,
                    "frac_cand_annotated": round(N / n_c_total, 4) if n_c_total else 0,
                    "n_bg_units": len(units) - n_c_total, "M_units_annotated": M,
                    "n_terms_tested": len(rows), "n_fdr_lt_0.05": n_sig, "n_fdr_lt_0.05_k_ge_3": n_sig_k3,
                    "status": "ok" if N > 0 else "skipped_no_annotated_candidates"})
                rows_by_vocab["go" if v == "go" else "domains"].extend(rows)
                log(f"{scope}\t{level}\t{v}\tN={N}/{n_c_total}\tM={M}\ttested={len(rows)}\tsig={n_sig}")
    TIMES["tests_s"] = round(time.time() - t, 1)

    for role, rows in rows_by_vocab.items():
        grp = defaultdict(list)
        for i, r in enumerate(rows):
            if r["scope"] != "whole_set":
                grp[(r["level"], r["vocab"])].append(i)
        for idx in grp.values():
            for i, q in zip(idx, bh_fdr([rows[i]["pvalue"] for i in idx])):
                rows[i]["fdr_bh_all_species"] = q
        for r in rows:
            r.setdefault("fdr_bh_all_species", "")
        rows.sort(key=lambda r: (r["scope"] != "whole_set", r["scope"], r["level"], r["vocab"], r["pvalue"], r["term"]))  # term: deterministic tie order
        fields = ["scope", "level", "vocab", "term", "name", "k_cand_with_term", "N_cand_annotated",
                  "n_bg_with_term", "M_bg_annotated", "expected", "fold", "pvalue", "fdr_bh",
                  "fdr_bh_all_species"]
        if role == "go":
            fields += ["go_evidence", "n_cand_units_non_IEA", "cand_evidence_codes"]
        fields += ["candidate_ids"]
        with gzip.open(out / f"enrichment_{role}.tsv.gz", "wt", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    with open(out / "tests.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(test_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(test_rows)

    # ---------- controls ----------
    if args.controls:
        inputs[str(args.controls)] = sha256(args.controls)
        ctrl = list(csv.DictReader(open(args.controls), delimiter="\t"))
        crow = []
        for c in ctrl:
            pf = set(split_ids(c.get("pfam_ids", "")))
            rx = re.compile(c["desc_regex"], re.I) if c.get("desc_regex") else None

            def hit(a):
                if not a:
                    return False
                if pf & set(split_ids(a["pfam"])):
                    return True
                return bool(rx and rx.search(a.get("desc", "") or ""))
            r = {"control": c["control"], "pfam_ids": c.get("pfam_ids", ""), "desc_regex": c.get("desc_regex", "")}
            # Column names follow the target group: IN_* / near_miss_out* (novelty),
            # OUT_* / near_miss_in* (loss).
            nm = f"near_miss_{oth_label.lower()}"
            in_hits = [k for k in universe if hit(ann_raw.get(k))]
            r[f"{tgt_label}_universe"] = len(in_hits)
            r[f"{tgt_label}_in_matrix"] = sum(1 for k in in_hits if k in matrix)
            r[f"{tgt_label}_pass_{'ingroup' if not loss else 'outgroup'}_rule"] = sum(1 for k in in_hits if k in in_rule)
            r["candidates"] = sum(1 for k in in_hits if k in candidates)
            r[f"{nm}1"] = sum(1 for k in in_hits if k in in_rule and matrix[k][1] == 1)
            r[f"{nm}2"] = sum(1 for k in in_hits if k in in_rule and matrix[k][1] == 2)
            for sp in in_sp + out_sp:
                if sp in all_sp_ann:
                    r[f"n_{sp}"] = sum(1 for a in all_sp_ann[sp].values() if hit(a))
                elif sp in in_sp:
                    r[f"n_{sp}"] = sum(1 for k in in_hits if sp_of[k] == sp)
                else:
                    r[f"n_{sp}"] = "NA"
            r["candidate_ids"] = ";".join(sorted(k for k in in_hits if k in candidates))
            r[f"{nm}1_ids"] = ";".join(sorted(k for k in in_hits if k in in_rule and matrix[k][1] == 1))
            crow.append(r)
        with open(out / "controls.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(crow[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(crow)

    # ---------- provenance ----------
    import goatools
    import scipy
    TIMES["total_s"] = round(time.time() - T0, 1)
    prov = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "study": args.study,
        "script": {"path": str(Path(sys.argv[0]).resolve()), "sha256": sha256(sys.argv[0])},
        "command": " ".join([sys.executable] + sys.argv),
        "class": "ephemeral analysis output (not tracked; regenerate with the command above)",
        "direction": args.direction,
        "target_group": {"label": tgt_label, "species": in_sp},
        "other_group": {"label": oth_label, "species": out_sp},
        "universe": args.universe, "annotations": args.annotations,
        "annotation_files": ann_files,
        "annotation_sources": ann_sources,
        "vocabulary_coverage": vocab_cov,
        "parameters": {"min_bg": args.min_bg, "exclude_iea": args.exclude_iea,
                       "report_ids_fdr": args.report_ids_fdr, "ingroup_min_frac": args.ingroup_min_frac,
                       "other_max_frac": args.other_max_frac,
                       "cluster_background": args.cluster_background},
        "statistics": "one-sided hypergeometric scipy.stats.hypergeom.sf(k-1,M,n,N); BH (lib/stats_util.bh_fdr) "
                      "within scope x level x vocabulary; scope-B fdr_bh_all_species = BH across species tests "
                      "of one level x vocabulary. goatools only parses the OBO and gives is_a ancestors.",
        "candidate_check": cand_check,
        "background_clustering": clust_meta,
        "pipeline_candidate_families": (f"{fam_source} (nf_NovInvenio MMSEQS_CLUSTER: mmseqs easy-cluster "
                                        "--min-seq-id 0.3 -c 0.8 --cov-mode 0 on "
                                        f"{'loss_candidates.fa' if loss else 'candidates.fa'})"),
        "candidates_without_family_row": len(nofam),
        "go_obo": {"path": str(args.obo), "data_version": obo_version, "sha256": sha256(args.obo)},
        "pfam_dat": ({"path": str(args.pfam_dat), "sha256": sha256(args.pfam_dat),
                      "release_note": "release taken from directory name; the .dat file has no version line"}
                     if args.pfam_dat else None),
        "interpro_entries": str(args.interpro_entries) if args.interpro_entries else None,
        "inputs_sha256": inputs,
        "tools": {"python": platform.python_version(), "scipy": scipy.__version__,
                  "goatools": getattr(goatools, "__version__", "unknown"),
                  "mmseqs": tool_version([args.mmseqs, "version"]) if args.cluster_background else None},
        "runtimes_s": TIMES,
        "host": platform.node(), "slurm_job": os.environ.get("SLURM_JOB_ID"),
    }
    try:
        import yaml
        with open(out / "PROVENANCE.yaml", "w") as fh:
            yaml.safe_dump(prov, fh, sort_keys=False, width=120)
    except ImportError:
        json.dump(prov, open(out / "PROVENANCE.yaml", "w"), indent=1)  # JSON is valid YAML
    log(f"done; outputs in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
