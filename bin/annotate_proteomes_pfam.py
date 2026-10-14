#!/usr/bin/env python3
"""Pfam-annotate the proteomes of a study that have no annotation file (NCBI/JGI sources).

Why: bin/build_study_config.py writes <study>/annotations/<Protein_Accession>.tsv.gz only
for UniProt proteomes (bin/extract_dat_annotations.py). Proteomes from NCBI or JGI get no
file, so the enrichment analysis has no terms for them. This script runs a Pfam-A
hmmsearch (HMMER 3, --cut_ga) on each of these proteomes on SLURM and converts the result
with bin/hmm_domtbl_to_annotations.py into the same file format. These files have Pfam
only: GO and InterPro stay empty.

hmmsearch, not hmmscan: on 500 PiroE2sp proteins (4 cores, job 29943954) hmmsearch took
81 s and hmmscan 362 s; both gave the same 525 (protein, Pfam accession) pairs.

Selection: every species.csv row with Protein_Source != "uniprot" whose
annotations/<Protein_Accession>.tsv.gz does not exist (--force: also when it exists).
UniProt rows are never selected; their file comes from bin/build_study_config.py.
--species limits the selection to the given Short names.

Job sizing: measured hmmsearch --cpu 8 --cut_ga wall-clock time is about 0.1 s per
protein (20,214 proteins 1,911 s; 23,663 proteins 2,448 s; 27,675 proteins 2,765 s;
milan nodes). One job runs at most --max-proteins-per-job proteins (default 35,000,
about 1 h). A larger proteome is split by sequence count into equal contiguous chunks;
the finalize step merges the chunk tables.

Modes (one of):
  (default)   write the job scripts and plan.json; print the sbatch commands
  --dry-run   print the plan; write nothing
  --submit    write the scripts, submit the chunk jobs, then submit one finalize job per
              species with --dependency=afterok on its chunk jobs
  --local     write the scripts and run them here, one after the other (needs $SCRATCH,
              hmmsearch and zstd; use inside an interactive SLURM job only)
  --finalize  merge the chunk tables of each selected species and write
              annotations/<Protein_Accession>.tsv.gz + .provenance.yaml

Files, per species, under --work-dir (default <study>/analysis/pfam_hmmsearch/):
  <Short>/plan.json             proteome path, protein count, chunks, Pfam path, job ids
  <Short>/chunk_NNN.sbatch      one hmmsearch job per chunk
  <Short>/finalize.sbatch       the dependent finalize job
  <Short>/chunk_NNN.domtbl.zst  hmmsearch --domtblout per chunk (compressed)
  <Short>/logs/                 SLURM logs
  <Short>.pfam.hmmsearch.domtbl.zst   merged domain table (finalize)

Job scripts use absolute paths only, --constraint=milan, and node-local $SCRATCH for the
chunk FASTA and the raw domain table.

Usage:
  bin/annotate_proteomes_pfam.py --study fungi/neocallimastigomycota_vs_early_fungi --dry-run
  bin/annotate_proteomes_pfam.py --study fungi/x --submit
  bin/annotate_proteomes_pfam.py --study fungi/x --species Neocalif --finalize
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))
sys.path.insert(0, str(REPO_ROOT / "lib"))

DEFAULT_PFAM_HMM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
DEFAULT_PYTHON = str(REPO_ROOT / ".pixi" / "envs" / "default" / "bin" / "python3")
PIXI_BIN = str(REPO_ROOT / ".pixi" / "envs" / "default" / "bin")
SEC_PER_PROTEIN = 0.1      # measured, hmmsearch --cpu 8 --cut_ga on milan nodes (see docstring)
MAX_PROTEINS_PER_JOB = 35000


@dataclass
class Species:
    short: str
    protein_accession: str
    protein_source: str
    taxon_id: str
    pep: Path
    annotation: Path
    n_proteins: int | None = None
    chunks: list[tuple[int, int]] = field(default_factory=list)


# --- pure logic (tested) ----------------------------------------------------------


def select_species(rows: list[dict], annotations_dir: Path, pep_of: dict[str, Path],
                   only: list[str] | None = None, force: bool = False
                   ) -> tuple[list[Species], list[tuple[str, str]]]:
    """Return (selected, skipped). skipped = [(Short, reason)]."""
    selected, skipped = [], []
    known = {r["Short"] for r in rows}
    for name in only or []:
        if name not in known:
            raise SystemExit(f"ERROR: --species {name} is not in species.csv")
    for r in rows:
        short = r["Short"]
        if only and short not in only:
            continue
        acc = (r.get("Protein_Accession") or "").strip()
        src = (r.get("Protein_Source") or "").strip().lower()
        ann = annotations_dir / f"{acc}.tsv.gz"
        if src == "uniprot":
            skipped.append((short, "UniProt source: annotation comes from bin/build_study_config.py"))
            continue
        if not acc:
            skipped.append((short, "no Protein_Accession in species.csv"))
            continue
        if ann.exists() and not force:
            skipped.append((short, f"{ann.name} exists (use --force to redo)"))
            continue
        selected.append(Species(short, acc, src, (r.get("Taxon_ID") or "").strip(),
                                pep_of.get(short, Path()), ann))
    return selected, skipped


def plan_chunks(n_proteins: int, max_per_job: int = MAX_PROTEINS_PER_JOB) -> list[tuple[int, int]]:
    """Split n proteins into the fewest equal contiguous chunks of <= max_per_job.
    Return [(start, end)] as 0-based, end-exclusive protein (record) indices."""
    if n_proteins <= 0:
        return []
    k = math.ceil(n_proteins / max_per_job)
    base, extra = divmod(n_proteins, k)
    out, start = [], 0
    for i in range(k):
        size = base + (1 if i < extra else 0)
        out.append((start, start + size))
        start += size
    return out


def time_limit(n_proteins: int, sec_per_protein: float = SEC_PER_PROTEIN) -> str:
    """SLURM -t for one chunk: 2x the estimate plus 1 h, at least 2 h. HH:MM:SS."""
    est = n_proteins * sec_per_protein
    hours = max(2, math.ceil(2 * est / 3600) + 1)
    return f"{hours:02d}:00:00"


def chunk_script(sp: Species, idx: int, start: int, end: int, work: Path, pfam_hmm: str,
                 cpus: int = 8, partition: str = "stajichlab", constraint: str = "milan",
                 sec_per_protein: float = SEC_PER_PROTEIN) -> str:
    n = end - start
    name = f"pfam_{sp.short}_c{idx:03d}"
    out = work / sp.short / f"chunk_{idx:03d}.domtbl.zst"
    q = shlex.quote
    return f"""#!/usr/bin/bash
#SBATCH -p {partition} -N 1 -n {cpus} --mem 16G -t {time_limit(n, sec_per_protein)} --constraint={constraint}
#SBATCH -J {name}
#SBATCH -o {work / sp.short / 'logs' / (name + '.%j.out')}
# Written by bin/annotate_proteomes_pfam.py. Pfam-A hmmsearch (--cut_ga) of proteins
# {start + 1}-{end} ({n} proteins) of {sp.short}. Estimated {n * sec_per_protein / 60:.0f} min at
# {sec_per_protein} s/protein. --constraint={constraint}: AVX2 nodes only.
# Absolute paths only (no BASH_SOURCE). Raw files stay in node-local $SCRATCH.
set -euo pipefail
export PATH={q(PIXI_BIN)}:"$PATH"
W="${{SCRATCH:?SCRATCH must be set (run under SLURM)}}/{name}"
mkdir -p "$W"
awk -v s={start} -v e={end} '/^>/{{n++}} n>s && n<=e' {q(str(sp.pep))} > "$W/in.fa"
got=$(grep -c '^>' "$W/in.fa" || true)
echo "proteins: $got (expected {n})"
if [ "$got" -ne {n} ]; then echo "ERROR: chunk has $got proteins, expected {n}" >&2; exit 1; fi
hmmsearch -h | sed -n 2p
t0=$(date +%s)
hmmsearch --cpu "${{SLURM_CPUS_ON_NODE:-{cpus}}}" --cut_ga --noali -o /dev/null \\
    --domtblout "$W/out.domtbl" {q(pfam_hmm)} "$W/in.fa"
echo "hmmsearch seconds: $(( $(date +%s) - t0 ))"
tail -n 1 "$W/out.domtbl" | grep -q '^# \\[ok\\]' || {{ echo "ERROR: domain table has no [ok] line" >&2; exit 1; }}
zstd -T0 -q -f "$W/out.domtbl" -o {q(str(out) + '.tmp')}
mv -f {q(str(out) + '.tmp')} {q(str(out))}
rm -rf "$W"
echo "wrote {out}"
"""


def finalize_script(sp: Species, study: str, work: Path, annotations_dir: Path, pfam_hmm: str,
                    force: bool, python: str = DEFAULT_PYTHON, partition: str = "stajichlab",
                    constraint: str = "milan") -> str:
    name = f"pfam_{sp.short}_final"
    q = shlex.quote
    cmd = [python, "-I", str(REPO_ROOT / "bin" / "annotate_proteomes_pfam.py"), "--study", study,
           "--species", sp.short, "--finalize", "--work-dir", str(work),
           "--annotations-dir", str(annotations_dir), "--pfam-hmm", pfam_hmm]
    if force:
        cmd.append("--force")
    return f"""#!/usr/bin/bash
#SBATCH -p {partition} -N 1 -n 1 --mem 4G -t 00:30:00 --constraint={constraint}
#SBATCH -J {name}
#SBATCH -o {work / sp.short / 'logs' / (name + '.%j.out')}
# Written by bin/annotate_proteomes_pfam.py. Merges the chunk domain tables of {sp.short} and
# writes {sp.annotation} (+ .provenance.yaml). Submitted with --dependency=afterok on the chunk jobs.
set -euo pipefail
export PATH={q(PIXI_BIN)}:"$PATH"
{' '.join(q(c) for c in cmd)}
"""


def merge_domtbls(texts: list[str]) -> str:
    """Merge chunk domain tables: header lines of the first table, all data lines in chunk
    order, then the trailer of the first table. Each table must end with '# [ok]'."""
    head, data = [], []
    for i, t in enumerate(texts):
        lines = t.splitlines()
        if not lines or lines[-1].strip() != "# [ok]":
            raise SystemExit(f"ERROR: chunk {i} domain table does not end with '# [ok]' (incomplete run?)")
        first_data = next((j for j, ln in enumerate(lines) if not ln.startswith("#")), None)
        rows = [ln for ln in lines if ln.strip() and not ln.startswith("#")]
        data.extend(rows)
        if i == 0:
            if first_data is None:  # no hits in the first chunk: header = the 3 column lines
                head = [ln for ln in lines[:3] if ln.startswith("#")]
                trailer = lines[3:]
            else:
                head = lines[:first_data]
                last_data = max(j for j, ln in enumerate(lines) if ln.strip() and not ln.startswith("#"))
                trailer = lines[last_data + 1:]
    note = [f"# Merged by bin/annotate_proteomes_pfam.py from {len(texts)} chunk table(s); "
            "trailer below is from chunk 0."] if len(texts) > 1 else []
    return "\n".join(head + data + note + trailer) + "\n"


# --- I/O helpers ------------------------------------------------------------------


def count_fasta(path: Path) -> int:
    n = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def pep_paths(study_dir: Path, rows: list[dict]) -> dict[str, Path]:
    """Short -> protein FASTA. config.csv Protein column if present, else data_dir/pep/<Short>.pep.fa."""
    cfg = study_dir / "config.csv"
    by_cfg = {}
    if cfg.is_file():
        by_cfg = {r["Short"]: r.get("Protein", "") for r in read_csv(cfg)}
    out = {}
    for r in rows:
        name = by_cfg.get(r["Short"]) or f"{r['Short']}.pep.fa"
        out[r["Short"]] = study_dir / "data_dir" / "pep" / name
    return out


def zstd_read(path: Path) -> str:
    res = subprocess.run([shutil.which("zstd") or "zstd", "-dc", str(path)], capture_output=True, check=True)
    return res.stdout.decode()


def zstd_write(text: str, path: Path) -> None:
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "wb") as fh:
        subprocess.run([shutil.which("zstd") or "zstd", "-q", "-T0", "-c"], input=text.encode(),
                       stdout=fh, check=True)
    tmp.replace(path)


def finalize(sp: Species, work: Path, pfam_hmm: str, force: bool) -> dict:
    d = work / sp.short
    plan_path = d / "plan.json"
    if not plan_path.is_file():
        raise SystemExit(f"ERROR: {plan_path} not found; write the jobs first")
    plan = json.loads(plan_path.read_text())
    if sp.annotation.exists() and not force:
        raise SystemExit(f"ERROR: {sp.annotation} exists; use --force to overwrite")
    n_now = count_fasta(sp.pep)
    if n_now != plan["n_proteins"]:
        raise SystemExit(f"ERROR: {sp.pep} has {n_now} proteins; plan.json says {plan['n_proteins']}")
    chunk_files = [d / f"chunk_{i:03d}.domtbl.zst" for i in range(len(plan["chunks"]))]
    missing = [str(p) for p in chunk_files if not p.is_file()]
    if missing:
        raise SystemExit(f"ERROR: missing chunk tables: {missing}")
    merged = work / f"{sp.short}.pfam.hmmsearch.domtbl.zst"
    zstd_write(merge_domtbls([zstd_read(p) for p in chunk_files]), merged)
    import hmm_domtbl_to_annotations as h2a
    argv = ["--domtbl", str(merged), "--out", str(sp.annotation), "--fasta", str(sp.pep),
            "--taxon-id", sp.taxon_id, "--pfam-hmm", pfam_hmm]
    summary = h2a.convert([merged], sp.annotation, fasta=sp.pep, taxon_id=sp.taxon_id, pfam_hmm=pfam_hmm,
                          argv=[str(REPO_ROOT / "bin" / "hmm_domtbl_to_annotations.py")] + argv)
    plan["finalized"] = {"date": datetime.datetime.now().isoformat(timespec="seconds"),
                         "merged_domtbl": str(merged), "annotation": str(sp.annotation), "counts": summary}
    plan_path.write_text(json.dumps(plan, indent=1) + "\n")
    return summary


# --- main -------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="<domain>/<study> under studies/")
    ap.add_argument("--species", nargs="+", help="Short names (default: all eligible species)")
    ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM_HMM, help=f"Pfam-A.hmm (default {DEFAULT_PFAM_HMM})")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="print the plan only")
    mode.add_argument("--submit", action="store_true", help="submit chunk jobs + dependent finalize jobs")
    mode.add_argument("--local", action="store_true", help="run the job scripts here, in order")
    mode.add_argument("--finalize", action="store_true", help="merge chunk tables and write annotations")
    ap.add_argument("--force", action="store_true", help="redo species that already have an annotation file")
    ap.add_argument("--max-proteins-per-job", type=int, default=MAX_PROTEINS_PER_JOB)
    ap.add_argument("--sec-per-protein", type=float, default=SEC_PER_PROTEIN,
                    help="wall-clock seconds per protein at --cpus 8 (job-time estimate)")
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--partition", default="stajichlab")
    ap.add_argument("--constraint", default="milan")
    ap.add_argument("--work-dir", type=Path, help="default <study>/analysis/pfam_hmmsearch")
    ap.add_argument("--annotations-dir", type=Path, help="default <study>/annotations")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    study_dir = (args.repo_root / "studies" / args.study).resolve()
    if not (study_dir / "species.csv").is_file():
        raise SystemExit(f"ERROR: {study_dir}/species.csv not found")
    work = (args.work_dir or study_dir / "analysis" / "pfam_hmmsearch").resolve()
    ann_dir = (args.annotations_dir or study_dir / "annotations").resolve()
    rows = read_csv(study_dir / "species.csv")
    # --finalize runs after the jobs; the annotation file does not exist yet, so the selection
    # itself must not skip on it (finalize() checks --force before it writes).
    selected, skipped = select_species(rows, ann_dir, pep_paths(study_dir, rows), args.species,
                                       force=args.force or args.finalize)
    for short, why in skipped:
        print(f"skip {short}: {why}", file=sys.stderr)
    if not selected:
        print("nothing to do", file=sys.stderr)
        return 0

    if args.finalize:
        for sp in selected:
            if sp.annotation.exists() and not args.force:
                print(f"skip {sp.short}: {sp.annotation.name} exists (use --force to redo)", file=sys.stderr)
                continue
            s = finalize(sp, work, args.pfam_hmm, args.force)
            pct = 100 * s["n_proteins_with_pfam"] / s["n_rows"] if s["n_rows"] else 0
            print(f"{sp.short}: wrote {sp.annotation} ({s['n_rows']} proteins, "
                  f"{s['n_proteins_with_pfam']} with >= 1 Pfam, {pct:.1f}%)", file=sys.stderr)
        return 0

    for sp in selected:
        if not sp.pep.is_file():
            raise SystemExit(f"ERROR: {sp.short}: protein FASTA {sp.pep} not found "
                             "(run bin/build_study_config.py first)")
        sp.n_proteins = count_fasta(sp.pep)
        sp.chunks = plan_chunks(sp.n_proteins, args.max_proteins_per_job)
        est = max((e - s for s, e in sp.chunks), default=0) * args.sec_per_protein / 60
        print(f"{sp.short}: {sp.n_proteins} proteins -> {len(sp.chunks)} job(s), "
              f"largest about {est:.0f} min; output {sp.annotation}", file=sys.stderr)
    if args.dry_run:
        return 0

    for sp in selected:
        d = work / sp.short
        (d / "logs").mkdir(parents=True, exist_ok=True)
        scripts = []
        for i, (s, e) in enumerate(sp.chunks):
            p = d / f"chunk_{i:03d}.sbatch"
            p.write_text(chunk_script(sp, i, s, e, work, args.pfam_hmm, args.cpus, args.partition,
                                      args.constraint, args.sec_per_protein))
            scripts.append(p)
        fin = d / "finalize.sbatch"
        fin.write_text(finalize_script(sp, args.study, work, ann_dir, args.pfam_hmm, args.force,
                                       partition=args.partition, constraint=args.constraint))
        plan = {"study": args.study, "short": sp.short, "protein_accession": sp.protein_accession,
                "protein_source": sp.protein_source, "pep": str(sp.pep), "n_proteins": sp.n_proteins,
                "chunks": sp.chunks, "pfam_hmm": args.pfam_hmm, "max_proteins_per_job": args.max_proteins_per_job,
                "sec_per_protein": args.sec_per_protein, "cpus": args.cpus, "annotation": str(sp.annotation),
                "created": datetime.datetime.now().isoformat(timespec="seconds"),
                "chunk_scripts": [str(p) for p in scripts], "finalize_script": str(fin), "job_ids": []}
        if args.submit:
            ids = []
            for p in scripts:
                res = subprocess.run(["sbatch", "--parsable", str(p)], capture_output=True, text=True, check=True)
                ids.append(res.stdout.strip().split(";")[0])
            res = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{':'.join(ids)}", str(fin)],
                                 capture_output=True, text=True, check=True)
            plan["job_ids"] = ids
            plan["finalize_job_id"] = res.stdout.strip().split(";")[0]
            print(f"{sp.short}: submitted chunk jobs {ids}, finalize {plan['finalize_job_id']}", file=sys.stderr)
        (d / "plan.json").write_text(json.dumps(plan, indent=1) + "\n")
        if args.local:
            if not os.environ.get("SCRATCH"):
                raise SystemExit("ERROR: --local needs $SCRATCH (run inside a SLURM job)")
            for p in scripts:
                subprocess.run(["bash", str(p)], check=True)
            finalize(sp, work, args.pfam_hmm, args.force)
        elif not args.submit:
            for p in scripts:
                print(f"  sbatch {p}", file=sys.stderr)
            print(f"  then: sbatch --dependency=afterok:<chunk job ids> {fin}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
