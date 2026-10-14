#!/usr/bin/env python3
"""Pfam-annotate the proteomes of a study that have no annotation file (NCBI/JGI/local).

Why: bin/build_study_config.py writes <study>/annotations/<name>.tsv.gz only for UniProt
proteomes (bin/extract_dat_annotations.py). Proteomes from NCBI, JGI or a local FASTA
get no file, so the enrichment analysis has no terms for them. This script runs a Pfam-A
hmmsearch (HMMER 3, --cut_ga) on each of these proteomes on SLURM and converts the
result with bin/hmm_domtbl_to_annotations.py into the same file format. These files
have Pfam only: GO and InterPro stay empty.

hmmsearch, not hmmscan: on 500 PiroE2sp proteins (4 cores, job 29943954) hmmsearch took
81 s and hmmscan 362 s; both gave the same 525 (protein, Pfam accession) pairs.

File name: <name> comes from lib/annotation_names.py (annotation_names): the
Protein_Accession for UniProt/NCBI accessions (unchanged: UP000193719, GCA_002104975.1);
for a local FASTA path, its basename without the proteome/compression suffix; Short as
the fallback. Hard guard: this script refuses to write any file outside the study's
annotations dir, the work dir, or <study>/annotations_orig/ (ValueError).

Selection: every species.csv row with Protein_Source != "uniprot" whose
annotations/<name>.tsv.gz does not exist (--force: also when it exists). UniProt rows
are never selected, except with --fill-empty-uniprot (below). --species limits the
selection to the given Short names.

--fill-empty-uniprot (opt-in): also select UniProt rows whose annotation file exists
and has no GO, Pfam or InterPro value in any row (a new UniProt proteome whose entries
have no cross-references yet). The finalize step copies the untouched file to
<study>/annotations_orig/<name>.tsv.gz (outside annotations/, because
bin/sync_reports.sh merges every annotations/*.tsv.gz), then fills only the pfam_ids and
pfam_names columns of the rows whose accession has a Pfam hit. All other columns and
rows stay as they are. A .provenance.yaml sidecar records both sources. A file that has
any term is never changed. Note: bin/build_study_config.py rewrites UniProt annotation
files from the .dat.gz each time it runs; that replaces a filled file with the
original empty one, and the sidecar checksum then no longer matches.

Job sizing: measured hmmsearch --cpu 8 --cut_ga wall-clock time is about 0.1 s per
protein (20,214 proteins 1,911 s; 23,663 proteins 2,448 s; 27,675 proteins 2,765 s;
milan nodes). One hmmsearch run (a "chunk") has at most --max-proteins-per-job proteins
(default 35,000, about 1 h). A larger proteome is split by sequence count into equal
contiguous chunks; the finalize step merges the chunk tables. Chunks are then packed
into batch jobs (first-fit decreasing by protein count) so that one batch job runs
about --batch-minutes (default 75) of estimated hmmsearch time: several small
proteomes run one after the other in one job, each writing its own domain table. One
finalize job (--dependency=afterany on all batch jobs) converts every species and
writes the provenance sidecars; a species with a missing chunk table is reported and
the job exits non-zero after it has done the others.

Modes (one of):
  (default)   write the job scripts and plan.json files; print the sbatch commands
  --dry-run   print the plan; write nothing
  --submit    write the scripts, submit the batch jobs, then one finalize job
  --local     write the scripts and run them here, one after the other (needs $SCRATCH,
              hmmsearch and zstd; use inside an interactive SLURM job only)
  --finalize  merge the chunk tables of each selected species and write
              annotations/<name>.tsv.gz + .provenance.yaml

Files under --work-dir (default <study>/analysis/pfam_hmmsearch/):
  <Short>/plan.json             proteome path, protein count, chunks, batches, job ids
  <Short>/chunk_NNN.sbatch      one hmmsearch run per chunk (run by a batch job; can
                                also be submitted on its own)
  <Short>/chunk_NNN.domtbl.zst  hmmsearch --domtblout per chunk (compressed)
  <Short>.pfam.hmmsearch.domtbl.zst   merged domain table (finalize)
  batches/<run tag>/batch_NNN.sbatch  batch jobs; finalize.sbatch; logs/; batches.json

Job scripts use absolute paths only, --constraint=milan, and node-local $SCRATCH for the
chunk FASTA and the raw domain table.

Usage:
  bin/annotate_proteomes_pfam.py --study fungi/neocallimastigomycota_vs_early_fungi --dry-run
  bin/annotate_proteomes_pfam.py --study fungi/x --submit
  bin/annotate_proteomes_pfam.py --study fungi/x --species Neocalif --finalize
  bin/annotate_proteomes_pfam.py --study fungi/pezizo_set1 --fill-empty-uniprot --dry-run
"""
from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import io
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
from annotation_names import annotation_names, annotation_path, ensure_within, is_safe_name  # noqa: E402

DEFAULT_PFAM_HMM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
DEFAULT_PYTHON = str(REPO_ROOT / ".pixi" / "envs" / "default" / "bin" / "python3")
PIXI_BIN = str(REPO_ROOT / ".pixi" / "envs" / "default" / "bin")
SEC_PER_PROTEIN = 0.1      # measured, hmmsearch --cpu 8 --cut_ga on milan nodes (see docstring)
MAX_PROTEINS_PER_JOB = 35000
BATCH_MINUTES = 75
TERM_COLUMNS = ("go_ids", "pfam_ids", "interpro_ids")


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
    mode: str = "new"          # "new" (non-UniProt, new file) or "fill_uniprot"


# --- pure logic (tested) ----------------------------------------------------------


def has_any_terms(path: Path) -> bool:
    """True if any row of an annotation TSV (.gz) has a GO, Pfam or InterPro value."""
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if any((row.get(c) or "").strip() for c in TERM_COLUMNS):
                return True
    return False


def select_species(rows: list[dict], annotations_dir: Path, pep_of: dict[str, Path],
                   only: list[str] | None = None, force: bool = False,
                   fill_empty_uniprot: bool = False
                   ) -> tuple[list[Species], list[tuple[str, str]]]:
    """Return (selected, skipped). skipped = [(Short, reason)].
    Annotation file names come from lib/annotation_names.py over ALL rows."""
    selected, skipped = [], []
    known = {r["Short"] for r in rows}
    for name in only or []:
        if name not in known:
            raise SystemExit(f"ERROR: --species {name} is not in species.csv")
    names = annotation_names(rows)
    for r in rows:
        short = r["Short"]
        if only and short not in only:
            continue
        if not is_safe_name(short):
            raise SystemExit(f"ERROR: Short {short!r} is not a safe file name (used in work paths)")
        acc = (r.get("Protein_Accession") or "").strip()
        src = (r.get("Protein_Source") or "").strip().lower()
        ann = annotation_path(annotations_dir, names[short])
        taxon = (r.get("Taxon_ID") or "").strip()
        if src == "uniprot":
            if not fill_empty_uniprot:
                skipped.append((short, "UniProt source: annotation comes from bin/build_study_config.py"))
            elif not ann.is_file():
                skipped.append((short, f"UniProt source without {ann.name} (run bin/build_study_config.py)"))
            elif has_any_terms(ann):
                skipped.append((short, f"UniProt source: {ann.name} has GO/Pfam/InterPro terms; not changed"))
            else:
                selected.append(Species(short, acc, src, taxon, pep_of.get(short, Path()), ann,
                                        mode="fill_uniprot"))
            continue
        if not acc:
            skipped.append((short, "no Protein_Accession in species.csv"))
            continue
        if ann.exists() and not force:
            skipped.append((short, f"{ann.name} exists (use --force to redo)"))
            continue
        selected.append(Species(short, acc, src, taxon, pep_of.get(short, Path()), ann))
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


def pack_batches(tasks: list[tuple[str, int, int]], max_proteins: int) -> list[list[tuple[str, int, int]]]:
    """Pack chunk tasks (Short, chunk index, n proteins) into batches of <= max_proteins total
    (first-fit decreasing; ties by Short, chunk). A task larger than max_proteins gets a
    batch of its own. Batches and the tasks inside them are in a deterministic order."""
    order = sorted(tasks, key=lambda t: (-t[2], t[0], t[1]))
    batches: list[list[tuple[str, int, int]]] = []
    loads: list[int] = []
    for t in order:
        for i, load in enumerate(loads):
            if load + t[2] <= max_proteins:
                batches[i].append(t)
                loads[i] += t[2]
                break
        else:
            batches.append([t])
            loads.append(t[2])
    return batches


def batch_capacity(batch_minutes: float, sec_per_protein: float) -> int:
    """Proteins per batch job for the target wall-clock time."""
    return max(1, int(batch_minutes * 60 / sec_per_protein))


def time_limit(n_proteins: int, sec_per_protein: float = SEC_PER_PROTEIN) -> str:
    """SLURM -t for one job: 2x the estimate plus 1 h, at least 2 h. HH:MM:SS."""
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
# Normally run by a batch job (bash <this file>); can also be submitted on its own.
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


def batch_script(idx: int, items: list[tuple[str, int, int, Path]], batch_dir: Path, cpus: int = 8,
                 partition: str = "stajichlab", constraint: str = "milan",
                 sec_per_protein: float = SEC_PER_PROTEIN) -> str:
    """items: [(Short, chunk index, n proteins, chunk script path)]. The batch runs the chunk
    scripts one after the other; a failed chunk does not stop the others; the job exits 1
    if any chunk failed."""
    total = sum(n for _, _, n, _ in items)
    name = f"pfam_b{idx:03d}"
    q = shlex.quote
    lines = "\n".join(f"#   {s} chunk {c:03d}: {n} proteins" for s, c, n, _ in items)
    scripts = " \\\n    ".join(q(str(p)) for *_, p in items)
    return f"""#!/usr/bin/bash
#SBATCH -p {partition} -N 1 -n {cpus} --mem 16G -t {time_limit(total, sec_per_protein)} --constraint={constraint}
#SBATCH -J {name}
#SBATCH -o {batch_dir / 'logs' / (name + '.%j.out')}
# Written by bin/annotate_proteomes_pfam.py. Batch {idx}: {len(items)} hmmsearch run(s), {total} proteins,
# estimated {total * sec_per_protein / 60:.0f} min at {sec_per_protein} s/protein. Runs in order:
{lines}
# Each chunk script writes its own domain table. --constraint={constraint}: AVX2 nodes only.
# Absolute paths only (no BASH_SOURCE).
set -uo pipefail
: "${{SCRATCH:?SCRATCH must be set (run under SLURM)}}"
failed=()
for s in {scripts}; do
    echo "== $(date +%T) start $s"
    if bash "$s"; then echo "== $(date +%T) done $s"; else echo "== $(date +%T) FAILED $s" >&2; failed+=("$s"); fi
done
if [ "${{#failed[@]}}" -gt 0 ]; then echo "ERROR: ${{#failed[@]}} chunk(s) failed: ${{failed[*]}}" >&2; exit 1; fi
echo "batch {idx}: all {len(items)} chunk(s) done"
"""


def finalize_script(shorts: list[str], study: str, work: Path, annotations_dir: Path, pfam_hmm: str,
                    force: bool, batch_dir: Path, python: str = DEFAULT_PYTHON,
                    partition: str = "stajichlab", constraint: str = "milan",
                    fill_empty_uniprot: bool = False) -> str:
    name = "pfam_final"
    q = shlex.quote
    cmd = [python, "-I", str(REPO_ROOT / "bin" / "annotate_proteomes_pfam.py"), "--study", study,
           "--species", *shorts, "--finalize", "--work-dir", str(work),
           "--annotations-dir", str(annotations_dir), "--pfam-hmm", pfam_hmm]
    if force:
        cmd.append("--force")
    if fill_empty_uniprot:
        cmd.append("--fill-empty-uniprot")
    hours = max(1, math.ceil(len(shorts) * 60 / 3600) + 1)   # about 1 min per species, + 1 h
    return f"""#!/usr/bin/bash
#SBATCH -p {partition} -N 1 -n 1 --mem 8G -t {hours:02d}:00:00 --constraint={constraint}
#SBATCH -J {name}
#SBATCH -o {batch_dir / 'logs' / (name + '.%j.out')}
# Written by bin/annotate_proteomes_pfam.py. Merges the chunk domain tables of {len(shorts)} species and
# writes annotations/<name>.tsv.gz (+ .provenance.yaml) for each. Submitted with
# --dependency=afterany on the batch jobs: a species with a missing chunk table is reported,
# the others are converted, and the job then exits non-zero.
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


def fill_rows(rows: list[dict], pfam_of: dict[str, tuple[str, str]]) -> dict:
    """Fill pfam_ids/pfam_names of rows in place from pfam_of {accession: (ids, names)}.
    Only rows whose accession is in pfam_of with a non-empty ids value change. Return counts."""
    n_filled = 0
    for r in rows:
        ids, names = pfam_of.get(r["accession"], ("", ""))
        if ids:
            r["pfam_ids"], r["pfam_names"] = ids, names
            n_filled += 1
    return {"n_rows": len(rows), "n_rows_filled": n_filled}


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


class Guard:
    """Hard guard: every file this script writes must be under one of the allowed roots."""

    def __init__(self, roots: list[Path]):
        self.roots = [Path(r).resolve() for r in roots]

    def __call__(self, path: Path) -> Path:
        return ensure_within(path, self.roots)


def write_text_guarded(guard: Guard, path: Path, text: str) -> None:
    guard(path)
    path.write_text(text)


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict]]:
    with gzip.open(path, "rt", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
        return list(rdr.fieldnames or []), rows


def write_tsv_gz(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
        gz.write(buf.getvalue().encode())
    tmp.replace(path)


def fill_uniprot(sp: Species, merged: Path, orig_dir: Path, pfam_hmm: str, guard: Guard) -> dict:
    """--fill-empty-uniprot finalize: keep the original UniProt file in orig_dir, fill
    pfam_ids/pfam_names of annotations/<name>.tsv.gz, write a sidecar."""
    import hmm_domtbl_to_annotations as h2a
    from provenance import build_record, sha256_of, write_record
    from uniprot_ids import bare_accession

    fieldnames, rows = read_tsv_gz(sp.annotation)
    for c in ("accession", "pfam_ids", "pfam_names"):
        if c not in fieldnames:
            raise SystemExit(f"ERROR: {sp.annotation} has no {c} column")
    if has_any_terms(sp.annotation):
        raise SystemExit(f"ERROR: {sp.annotation} has GO/Pfam/InterPro terms; --fill-empty-uniprot does not change it")
    orig = guard(orig_dir / sp.annotation.name)
    cur_sha = sha256_of(sp.annotation)
    if orig.exists():
        if sha256_of(orig) != cur_sha:
            raise SystemExit(f"ERROR: {orig} exists and differs from {sp.annotation}; not overwritten")
    else:
        orig.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(sp.annotation, orig)

    lines = h2a.read_text(merged)
    trailer = h2a.parse_trailer(lines)
    layout = h2a.detect_layout(lines, trailer, "auto")
    domains = h2a.parse_domtbl(lines, layout)
    fasta = h2a.read_fasta_headers(sp.pep)
    built = h2a.build_rows(domains, fasta, "")
    pfam_of = {bare_accession(b["accession"]): (b["pfam_ids"], b["pfam_names"]) for b in built}
    row_accs = {r["accession"] for r in rows}
    hit_accs = {a for a, (ids, _) in pfam_of.items() if ids}
    counts = fill_rows(rows, pfam_of)
    counts.update({
        "n_fasta_proteins": len(fasta),
        "n_fasta_proteins_with_row": len(set(pfam_of) & row_accs),
        "n_fasta_proteins_with_pfam": len(hit_accs),
        "n_proteins_with_pfam_but_no_row": len(hit_accs - row_accs),
        "n_rows_not_in_fasta": len(row_accs - set(pfam_of)),
        "n_domain_rows": len(domains), "n_pfam_families": len({d["pfam"] for d in domains}),
    })
    if hit_accs and counts["n_rows_filled"] == 0:
        raise SystemExit(f"ERROR: {sp.short}: {len(hit_accs)} proteins have Pfam hits but no row of "
                         f"{sp.annotation} matches their accession; file not changed")
    guard(sp.annotation)
    write_tsv_gz(sp.annotation, fieldnames, rows)
    release = h2a.pfam_release_from_path(pfam_hmm)
    cmd = " ".join(shlex.quote(x) for x in [sys.executable] + sys.argv)
    record = build_record(
        source_url=f"file://{orig}; file://{merged.resolve()}",
        source_release=f"UniProt annotation as extracted by bin/extract_dat_annotations.py (see "
                       f"DATA_MANIFEST.yaml for the proteome release); {release}; "
                       f"HMMER {trailer.get('version', 'not recorded')} ({layout})",
        license="UniProt: CC BY 4.0; " + h2a.PFAM_LICENSE,
        local_path=sp.annotation,
        derived_by=f"bin/annotate_proteomes_pfam.py --fill-empty-uniprot: {cmd}",
        extra={
            "annotation_source": "UniProt reference proteome (no GO/Pfam/InterPro in the UniProt "
                                 "entries) + Pfam-A HMM search (hmmsearch) for pfam_ids/pfam_names",
            "vocabularies_available": ["pfam", "pfam_clan (derived from Pfam-A.hmm.dat by the reader)"],
            "vocabularies_not_available": ["go", "interpro"],
            "statement": "The UniProt file had no GO, Pfam or InterPro value in any row. Only "
                         "pfam_ids and pfam_names were filled, from a Pfam-A hmmsearch (--cut_ga); "
                         "all other columns and rows are the UniProt values. GO and InterPro stay empty.",
            "original_file": {"path": str(orig), "sha256": cur_sha},
            "source_domtbl": {"path": str(merged.resolve()), "sha256": sha256_of(merged),
                              "hmmer_command": trailer.get("option_settings", "not recorded")},
            "pfam_hmm": pfam_hmm, "pfam_release": release,
            "protein_fasta": {"path": str(sp.pep.resolve()), "sha256": sha256_of(sp.pep)},
            "counts": counts,
            "note": "bin/build_study_config.py rewrites this file from the .dat.gz when it runs; "
                    "the checksum above then no longer matches the file.",
        },
    )
    side = guard(Path(str(sp.annotation) + ".provenance.yaml"))
    write_record(record, side)
    return {"n_rows": counts["n_rows"], "n_proteins_with_pfam": counts["n_rows_filled"], **counts}


def finalize(sp: Species, work: Path, pfam_hmm: str, force: bool, guard: Guard | None = None,
             orig_dir: Path | None = None) -> dict:
    guard = guard or Guard([work, sp.annotation.parent])
    d = work / sp.short
    plan_path = d / "plan.json"
    if not plan_path.is_file():
        raise SystemExit(f"ERROR: {plan_path} not found; write the jobs first")
    plan = json.loads(plan_path.read_text())
    if sp.mode == "new" and sp.annotation.exists() and not force:
        raise SystemExit(f"ERROR: {sp.annotation} exists; use --force to overwrite")
    n_now = count_fasta(sp.pep)
    if n_now != plan["n_proteins"]:
        raise SystemExit(f"ERROR: {sp.pep} has {n_now} proteins; plan.json says {plan['n_proteins']}")
    chunk_files = [d / f"chunk_{i:03d}.domtbl.zst" for i in range(len(plan["chunks"]))]
    missing = [str(p) for p in chunk_files if not p.is_file()]
    if missing:
        raise SystemExit(f"ERROR: {sp.short}: missing chunk tables: {missing}")
    merged = guard(work / f"{sp.short}.pfam.hmmsearch.domtbl.zst")
    zstd_write(merge_domtbls([zstd_read(p) for p in chunk_files]), merged)
    if sp.mode == "fill_uniprot":
        summary = fill_uniprot(sp, merged, orig_dir or sp.annotation.parent.parent / "annotations_orig",
                               pfam_hmm, guard)
    else:
        import hmm_domtbl_to_annotations as h2a
        guard(sp.annotation)
        guard(Path(str(sp.annotation) + ".provenance.yaml"))
        argv = ["--domtbl", str(merged), "--out", str(sp.annotation), "--fasta", str(sp.pep),
                "--taxon-id", sp.taxon_id, "--pfam-hmm", pfam_hmm]
        summary = h2a.convert([merged], sp.annotation, fasta=sp.pep, taxon_id=sp.taxon_id, pfam_hmm=pfam_hmm,
                              argv=[str(REPO_ROOT / "bin" / "hmm_domtbl_to_annotations.py")] + argv)
    plan["finalized"] = {"date": datetime.datetime.now().isoformat(timespec="seconds"),
                         "mode": sp.mode, "merged_domtbl": str(merged), "annotation": str(sp.annotation),
                         "counts": summary}
    write_text_guarded(guard, plan_path, json.dumps(plan, indent=1) + "\n")
    return summary


# --- main -------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="<domain>/<study> under studies/")
    ap.add_argument("--species", nargs="+", help="Short names (default: all eligible species)")
    ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM_HMM, help=f"Pfam-A.hmm (default {DEFAULT_PFAM_HMM})")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="print the plan only")
    mode.add_argument("--submit", action="store_true", help="submit batch jobs + one dependent finalize job")
    mode.add_argument("--local", action="store_true", help="run the chunk scripts here, in order")
    mode.add_argument("--finalize", action="store_true", help="merge chunk tables and write annotations")
    ap.add_argument("--force", action="store_true", help="redo non-UniProt species that already have a file")
    ap.add_argument("--fill-empty-uniprot", action="store_true",
                    help="also fill pfam_ids/pfam_names of UniProt annotation files that have no "
                         "GO/Pfam/InterPro value in any row (original kept in <study>/annotations_orig/)")
    ap.add_argument("--max-proteins-per-job", type=int, default=MAX_PROTEINS_PER_JOB,
                    help="max proteins per hmmsearch run (chunk)")
    ap.add_argument("--batch-minutes", type=float, default=BATCH_MINUTES,
                    help="target estimated hmmsearch minutes per batch job (default 75)")
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
    orig_dir = study_dir / "annotations_orig"
    guard = Guard([ann_dir, work, orig_dir])
    rows = read_csv(study_dir / "species.csv")
    # --finalize runs after the jobs; the annotation file does not exist yet, so the selection
    # itself must not skip on it (finalize() checks --force before it writes).
    selected, skipped = select_species(rows, ann_dir, pep_paths(study_dir, rows), args.species,
                                       force=args.force or args.finalize,
                                       fill_empty_uniprot=args.fill_empty_uniprot)
    for short, why in skipped:
        print(f"skip {short}: {why}", file=sys.stderr)
    if not selected:
        print("nothing to do", file=sys.stderr)
        return 0
    for sp in selected:                      # hard guard before anything is planned or written
        guard(sp.annotation)

    if args.finalize:
        errors = []
        for sp in selected:
            if sp.mode == "new" and sp.annotation.exists() and not args.force:
                print(f"skip {sp.short}: {sp.annotation.name} exists (use --force to redo)", file=sys.stderr)
                continue
            try:
                s = finalize(sp, work, args.pfam_hmm, args.force, guard, orig_dir)
            except SystemExit as e:
                print(str(e), file=sys.stderr)
                errors.append(f"{sp.short}: {e}")
                continue
            pct = 100 * s["n_proteins_with_pfam"] / s["n_rows"] if s["n_rows"] else 0
            print(f"{sp.short}: wrote {sp.annotation} ({sp.mode}; {s['n_rows']} rows, "
                  f"{s['n_proteins_with_pfam']} with >= 1 Pfam, {pct:.1f}%)", file=sys.stderr)
        if errors:
            raise SystemExit(f"ERROR: finalize failed for {len(errors)} species:\n" + "\n".join(errors))
        return 0

    tasks = []
    for sp in selected:
        if not sp.pep.is_file():
            raise SystemExit(f"ERROR: {sp.short}: protein FASTA {sp.pep} not found "
                             "(run bin/build_study_config.py first)")
        sp.n_proteins = count_fasta(sp.pep)
        sp.chunks = plan_chunks(sp.n_proteins, args.max_proteins_per_job)
        tasks += [(sp.short, i, e - s) for i, (s, e) in enumerate(sp.chunks)]
        print(f"{sp.short}: {sp.n_proteins} proteins -> {len(sp.chunks)} chunk(s) ({sp.mode}); "
              f"output {sp.annotation}", file=sys.stderr)
    cap = batch_capacity(args.batch_minutes, args.sec_per_protein)
    batches = pack_batches(tasks, cap)
    total = sum(t[2] for t in tasks)
    for b, items in enumerate(batches):
        n = sum(t[2] for t in items)
        print(f"batch {b:03d}: {len(items)} run(s), {n} proteins, about {n * args.sec_per_protein / 60:.0f} min: "
              + ", ".join(f"{s}#{c}" for s, c, _ in items), file=sys.stderr)
    print(f"total: {len(selected)} species, {total} proteins, {len(tasks)} hmmsearch run(s), "
          f"{len(batches)} batch job(s) (cap {cap} proteins = {args.batch_minutes:g} min), "
          f"about {total * args.sec_per_protein / 3600:.1f} h hmmsearch + 1 finalize job", file=sys.stderr)
    if args.dry_run:
        return 0

    tag = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    batch_dir = work / "batches" / tag
    guard(batch_dir)
    (batch_dir / "logs").mkdir(parents=True, exist_ok=True)
    by_short = {sp.short: sp for sp in selected}
    chunk_paths: dict[tuple[str, int], Path] = {}
    for sp in selected:
        d = guard(work / sp.short)
        (d / "logs").mkdir(parents=True, exist_ok=True)
        for i, (s, e) in enumerate(sp.chunks):
            p = d / f"chunk_{i:03d}.sbatch"
            write_text_guarded(guard, p, chunk_script(sp, i, s, e, work, args.pfam_hmm, args.cpus,
                                                      args.partition, args.constraint, args.sec_per_protein))
            chunk_paths[(sp.short, i)] = p
    batch_paths = []
    batch_of: dict[tuple[str, int], int] = {}
    for b, items in enumerate(batches):
        p = batch_dir / f"batch_{b:03d}.sbatch"
        write_text_guarded(guard, p, batch_script(
            b, [(s, c, n, chunk_paths[(s, c)]) for s, c, n in items], batch_dir, args.cpus,
            args.partition, args.constraint, args.sec_per_protein))
        batch_paths.append(p)
        for s, c, _ in items:
            batch_of[(s, c)] = b
    shorts = [sp.short for sp in selected]
    fin = batch_dir / "finalize.sbatch"
    write_text_guarded(guard, fin, finalize_script(shorts, args.study, work, ann_dir, args.pfam_hmm,
                                                   args.force, batch_dir, partition=args.partition,
                                                   constraint=args.constraint,
                                                   fill_empty_uniprot=args.fill_empty_uniprot))
    job_ids: list[str] = []
    fin_id = None
    if args.submit:
        for p in batch_paths:
            res = subprocess.run(["sbatch", "--parsable", str(p)], capture_output=True, text=True, check=True)
            job_ids.append(res.stdout.strip().split(";")[0])
        res = subprocess.run(["sbatch", "--parsable", f"--dependency=afterany:{':'.join(job_ids)}", str(fin)],
                             capture_output=True, text=True, check=True)
        fin_id = res.stdout.strip().split(";")[0]
        print(f"submitted batch jobs {job_ids}, finalize {fin_id}", file=sys.stderr)
    created = datetime.datetime.now().isoformat(timespec="seconds")
    write_text_guarded(guard, batch_dir / "batches.json", json.dumps({
        "study": args.study, "created": created, "species": shorts, "batch_capacity_proteins": cap,
        "batch_minutes": args.batch_minutes, "sec_per_protein": args.sec_per_protein,
        "batches": [{"script": str(p), "job_id": job_ids[b] if job_ids else None,
                     "chunks": [[s, c, n] for s, c, n in batches[b]]} for b, p in enumerate(batch_paths)],
        "finalize_script": str(fin), "finalize_job_id": fin_id}, indent=1) + "\n")
    for sp in selected:
        plan = {"study": args.study, "short": sp.short, "protein_accession": sp.protein_accession,
                "protein_source": sp.protein_source, "mode": sp.mode, "pep": str(sp.pep),
                "n_proteins": sp.n_proteins, "chunks": sp.chunks, "pfam_hmm": args.pfam_hmm,
                "max_proteins_per_job": args.max_proteins_per_job, "sec_per_protein": args.sec_per_protein,
                "cpus": args.cpus, "annotation": str(sp.annotation), "created": created,
                "chunk_scripts": [str(chunk_paths[(sp.short, i)]) for i in range(len(sp.chunks))],
                "batch_dir": str(batch_dir),
                "chunk_batch": [batch_of[(sp.short, i)] for i in range(len(sp.chunks))],
                "job_ids": sorted({job_ids[batch_of[(sp.short, i)]] for i in range(len(sp.chunks))})
                if job_ids else [],
                "finalize_script": str(fin), "finalize_job_id": fin_id}
        write_text_guarded(guard, work / sp.short / "plan.json", json.dumps(plan, indent=1) + "\n")
    if args.local:
        if not os.environ.get("SCRATCH"):
            raise SystemExit("ERROR: --local needs $SCRATCH (run inside a SLURM job)")
        for sp in selected:
            for i in range(len(sp.chunks)):
                subprocess.run(["bash", str(chunk_paths[(sp.short, i)])], check=True)
            finalize(by_short[sp.short], work, args.pfam_hmm, args.force, guard, orig_dir)
    elif not args.submit:
        for p in batch_paths:
            print(f"  sbatch {p}", file=sys.stderr)
        print(f"  then: sbatch --dependency=afterany:<batch job ids> {fin}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
