#!/usr/bin/env python3
"""Pull many genome assemblies (DNA + GFF3[+protein]) in one `datasets download`
call per chunk, instead of bin/fetch_genome_assembly.py's one-accession-per-call.

Why this exists: a study with many NCBI-sourced rows (e.g. a pangenome study
enumerated by `bin/ni discover`) drove bin/build_study_config.py to shell out to
`datasets download genome accession <single accession>` once per row -- 293
separate NCBI Datasets API round-trips for a 293-strain study. The `datasets` CLI
accepts multiple accessions in one invocation and returns one zip with a
per-accession `ncbi_dataset/data/<accession>/` subfolder each, identical in shape
to a single-accession package -- so batching changes only how many network round
trips happen, not the on-disk layout bin/build_study_config.py's glob-based
lookup depends on.

Chunked at --batch-size (default 75): NCBI does not document a hard limit on
accessions-per-call, but an unbounded single call risks a very large zip and a
correspondingly long single point of failure (one bad accession or one timeout
loses the whole batch). 75 keeps each call's zip in a modest size range for
typical fungal-genome studies while still cutting round-trips by ~75x for a
293-genome study (4 calls instead of 293).

This is a separate script from bin/fetch_genome_assembly.py, not an extension of
it: fetch_genome_assembly.py's single-accession contract (one accession in, one
--short label for its log line, one call site in
bin/build_study_config.py's per-row path for local_genome/local_faa rows that
never call it) stays untouched for every existing caller. This script is only
ever invoked once per build_study_config.py run, for the whole set of
Genome_Source=ncbi / Protein_Source=ncbi accessions at once.

Requires the `datasets` CLI (bioconda: ncbi-datasets-cli) on PATH, same as
fetch_genome_assembly.py.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from provenance import build_record, sha256_of, write_record  # noqa: E402

LICENSE = "Public Domain (NCBI, https://www.ncbi.nlm.nih.gov/home/about/policies/)"


def check_datasets_on_path() -> str:
    exe = shutil.which("datasets")
    if not exe:
        sys.exit(
            "ERROR: `datasets` (ncbi-datasets-cli) not found on PATH.\n"
            "Install it (bioconda: ncbi-datasets-cli) -- see NII's pixi.toml."
        )
    return exe


def datasets_version(exe: str) -> str:
    out = subprocess.run([exe, "--version"], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def chunked(items: list[str], size: int) -> list[list[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def fetch_batch(
    accessions: list[str],
    *,
    outdir: Path,
    include_protein: bool,
    batch_size: int,
) -> None:
    """Download every accession in `accessions` into outdir/<accession>/extracted/...,
    the same per-accession layout bin/fetch_genome_assembly.py produces -- callers
    (bin/build_study_config.py) glob `outdir/<accession>/extracted/ncbi_dataset/data/<accession>/*.fna`
    identically regardless of which script fetched it.
    """
    exe = check_datasets_on_path()
    include_list = "genome,gff3,protein" if include_protein else "genome,gff3"
    version = datasets_version(exe)

    for batch_num, batch in enumerate(chunked(accessions, batch_size), start=1):
        print(
            f"[batch {batch_num}] datasets download genome accession "
            f"{' '.join(batch)} --include {include_list} ({len(batch)} accessions)",
            file=sys.stderr,
        )
        batch_zip = outdir / f"_batch_{batch_num:03d}.zip"
        batch_zip.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                exe, "download", "genome", "accession", *batch,
                "--include", include_list,
                "--filename", str(batch_zip),
            ],
            check=True,
        )

        extract_dir = outdir / f"_batch_{batch_num:03d}_extracted"
        with zipfile.ZipFile(batch_zip) as zf:
            zf.extractall(extract_dir)
        batch_zip.unlink()

        for accession in batch:
            src_data_dir = extract_dir / "ncbi_dataset" / "data" / accession
            if not src_data_dir.exists():
                print(
                    f"[{accession}] WARNING: not found in batch {batch_num}'s package "
                    "-- accession may be invalid, suppressed, or superseded",
                    file=sys.stderr,
                )
                continue

            # Re-shape into the same per-accession layout
            # fetch_genome_assembly.py produces, so build_study_config.py's glob
            # (outdir/<accession>/extracted/ncbi_dataset/data/<accession>/*.fna)
            # finds it identically regardless of which script fetched it.
            dest_data_dir = outdir / accession / "extracted" / "ncbi_dataset" / "data" / accession
            dest_data_dir.parent.mkdir(parents=True, exist_ok=True)
            if dest_data_dir.exists():
                shutil.rmtree(dest_data_dir)
            shutil.move(str(src_data_dir), str(dest_data_dir))

            fna_candidates = sorted(dest_data_dir.glob("*.fna"))
            gff_candidates = sorted(dest_data_dir.glob("*.gff"))
            faa_candidates = sorted(dest_data_dir.glob("*.faa")) if include_protein else []
            if not fna_candidates:
                sys.exit(f"[{accession}] ERROR: no .fna found under {dest_data_dir} -- datasets output layout may have changed.")
            if include_protein and not faa_candidates:
                print(
                    f"[{accession}] WARNING: --include-protein requested but no protein "
                    "in this assembly's Datasets package",
                    file=sys.stderr,
                )

            for f, kind in (
                (fna_candidates[0], "genome"),
                (gff_candidates[0] if gff_candidates else None, "gff3"),
                (faa_candidates[0] if faa_candidates else None, "protein"),
            ):
                if f is None:
                    continue
                checksum = sha256_of(f)
                record = build_record(
                    source_url=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/",
                    source_release=f"NCBI Datasets CLI {version}, accession {accession}",
                    license=LICENSE,
                    local_path=f,
                    checksum=checksum,
                    extra={"kind": kind, "accession": accession},
                )
                write_record(record, f.with_suffix(f.suffix + ".provenance.yaml"))
                print(f"[{accession}] OK  {kind}: {f.name}  sha256={checksum[:12]}...", file=sys.stderr)

        shutil.rmtree(extract_dir)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", action="append", required=True, dest="accessions",
                     help="GCA/GCF accession; repeat for each genome (e.g. --accession GCA_1 --accession GCA_2)")
    ap.add_argument("--outdir", default="data/ncbi", help="Root output directory (default: data/ncbi)")
    ap.add_argument("--include-protein", action="store_true", help="Also download each assembly's protein FASTA, if available")
    ap.add_argument("--batch-size", type=int, default=75, help="Accessions per `datasets download` call (default: 75)")
    args = ap.parse_args()

    if args.batch_size < 1:
        sys.exit("ERROR: --batch-size must be >= 1")

    # De-dup while preserving order -- a study can list the same accession for
    # both a Protein_Source=ncbi and Genome_Source=ncbi row (the common case),
    # and re-fetching it twice in the same run would be pure waste.
    seen = set()
    accessions = [a for a in args.accessions if not (a in seen or seen.add(a))]

    fetch_batch(
        accessions,
        outdir=Path(args.outdir),
        include_protein=args.include_protein,
        batch_size=args.batch_size,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
