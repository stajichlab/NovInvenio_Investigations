#!/usr/bin/env python3
"""Pull one genome assembly (DNA + GFF3) via NCBI Datasets, keyed by GCA/GCF accession.

Added per the independent (Fable-model) design review: UniProt's per-proteome
.fasta.gz/.dat.gz pull has no DNA, but TBLASTN validation (and GFF3 chrom/start
columns) need the genome assembly. Using the *same* GCA accession UniProt's own
proteome metadata reports (genomeAssembly.assemblyId) guarantees the genome and the
protein annotation are version-matched -- a mismatch here would manufacture false
"genome hit, no protein call" TBLASTN signals, which is exactly the failure mode
this script exists to avoid. See DESIGN.md Sec 5.

Ephemeral/recipe-driven pull (DESIGN.md Sec 4): downloads land under --outdir
(default data/ncbi/, gitignored) and are never archived.

Requires the `datasets` CLI (bioconda: ncbi-datasets-cli) on PATH.
"""
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accession", required=True, help="GenBank/RefSeq assembly accession, e.g. GCA_000182925.2")
    ap.add_argument("--outdir", default="data/ncbi", help="Root output directory (default: data/ncbi)")
    ap.add_argument("--short", help="Short species code, used only for a human-readable log line")
    ap.add_argument("--include-protein", action="store_true", help="Also download the protein FASTA from this assembly (if available)")
    args = ap.parse_args()

    exe = check_datasets_on_path()
    label = args.short or args.accession
    outdir = Path(args.outdir) / args.accession
    outdir.mkdir(parents=True, exist_ok=True)
    zip_path = outdir / f"{args.accession}.zip"

    include_list = "genome,gff3"
    if args.include_protein:
        include_list = "genome,gff3,protein"
    print(f"[{label}] datasets download genome accession {args.accession} --include {include_list}", file=sys.stderr)
    subprocess.run(
        [
            exe, "download", "genome", "accession", args.accession,
            "--include", include_list,
            "--filename", str(zip_path),
        ],
        check=True,
    )

    extract_dir = outdir / "extracted"
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(extract_dir)

    data_dir = extract_dir / "ncbi_dataset" / "data" / args.accession
    fna_candidates = sorted(data_dir.glob("*.fna"))
    gff_candidates = sorted(data_dir.glob("*.gff"))
    faa_candidates = sorted(data_dir.glob("*.faa")) if args.include_protein else []
    if not fna_candidates:
        sys.exit(f"[{label}] ERROR: no .fna found under {data_dir} -- datasets output layout may have changed.")
    if args.include_protein and not faa_candidates:
        sys.exit(f"[{label}] ERROR: --include-protein was requested but no .faa found under {data_dir} -- datasets output layout may have changed.")

    fna = fna_candidates[0]
    gff = gff_candidates[0] if gff_candidates else None
    faa = faa_candidates[0] if faa_candidates else None

    records = []
    for f, kind in ((fna, "genome"), (gff, "gff3"), (faa, "protein")):
        if f is None:
            if kind == "protein" and args.include_protein:
                print(f"[{label}] WARNING: --include-protein requested but no protein in this assembly's Datasets package", file=sys.stderr)
            elif kind == "gff3":
                print(f"[{label}] WARNING: no GFF3 in this assembly's Datasets package", file=sys.stderr)
            continue
        checksum = sha256_of(f)
        record = build_record(
            source_url=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{args.accession}/",
            source_release=f"NCBI Datasets CLI {datasets_version(exe)}, accession {args.accession}",
            license=LICENSE,
            local_path=f,
            checksum=checksum,
            extra={"kind": kind, "accession": args.accession},
        )
        write_record(record, f.with_suffix(f.suffix + ".provenance.yaml"))
        records.append(record)
        print(f"[{label}] OK  {kind}: {f.name}  sha256={checksum[:12]}...", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
