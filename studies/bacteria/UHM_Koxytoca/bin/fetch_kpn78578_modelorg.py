#!/usr/bin/env python3
"""One-off pull for the Klebsiella pneumoniae MGH 78578 model organism, used by
studies/bacteria/UHM_Koxytoca/modelorgs.yaml (Kpn78578 entry).

Fetches both sources agreed on for this model organism (see DESIGN discussion):
  1. NCBI RefSeq assembly GCF_000016305.1 (genome + GFF3 + protein.faa, all 6
     replicons: chromosome NC_009648.1 + 5 plasmids incl. NC_009653.1) via
     `datasets download genome ... --include genome,gff3,protein` -- the plain
     bin/fetch_genome_assembly.py recipe only pulls genome+gff3 (no protein), so
     this is its own small pull rather than an addition to that shared script.
  2. UniProt reference proteome UP000000265 (taxid 272620) -- confirmed via
     https://rest.uniprot.org/proteomes/search?query=taxonomy_id:573 to be the
     exact MGH 78578 strain proteome (not just the species-level K. pneumoniae
     taxon 573), fetched via the existing bin/fetch_uniprot_proteome.py.

Both land in the ephemeral, gitignored data/ caches with provenance sidecars (never
committed raw -- DESIGN.md Sec 4). NCBI's GCF_000016305.1 pull (genome+GFF3+protein,
all 6 replicons) is kept as reference data only -- it is NOT what modelorgs.yaml's
diamond_fasta lookup uses, because its protein_fasta and gene_names_csv would then
live in two different ID spaces (NCBI WP_ accessions vs. UniProt accessions from
bin/extract_dat_annotations.py's own .dat-derived TSV). Instead this script builds
the diamond_fasta lookup entirely from the UniProt UP000000265 pull, single ID space
throughout:
  - config_support/modelorgs/Kpn78578_protein.faa -- the UniProt fasta with headers
    rewritten to ">ACCESSION gene=ACCESSION" (UniProt's own header has no "gene="
    key=value field lib/model_organisms.py's diamond_fasta id_transform expects --
    this makes the lookup a same-ID round-trip instead of a no-op failure)
  - config_support/modelorgs/Kpn78578_gene_names_UniProt.tsv (UniProt DR-derived
    gene_name/description/GO/Pfam, via bin/extract_dat_annotations.py, keyed by the
    same UniProt accession -- modelorgs.yaml's gene_names_csv)

Requires the `datasets` CLI (module load ncbi_datasets, or pixi's ncbi-datasets-cli).

This is a study-specific script (lives under this study's own bin/, not NII's shared
bin/ -- see CLAUDE.md's "Where new code goes"), so it hardcodes NII_ROOT below rather
than deriving paths from __file__/cwd; NII_BIN is where the shared, reused-as-is
fetch_uniprot_proteome.py / extract_dat_annotations.py scripts still live.
"""
import gzip
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
NII_BIN = NII_ROOT / "bin"
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import build_record, sha256_of, write_record  # noqa: E402

NCBI_ACCESSION = "GCF_000016305.1"
UNIPROT_PROTEOME_ID = "UP000000265"
UNIPROT_TAXID = "272620"
NCBI_LICENSE = "Public Domain (NCBI, https://www.ncbi.nlm.nih.gov/home/about/policies/)"


def check_datasets_on_path() -> str:
    exe = shutil.which("datasets")
    if not exe:
        sys.exit(
            "ERROR: `datasets` (ncbi-datasets-cli) not found on PATH.\n"
            "Run: source /etc/profile.d/modules.sh && module load ncbi_datasets/18.30.1"
        )
    return exe


def datasets_version(exe: str) -> str:
    out = subprocess.run([exe, "--version"], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def fetch_ncbi(outdir: Path) -> Path:
    exe = check_datasets_on_path()
    acc_dir = outdir / NCBI_ACCESSION
    acc_dir.mkdir(parents=True, exist_ok=True)
    zip_path = acc_dir / f"{NCBI_ACCESSION}.zip"

    print(f"[Kpn78578] datasets download genome accession {NCBI_ACCESSION} --include genome,gff3,protein", file=sys.stderr)
    subprocess.run(
        [exe, "download", "genome", "accession", NCBI_ACCESSION,
         "--include", "genome,gff3,protein", "--filename", str(zip_path)],
        check=True,
    )

    extract_dir = acc_dir / "extracted"
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(extract_dir)

    data_dir = extract_dir / "ncbi_dataset" / "data" / NCBI_ACCESSION
    protein_faa = data_dir / "protein.faa"
    if not protein_faa.exists():
        sys.exit(f"ERROR: no protein.faa found under {data_dir}")

    for f, kind in (
        (next(data_dir.glob("*.fna")), "genome"),
        (data_dir / "genomic.gff", "gff3"),
        (protein_faa, "protein"),
    ):
        checksum = sha256_of(f)
        record = build_record(
            source_url=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{NCBI_ACCESSION}/",
            source_release=f"NCBI Datasets CLI {datasets_version(exe)}, accession {NCBI_ACCESSION}",
            license=NCBI_LICENSE,
            local_path=f,
            checksum=checksum,
            extra={"kind": kind, "accession": NCBI_ACCESSION},
        )
        write_record(record, f.with_suffix(f.suffix + ".provenance.yaml"))
        print(f"[Kpn78578] OK  {kind}: {f.name}  sha256={checksum[:12]}...", file=sys.stderr)

    return protein_faa


def fetch_uniprot(outdir: Path) -> Path:
    subprocess.run(
        [sys.executable, str(NII_BIN / "fetch_uniprot_proteome.py"),
         "--proteome-id", UNIPROT_PROTEOME_ID, "--taxid", UNIPROT_TAXID,
         "--outdir", str(outdir), "--short", "Kpn78578"],
        check=True,
    )
    stem = f"{UNIPROT_PROTEOME_ID}_{UNIPROT_TAXID}"
    dat_gz = outdir / UNIPROT_PROTEOME_ID / f"{stem}.dat.gz"
    if not dat_gz.exists():
        sys.exit(f"ERROR: expected {dat_gz} not found after UniProt fetch")
    return dat_gz


UNIPROT_HEADER_RE = re.compile(r"^>(?:sp|tr)\|([A-Za-z0-9]+)\|\S+")


def write_diamond_ready_fasta(uniprot_fasta_gz: Path, out_path: Path) -> None:
    """Rewrite UniProt's '>tr|ACCESSION|ENTRY_NAME description...' headers to
    '>ACCESSION gene=ACCESSION' -- see module docstring for why."""
    n = 0
    with gzip.open(uniprot_fasta_gz, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                m = UNIPROT_HEADER_RE.match(line)
                if not m:
                    sys.exit(f"ERROR: unrecognized UniProt FASTA header: {line!r}")
                acc = m.group(1)
                fout.write(f">{acc} gene={acc}\n")
                n += 1
            else:
                fout.write(line)
    print(f"[Kpn78578] wrote {out_path} ({n} proteins, headers rewritten for diamond_fasta)", file=sys.stderr)


def main() -> int:
    ncbi_cache = NII_ROOT / "data" / "ncbi"
    uniprot_cache = NII_ROOT / "data" / "uniprot"
    tracked_dir = NII_ROOT / "config_support" / "modelorgs"
    tracked_dir.mkdir(parents=True, exist_ok=True)

    fetch_ncbi(ncbi_cache)  # reference data only, see module docstring
    dat_gz = fetch_uniprot(uniprot_cache)
    uniprot_fasta_gz = dat_gz.with_name(dat_gz.name.replace(".dat.gz", ".fasta.gz"))

    tracked_protein = tracked_dir / "Kpn78578_protein.faa"
    write_diamond_ready_fasta(uniprot_fasta_gz, tracked_protein)

    gene_names_tsv = tracked_dir / "Kpn78578_gene_names_UniProt.tsv"
    subprocess.run(
        [sys.executable, str(NII_BIN / "extract_dat_annotations.py"),
         "--dat-gz", str(dat_gz), "--output", str(gene_names_tsv)],
        check=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
