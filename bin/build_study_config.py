#!/usr/bin/env python3
"""Drive a study's species.csv through per-species source dispatch and emit a
nf_NovInvenio-ready run config.

Input: <study-dir>/species.csv, columns:
    Short,Species,Strain,Group,TaxonGroup,
    Protein_Source,Protein_Accession,Taxon_ID,
    Genome_Source,Genome_Accession,
    GFF3_Source,GFF3_Accession

Protein_Source/Genome_Source/GFF3_Source are independent per row -- a species'
protein and genome can come from unrelated places (e.g. an already-local protein
FASTA paired with a genome that still needs an NCBI fetch). Recognized values:

  Protein_Source:
    uniprot    -- Protein_Accession = UniProt proteome ID (Taxon_ID also required).
                  Fetched via fetch_uniprot_proteome.py; also drives GO/Pfam/
                  InterPro/gene-name annotation extraction (extract_dat_annotations.py)
                  -- the only Protein_Source that does.
    ncbi       -- Protein_Accession = GCA/GCF assembly accession. Protein FASTA
                  included in the same NCBI Datasets genome package (fetched via
                  fetch_genome_assembly.py --include-protein), with provenance record.
    local_faa  -- Protein_Accession = path to an existing protein FASTA. Copied in
                  directly, with a provenance record (no fetch).

  Genome_Source:
    ncbi          -- Genome_Accession = GCA/GCF assembly accession. Fetched via
                     fetch_genome_assembly.py (genome + GFF3 together).
    local_genome  -- Genome_Accession = path to an existing genome FASTA. Copied in
                     directly, with a provenance record (no fetch).
    (blank)       -- no genome for this row; DNA config.csv cell left empty.

  GFF3_Source (optional; blank means "auto" -- see below):
    local_gff3 -- GFF3_Accession = path to an existing GFF3. Always used, regardless
                  of Genome_Source.
    none       -- force the GFF3 config.csv cell empty even if Genome_Source=ncbi
                  would otherwise have supplied one.
    (blank)    -- auto: if Genome_Source=ncbi, use whatever GFF3 that NCBI package
                  included (today's behavior); otherwise empty.

Stem naming: a uniprot-sourced protein keeps the existing
"{Proteome_ID}_{Taxon_ID}" stem (self-documenting, collision-proof across
strains sharing a Short by mistake). Every other row uses Short as its stem --
Short is already required to be unique per study.

For each species this:
  1. Resolves Protein (fetch or local copy).
  2. Resolves Genome+GFF3 (fetch or local copy/copies).
  3. If Protein_Source == "uniprot", runs bin/extract_dat_annotations.py against
     the cached .dat.gz, writing <study-dir>/annotations/<Protein_Accession>.tsv
     (gitignored, regenerable, produced unconditionally -- even with --skip-fetch
     -- so bin/run_study.sh's post-run report sync always has current annotation
     to merge). Rows sourced any other way get no annotation file (unchanged from
     today's local-file studies).
  4. Writes <study-dir>/config.csv (GROUP,Species,Strain,Protein,DNA,GFF3,Short,
     TaxonGroup) with basenames only, and <study-dir>/DATA_MANIFEST.yaml
     aggregating every species' provenance record plus a record for this
     generation step itself.

Re-running with --skip-fetch rebuilds config.csv/data_dir from whatever's already
cached (uniprot/ncbi rows) without re-downloading; local_* rows are always
re-copied (there's nothing to "skip" -- they were never fetched).
"""
import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from provenance import append_manifest, build_record, sha256_of  # noqa: E402

BIN = Path(__file__).parent
LOCAL_LICENSE_DEFAULT = "Internal / unpublished (not yet released outside this project)"


def run(cmd: list[str]) -> None:
    print(f"$ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def gunzip_to(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def load_provenance(sidecar: Path) -> dict:
    with open(sidecar) as fh:
        return yaml.safe_load(fh)


def _local_copy_record(src: Path, dest: Path, license_: str) -> dict:
    return build_record(
        source_url=f"(internal -- {src})",
        source_release="local file, not independently versioned",
        license=license_,
        local_path=dest,
        checksum=sha256_of(dest),
        derived_by=f"copied from {src} via bin/build_study_config.py",
    )


def resolve_protein(row, args, pep_dir, manifest_records):
    """Returns (stem, protein_basename, dat_gz_path_or_None)."""
    short = row["Short"]
    source = row["Protein_Source"]
    accession = row["Protein_Accession"]

    if source == "uniprot":
        taxid = row["Taxon_ID"]
        if not args.skip_fetch:
            run([
                sys.executable, str(BIN / "fetch_uniprot_proteome.py"),
                "--proteome-id", accession, "--taxid", taxid,
                "--outdir", args.uniprot_cache, "--short", short,
            ])
        stem = f"{accession}_{taxid}"
        fasta_gz = Path(args.uniprot_cache) / accession / f"{stem}.fasta.gz"
        dat_gz = Path(args.uniprot_cache) / accession / f"{stem}.dat.gz"
        if not fasta_gz.exists():
            sys.exit(f"[{short}] ERROR: expected {fasta_gz} not found -- run without --skip-fetch first")
        pep_dir.mkdir(parents=True, exist_ok=True)
        pep_out = pep_dir / f"{stem}.pep.fa"
        gunzip_to(fasta_gz, pep_out)
        for sidecar in (fasta_gz, dat_gz):
            sc = sidecar.with_suffix(sidecar.suffix + ".provenance.yaml")
            if sc.exists():
                manifest_records.append(load_provenance(sc))
        return stem, pep_out.name, (dat_gz if dat_gz.exists() else None)

    if source == "ncbi":
        if not args.skip_fetch:
            run([
                sys.executable, str(BIN / "fetch_genome_assembly.py"),
                "--accession", accession, "--outdir", args.ncbi_cache, "--short", short,
                "--include-protein",
            ])
        genome_dir = Path(args.ncbi_cache) / accession / "extracted" / "ncbi_dataset" / "data" / accession
        faa_candidates = sorted(genome_dir.glob("*.faa")) if genome_dir.exists() else []
        if not faa_candidates:
            sys.exit(f"[{short}] ERROR: expected protein FASTA under {genome_dir} not found -- run without --skip-fetch first")
        stem = short
        pep_dir.mkdir(parents=True, exist_ok=True)
        pep_out = pep_dir / f"{stem}.pep.fa"
        shutil.copyfile(faa_candidates[0], pep_out)
        # Load provenance for the protein file
        sc = faa_candidates[0].with_suffix(faa_candidates[0].suffix + ".provenance.yaml")
        if sc.exists():
            manifest_records.append(load_provenance(sc))
        return stem, pep_out.name, None

    if source == "local_faa":
        src = Path(accession)
        if not src.exists():
            sys.exit(f"[{short}] ERROR: expected local protein FASTA {src} not found")
        stem = short
        pep_dir.mkdir(parents=True, exist_ok=True)
        pep_out = pep_dir / f"{stem}.pep.fa"
        shutil.copyfile(src, pep_out)
        manifest_records.append(_local_copy_record(src, pep_out, args.local_license))
        return stem, pep_out.name, None

    sys.exit(f"[{short}] ERROR: unknown Protein_Source {source!r}")


def resolve_genome_and_gff3(row, args, stem, dna_dir, gff3_dir, manifest_records):
    """Returns (dna_basename, gff3_basename) -- either may be ''."""
    short = row["Short"]
    gsource = row["Genome_Source"]
    gaccession = row["Genome_Accession"]
    gffsource = row["GFF3_Source"]
    gffaccession = row["GFF3_Accession"]

    dna_out_name = ""
    gff3_out_name = ""

    if gsource == "ncbi":
        if not args.skip_fetch:
            run([
                sys.executable, str(BIN / "fetch_genome_assembly.py"),
                "--accession", gaccession, "--outdir", args.ncbi_cache, "--short", short,
            ])
        genome_dir = Path(args.ncbi_cache) / gaccession / "extracted" / "ncbi_dataset" / "data" / gaccession
        fna_candidates = sorted(genome_dir.glob("*.fna")) if genome_dir.exists() else []
        gff_candidates = sorted(genome_dir.glob("*.gff")) if genome_dir.exists() else []
        if not fna_candidates:
            sys.exit(f"[{short}] ERROR: expected genome under {genome_dir} not found -- run without --skip-fetch first")
        dna_dir.mkdir(parents=True, exist_ok=True)
        dna_out = dna_dir / f"{stem}.dna.fa"
        shutil.copyfile(fna_candidates[0], dna_out)
        dna_out_name = dna_out.name
        # Load provenance for the genome file (always used)
        sc = fna_candidates[0].with_suffix(fna_candidates[0].suffix + ".provenance.yaml")
        if sc.exists():
            manifest_records.append(load_provenance(sc))
        # Determine which GFF3 to use, if any: gffsource=="none" suppresses NCBI GFF3
        if gffsource == "none":
            gff_candidates = []
        if gffsource == "" and gff_candidates:
            gff3_dir.mkdir(parents=True, exist_ok=True)
            gff3_out_name = f"{stem}.gff3"
            shutil.copyfile(gff_candidates[0], gff3_dir / gff3_out_name)
            # Load provenance for the GFF3 file (only if we actually use it)
            sc = gff_candidates[0].with_suffix(gff_candidates[0].suffix + ".provenance.yaml")
            if sc.exists():
                manifest_records.append(load_provenance(sc))
    elif gsource == "local_genome":
        src = Path(gaccession)
        if not src.exists():
            sys.exit(f"[{short}] ERROR: expected local genome FASTA {src} not found")
        dna_dir.mkdir(parents=True, exist_ok=True)
        dna_out = dna_dir / f"{stem}.dna.fa"
        shutil.copyfile(src, dna_out)
        dna_out_name = dna_out.name
        manifest_records.append(_local_copy_record(src, dna_out, args.local_license))
    elif gsource == "":
        pass
    else:
        sys.exit(f"[{short}] ERROR: unknown Genome_Source {gsource!r}")

    if gffsource == "local_gff3":
        src = Path(gffaccession)
        if not src.exists():
            sys.exit(f"[{short}] ERROR: expected local GFF3 {src} not found")
        gff3_dir.mkdir(parents=True, exist_ok=True)
        gff3_out = gff3_dir / f"{stem}.gff3"
        shutil.copyfile(src, gff3_out)
        gff3_out_name = gff3_out.name
        manifest_records.append(_local_copy_record(src, gff3_out, args.local_license))
    elif gffsource not in ("", "none"):
        sys.exit(f"[{short}] ERROR: unknown GFF3_Source {gffsource!r}")

    return dna_out_name, gff3_out_name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study-dir", required=True, help="e.g. studies/fungi/pezizo_set1")
    ap.add_argument("--uniprot-cache", default="data/uniprot")
    ap.add_argument("--ncbi-cache", default="data/ncbi")
    ap.add_argument("--skip-fetch", action="store_true", help="Rebuild config/data_dir from an already-populated cache, no downloads")
    ap.add_argument("--local-license", default=LOCAL_LICENSE_DEFAULT)
    args = ap.parse_args()

    study_dir = Path(args.study_dir)
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        sys.exit(f"ERROR: {species_csv} not found")

    link_dir = study_dir / "data_dir"
    pep_dir, dna_dir, gff3_dir = (link_dir / d for d in ("pep", "dna", "gff3"))

    config_rows = []
    manifest_records = []
    seen_stems: dict[str, str] = {}

    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            short = row["Short"]
            stem, protein_name, dat_gz = resolve_protein(row, args, pep_dir, manifest_records)

            if stem in seen_stems:
                sys.exit(
                    f"ERROR: {species_csv} resolves two rows to the same stem {stem!r} "
                    f"(Short={seen_stems[stem]!r} and Short={short!r})"
                )
            seen_stems[stem] = short

            dna_name, gff3_name = resolve_genome_and_gff3(row, args, stem, dna_dir, gff3_dir, manifest_records)

            if row["Protein_Source"] == "uniprot" and dat_gz is not None:
                annotations_dir = study_dir / "annotations"
                annotations_dir.mkdir(parents=True, exist_ok=True)
                run([
                    sys.executable, str(BIN / "extract_dat_annotations.py"),
                    "--dat-gz", str(dat_gz),
                    "--output", str(annotations_dir / f"{row['Protein_Accession']}.tsv"),
                ])

            config_rows.append({
                "GROUP": row["Group"],
                "Species": row["Species"],
                "Strain": row["Strain"],
                "Protein": protein_name,
                "DNA": dna_name,
                "GFF3": gff3_name,
                "Short": short,
                "TaxonGroup": row["TaxonGroup"],
            })

    config_csv = study_dir / "config.csv"
    with open(config_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["GROUP", "Species", "Strain", "Protein", "DNA", "GFF3", "Short", "TaxonGroup"])
        w.writeheader()
        w.writerows(config_rows)

    manifest_records.append(build_record(
        source_url="(derived, no external source)",
        source_release="n/a",
        license="n/a",
        local_path=config_csv,
        derived_by=f"bin/build_study_config.py --study-dir {study_dir}",
    ))
    append_manifest(manifest_records, study_dir / "DATA_MANIFEST.yaml")

    print(f"\nWrote {config_csv} ({len(config_rows)} species)", file=sys.stderr)
    print(f"Wrote {link_dir} (pep/dna/gff3 -- pass as --data_dir to nf_NovInvenio)", file=sys.stderr)
    print(f"Wrote {study_dir / 'DATA_MANIFEST.yaml'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
