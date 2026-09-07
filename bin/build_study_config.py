#!/usr/bin/env python3
"""Drive a study's species.csv through both pull scripts and emit a nf_NovInvenio-ready
run config.

Input: <study-dir>/species.csv, columns:
    Short,Species,Strain,Group,TaxonGroup,UniProt_Proteome_ID,Taxon_ID,GCA_Accession

For each species this:
  1. Runs bin/fetch_uniprot_proteome.py (protein FASTA + taxonomy/GO/Pfam/InterPro .dat)
  2. Runs bin/fetch_genome_assembly.py (DNA + GFF3, keyed off the same UniProt-reported
     GCA accession, so genome and protein annotation are version-matched -- see
     DESIGN.md Sec 5)
  3. Materializes plain (decompressed) files into <study-dir>/data_dir/{pep,dna,gff3}/,
     named by the UniProt proteome stem ("{Proteome_ID}_{Taxon_ID}", e.g.
     "UP001658139_2528406") rather than Short -- self-documenting back to the exact
     UniProt record regardless of what Short a study happens to assign, and
     inherently unique even across multiple strains of the same species (each strain
     is its own UniProt proteome, hence its own stem) without relying on the study
     author never reusing a Short. config.csv's Protein/DNA/GFF3 columns carry
     whatever this actually materializes to -- nf_NovInvenio's resolve_fa() only
     needs those basenames to exist under --data_dir, not to equal Short.
  4. Runs bin/extract_dat_annotations.py against each species' cached .dat.gz, writing
     <study-dir>/annotations/<UniProt_Proteome_ID>.tsv (gene_name/description/GO/Pfam/
     InterPro per accession) -- gitignored, regenerable from the cache in ~2s/species,
     but produced here unconditionally (even with --skip-fetch) so bin/run_study.sh's
     post-run report sync (bin/sync_reports.sh) always has current annotation to merge.
  5. Writes <study-dir>/config.csv (GROUP,Species,Strain,Protein,DNA,GFF3,Short,
     TaxonGroup) with basenames only, and <study-dir>/DATA_MANIFEST.yaml aggregating
     every species' provenance record plus a record for this generation step itself.

Re-running with --skip-fetch rebuilds config.csv/data_dir from whatever's already in
the ephemeral data/ cache, without re-downloading (useful after the first pull, or
when only the species list/grouping changed, not the source data).
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
from provenance import append_manifest, build_record, now_utc_iso  # noqa: E402

BIN = Path(__file__).parent


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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--study-dir", required=True, help="e.g. studies/fungi/pezizo_set1")
    ap.add_argument("--uniprot-cache", default="data/uniprot")
    ap.add_argument("--ncbi-cache", default="data/ncbi")
    ap.add_argument("--skip-fetch", action="store_true", help="Rebuild config/data_dir from an already-populated cache, no downloads")
    args = ap.parse_args()

    study_dir = Path(args.study_dir)
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        sys.exit(f"ERROR: {species_csv} not found")

    link_dir = study_dir / "data_dir"
    pep_dir, dna_dir, gff3_dir = (link_dir / d for d in ("pep", "dna", "gff3"))

    config_rows = []
    manifest_records = []
    seen_stems: dict[str, str] = {}  # stem -> Short, to catch a copy-paste duplicate row

    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            short = row["Short"]
            upid = row["UniProt_Proteome_ID"]
            taxid = row["Taxon_ID"]
            gca = row["GCA_Accession"]

            if not args.skip_fetch:
                run([
                    sys.executable, str(BIN / "fetch_uniprot_proteome.py"),
                    "--proteome-id", upid, "--taxid", taxid,
                    "--outdir", args.uniprot_cache, "--short", short,
                ])
                run([
                    sys.executable, str(BIN / "fetch_genome_assembly.py"),
                    "--accession", gca, "--outdir", args.ncbi_cache, "--short", short,
                ])

            # locate cached files
            stem = f"{upid}_{taxid}"
            if stem in seen_stems:
                sys.exit(
                    f"ERROR: {species_csv} lists proteome {stem} twice "
                    f"(Short={seen_stems[stem]!r} and Short={short!r}) -- likely a "
                    f"copy-paste mistake, since each row should be a distinct strain/proteome"
                )
            seen_stems[stem] = short
            fasta_gz = Path(args.uniprot_cache) / upid / f"{stem}.fasta.gz"
            dat_gz = Path(args.uniprot_cache) / upid / f"{stem}.dat.gz"
            genome_dir = Path(args.ncbi_cache) / gca / "extracted" / "ncbi_dataset" / "data" / gca
            fna_candidates = sorted(genome_dir.glob("*.fna")) if genome_dir.exists() else []
            gff_candidates = sorted(genome_dir.glob("*.gff")) if genome_dir.exists() else []
            if not fasta_gz.exists():
                sys.exit(f"[{short}] ERROR: expected {fasta_gz} not found -- run without --skip-fetch first")
            if not fna_candidates:
                sys.exit(f"[{short}] ERROR: expected genome under {genome_dir} not found -- run without --skip-fetch first")

            # per-species GO/Pfam/InterPro/gene-name/description extract (report-facing
            # annotation, distinct from pep_out/dna_out/gff3_out above)
            if dat_gz.exists():
                annotations_dir = study_dir / "annotations"
                annotations_dir.mkdir(parents=True, exist_ok=True)
                run([
                    sys.executable, str(BIN / "extract_dat_annotations.py"),
                    "--dat-gz", str(dat_gz),
                    "--output", str(annotations_dir / f"{upid}.tsv"),
                ])

            # materialize plain files into data_dir/{pep,dna,gff3}/, named by the
            # UniProt proteome stem (see module docstring point 3), not Short
            pep_out = pep_dir / f"{stem}.pep.fa"
            dna_out = dna_dir / f"{stem}.dna.fa"
            gunzip_to(fasta_gz, pep_out)
            dna_dir.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(fna_candidates[0], dna_out)
            gff3_out = ""
            if gff_candidates:
                gff3_dir.mkdir(parents=True, exist_ok=True)
                gff3_out = f"{stem}.gff3"
                shutil.copyfile(gff_candidates[0], gff3_dir / gff3_out)

            config_rows.append({
                "GROUP": row["Group"],
                "Species": row["Species"],
                "Strain": row["Strain"],
                "Protein": pep_out.name,
                "DNA": dna_out.name,
                "GFF3": gff3_out,
                "Short": short,
                "TaxonGroup": row["TaxonGroup"],
            })

            # fold this species' individual provenance sidecars into the study manifest
            for sidecar in (fasta_gz, dat_gz):
                sc = sidecar.with_suffix(sidecar.suffix + ".provenance.yaml")
                if sc.exists():
                    manifest_records.append(load_provenance(sc))
            for f in (fna_candidates[:1] + gff_candidates[:1]):
                sc = f.with_suffix(f.suffix + ".provenance.yaml")
                if sc.exists():
                    manifest_records.append(load_provenance(sc))

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
