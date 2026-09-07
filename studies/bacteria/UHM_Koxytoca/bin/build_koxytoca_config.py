#!/usr/bin/env python3
"""Build studies/bacteria/UHM_Koxytoca's config.csv + data_dir/ from its species.csv.

Unlike NII's shared bin/build_study_config.py (which pulls both protein and DNA per
species from UniProt/NCBI keyed off a UniProt proteome ID), this study's proteomes
already exist as local .faa files -- 16 NCBI RefSeq outgroup references and 18 UHM
metagenome-assembled (metashot binning + prodigal gene calling) ingroup genomes. This
script materializes those directly instead of fetching them, and only fetches DNA for
the outgroup (the only group nf_NovInvenio's VALIDATE workflow ever runs TBLASTN
against when cluster_tool=pairwise -- see fetch_koxytoca_outgroup_dna.sh's header).

This is a study-specific script (lives under this study's own bin/, not NII's shared
bin/ -- see CLAUDE.md's "Where new code goes"), so it hardcodes NII_ROOT below rather
than deriving paths from __file__/cwd.

Input: <study-dir>/species.csv, columns:
    Short,Species,Strain,Group,TaxonGroup,Source,Accession
  Source is "ncbi_refseq" (Accession = GCF_* assembly accession, matches the outgroup
  .faa/ncbi cache filenames) or "uhm_mag_metashot" (Accession = original bin filename
  stem, matches the ingroup .faa filenames).

Requires fetch_koxytoca_outgroup_dna.sh to have been run first (populates
--ncbi-cache with each outgroup accession's genome + provenance sidecar).

Writes:
  - <study-dir>/data_dir/pep/<Short>.pep.fa  (IN) or <Accession>.faa (OUT, same
    filename as the original koxytoca_outgroup_faa/<Accession>.faa -- keeps a
    file identifiable by its GCF_* accession instead of the display-only Short
    code)
  - <study-dir>/data_dir/dna/<Accession>.fna  (OUT only, copied from --ncbi-cache,
    same GCF_* naming as the protein file above)
  - <study-dir>/config.csv (GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup)
  - <study-dir>/DATA_MANIFEST.yaml (provenance: NCBI fetch sidecars for outgroup DNA,
    plus a synthesized record for each copied .faa -- NCBI Datasets provenance for
    outgroup proteins, and the metashot/prodigal + UHM_combined_MAG_set attribution
    for ingroup proteins, since neither is fetched by a pull script of its own)
"""
import argparse
import csv
import shutil
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import append_manifest, build_record, now_utc_iso, sha256_of  # noqa: E402

NCBI_LICENSE = "Public Domain (NCBI, https://www.ncbi.nlm.nih.gov/home/about/policies/)"
UHM_LICENSE = "Internal / unpublished (Stajich lab, UHM collaboration)"


def load_provenance(sidecar: Path) -> dict:
    import yaml
    with open(sidecar) as fh:
        return yaml.safe_load(fh)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study-dir", default=str(NII_ROOT / "studies/bacteria/UHM_Koxytoca"))
    ap.add_argument("--outgroup-faa-dir", default="/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_outgroup_faa")
    ap.add_argument("--ingroup-faa-dir", default="/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_ingroup_faa")
    ap.add_argument("--ncbi-cache", default=str(NII_ROOT / "data/ncbi"))
    args = ap.parse_args()

    study_dir = Path(args.study_dir)
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        sys.exit(f"ERROR: {species_csv} not found")

    outgroup_faa_dir = Path(args.outgroup_faa_dir)
    ingroup_faa_dir = Path(args.ingroup_faa_dir)
    ncbi_cache = Path(args.ncbi_cache)

    link_dir = study_dir / "data_dir"
    pep_dir, dna_dir = link_dir / "pep", link_dir / "dna"
    pep_dir.mkdir(parents=True, exist_ok=True)

    config_rows = []
    manifest_records = []

    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            short = row["Short"]
            group = row["Group"]
            source = row["Source"]
            accession = row["Accession"]

            if source == "ncbi_refseq":
                faa_src = outgroup_faa_dir / f"{accession}.faa"
            elif source == "uhm_mag_metashot":
                faa_src = ingroup_faa_dir / f"{accession}.faa"
            else:
                sys.exit(f"[{short}] ERROR: unknown Source {source!r}")

            if not faa_src.exists():
                sys.exit(f"[{short}] ERROR: expected protein FASTA {faa_src} not found")

            pep_out = pep_dir / (f"{accession}.faa" if group == "OUT" else f"{short}.pep.fa")
            shutil.copyfile(faa_src, pep_out)

            dna_out_name = ""
            if group == "OUT":
                genome_dir = ncbi_cache / accession / "extracted" / "ncbi_dataset" / "data" / accession
                fna_candidates = sorted(genome_dir.glob("*.fna")) if genome_dir.exists() else []
                if not fna_candidates:
                    sys.exit(f"[{short}] ERROR: expected genome under {genome_dir} -- run studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh first")
                dna_dir.mkdir(parents=True, exist_ok=True)
                dna_out = dna_dir / f"{accession}.fna"
                shutil.copyfile(fna_candidates[0], dna_out)
                dna_out_name = dna_out.name

                for f in fna_candidates[:1]:
                    sc = f.with_suffix(f.suffix + ".provenance.yaml")
                    if sc.exists():
                        manifest_records.append(load_provenance(sc))

                manifest_records.append(build_record(
                    source_url=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/",
                    source_release=f"NCBI RefSeq assembly {accession}",
                    license=NCBI_LICENSE,
                    local_path=pep_out,
                    checksum=sha256_of(pep_out),
                    derived_by=f"protein.faa from the {accession} NCBI Datasets genome package (jpere468/klebsiella_story/koxytoca_outgroup_faa), copied via studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py",
                ))
            else:
                manifest_records.append(build_record(
                    source_url="(internal -- UCR HPCC, /bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_ingroup_faa)",
                    source_release="UHM combined MAG set (binned/annotated 2025-2026, frozen for this study 2026-09-01)",
                    license=UHM_LICENSE,
                    local_path=pep_out,
                    checksum=sha256_of(pep_out),
                    derived_by=(
                        f"MAG {accession}: metashot binning + prodigal gene calling "
                        "(jpere468 klebsiella_story project); copied via studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py"
                    ),
                ))

            config_rows.append({
                "GROUP": group,
                "Species": row["Species"],
                "Strain": row["Strain"],
                "Protein": pep_out.name,
                "DNA": dna_out_name,
                "GFF3": "",
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
        derived_by=f"studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py --study-dir {study_dir}",
    ))
    append_manifest(manifest_records, study_dir / "DATA_MANIFEST.yaml")

    print(f"\nWrote {config_csv} ({len(config_rows)} species)", file=sys.stderr)
    print(f"Wrote {link_dir} (pep/dna -- pass as --data_dir to nf_NovInvenio)", file=sys.stderr)
    print(f"Wrote {study_dir / 'DATA_MANIFEST.yaml'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
