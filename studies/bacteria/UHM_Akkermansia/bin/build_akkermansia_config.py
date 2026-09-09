#!/usr/bin/env python3
"""Build studies/bacteria/UHM_Akkermansia's config.csv + data_dir/ from its species.csv.

Source data: Leila Shadmani's ch3-chitin-evolution project (--staging-dir, default
/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novinvenio_set) -- 6 candidate-
genus ingroup MAGs + 8 Akkermansia-clade outgroup representative genomes, each a
metashot-binned assembly with prodigal gene calling (see --staging-dir's own
groups.tsv for the per-proteome cluster/clade and note columns this script folds into
DATA_MANIFEST.yaml).

Proteins are read from --staging-dir's ingroup_prefixed/outgroup_prefixed subdirs
(not the original ingroup/outgroup dirs): Leila re-wrote every header there to lead
with its own source genome id (e.g. ">EHM034720|EHA02359_bin.2^_0_1",
">UHM1073.23039_R.bin.125|k141_222743_1") specifically to fix bare megahit contig ids
colliding across proteomes -- UHM_Koxytoca's same problem, fixed upstream this time
instead of by a header-prefixing step in this script. One header artifact is still
fixed here on copy: a subset of the ingroup proteomes (the ones whose bin id is
followed by prodigal coordinates) carry a "^_" artifact instead of a real separator
(">EHM034720|EHA02359_bin.2^_0_1"), rewritten to "_gene_"
(">EHM034720|EHA02359_bin.2_gene_0_1"). A no-op when the artifact is absent (every
other proteome).

DNA: --staging-dir's genome_map.tsv (new 2026-09-07) maps every proteome to its
source_genome and a local nucleotide_fasta -- present for 11 of 14 rows (the 6
ingroup MAGs + outgroup clades 01/02/04/05/06); the other 3 (clade03/07/08) are
listed there as "MISSING" because their source genome is instead an NCBI WGS
assembly, already resolved via species.csv's NCBI_Accession column (see
bin/fetch_akkermansia_outgroup_dna.sh) and fetched into --ncbi-cache. Rows with an
NCBI_Accession keep using that NCBI-fetched genome (already committed, already
checksummed) even where genome_map.tsv also has a redundant local copy (clade02,
clade04); genome_map.tsv's nucleotide_fasta is only used for rows with no
NCBI_Accession. Net effect: all 14 rows now get DNA, so --run_tool pairwise TBLASTN
validate/loss-search workflows are runnable against the full set.

GFF3: --gff-dir (new 2026-09-01, Leila Shadmani's pangenome/gff199/) holds one
prodigal-called <source_genome>.gff per genome, keyed by genome_map.tsv's
source_genome column (e.g. NOVEL_sp01_n18.faa -> EHM034720 -> EHM034720.gff) --
not species.csv's own Accession column, which is a proteome/study id
("NOVEL_sp01_n18"), not a genome accession. Only rows with a genome_map.tsv
source_genome entry that also has a match in this 199-genome set get a GFF3 (the 6
ingroup rows + outgroup clades 01/05/06); rows with no genome_map.tsv entry, or an
NCBI-accession outgroup genome not represented in gff199 (clade02/03/04/07/08), are
left with an empty GFF3 config.csv cell (no error). The same "^_" prodigal-header
artifact fixed in pep.fa headers (module docstring above) also appears in these
GFFs' seqid column and "# Sequence Data" seqhdr= comment (only in the EHM034* rows),
rewritten identically to "_gene_" so contig ids stay consistent between a row's
.pep.fa and .gff3.

This is a study-specific script (lives under this study's own bin/, not NII's shared
bin/ -- see CLAUDE.md's "Where new code goes"), so it hardcodes NII_ROOT below rather
than deriving paths from __file__/cwd.

Input: <study-dir>/species.csv, columns:
    Short,Species,Strain,Group,TaxonGroup,Source,Accession,NCBI_Accession
  Source is "ingroup" or "outgroup" (selects the --staging-dir *_prefixed subdirectory);
  Accession is the original proteome's filename stem (<Accession>.faa under that
  subdirectory, and the key into --staging-dir's genome_map.tsv); NCBI_Accession is a
  GCA_* assembly accession (only set for the 3 outgroup clades with one -- see module
  docstring) already fetched into --ncbi-cache by bin/fetch_akkermansia_outgroup_dna.sh.

Writes:
  - <study-dir>/data_dir/pep/<Short>.pep.fa (copy of the source *_prefixed .faa,
    with the '^_' artifact fix applied)
  - <study-dir>/data_dir/dna/<Short>.dna.fa (copied from --ncbi-cache for rows with an
    NCBI_Accession, else from genome_map.tsv's nucleotide_fasta when present)
  - <study-dir>/data_dir/gff/<Short>.gff3 (copy of --gff-dir's <Accession>.gff, for the
    6 ingroup rows that have one; '^_' artifact fix applied, see module docstring)
  - <study-dir>/config.csv (GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup)
  - <study-dir>/DATA_MANIFEST.yaml (provenance: Leila Shadmani, ch3-chitin-evolution,
    metashot assembly + prodigal gene calling, frozen for this study 2026-09-01, for
    the proteins, GFFs, and the genome_map.tsv-sourced DNA; NCBI Datasets provenance
    for the 3 NCBI-accession outgroup genomes' DNA)
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import yaml

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import append_manifest, build_record, sha256_of  # noqa: E402

LICENSE = "Internal / unpublished (Stajich lab, Leila Shadmani ch3-chitin-evolution project)"
SOURCE_RELEASE = (
    "Leila Shadmani ch3-chitin-evolution MAG/genome set (metashot binning + prodigal "
    "gene calling), frozen for this study 2026-09-01"
)
NCBI_LICENSE = "Public Domain (NCBI, https://www.ncbi.nlm.nih.gov/home/about/policies/)"

GFF_SOURCE_RELEASE = (
    "Leila Shadmani ch3-chitin-evolution pangenome/gff199 prodigal GFF3 set, "
    "shared for this study 2026-09-01"
)

SOURCE_SUBDIR = {"ingroup": "ingroup_prefixed", "outgroup": "outgroup_prefixed"}

CARET_ARTIFACT_RE = re.compile(r"\^_")


def copy_with_header_fix(src: Path, dst: Path) -> None:
    """Copy src to dst, rewriting any '^_' prodigal-header artifact to '_gene_'
    on header lines only (see module docstring). A no-op on non-header lines
    and on headers without the artifact."""
    with open(src) as fh_in, open(dst, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                line = CARET_ARTIFACT_RE.sub("_gene_", line, count=1)
            fh_out.write(line)


def copy_gff_with_artifact_fix(src: Path, dst: Path) -> None:
    """Copy a prodigal GFF3 src to dst, rewriting the same '^_' artifact (see
    module docstring) wherever it appears in the seqid column (col 1) or the
    '# Sequence Data: ...;seqhdr="..."' comment line. No-op elsewhere/when absent."""
    with open(src) as fh_in, open(dst, "w") as fh_out:
        for line in fh_in:
            if line.startswith("# Sequence Data:"):
                line = CARET_ARTIFACT_RE.sub("_gene_", line)
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if fields and fields[0]:
                    fields[0] = CARET_ARTIFACT_RE.sub("_gene_", fields[0])
                line = "\t".join(fields) + "\n"
            fh_out.write(line)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study-dir", default=str(NII_ROOT / "studies/bacteria/UHM_Akkermansia"))
    ap.add_argument("--staging-dir", default="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novinvenio_set")
    ap.add_argument("--gff-dir", default="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/gff199")
    ap.add_argument("--ncbi-cache", default=str(NII_ROOT / "data/ncbi"))
    args = ap.parse_args()

    study_dir = Path(args.study_dir)
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        sys.exit(f"ERROR: {species_csv} not found")

    staging_dir = Path(args.staging_dir)
    gff_dir = Path(args.gff_dir)
    ncbi_cache = Path(args.ncbi_cache)
    groups_tsv = staging_dir / "groups.tsv"
    groups_by_proteome = {}
    if groups_tsv.exists():
        with open(groups_tsv, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                groups_by_proteome[row["proteome"]] = row

    genome_map_tsv = staging_dir / "genome_map.tsv"
    genome_map_by_proteome = {}
    if genome_map_tsv.exists():
        with open(genome_map_tsv, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                genome_map_by_proteome[row["proteome"]] = row

    pep_dir = study_dir / "data_dir" / "pep"
    dna_dir = study_dir / "data_dir" / "dna"
    pep_dir.mkdir(parents=True, exist_ok=True)

    config_rows = []
    manifest_records = []

    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            short = row["Short"]
            source = row["Source"]
            accession = row["Accession"]

            if source not in ("ingroup", "outgroup"):
                sys.exit(f"[{short}] ERROR: unknown Source {source!r} (expected ingroup/outgroup)")

            faa_src = staging_dir / SOURCE_SUBDIR[source] / f"{accession}.faa"
            if not faa_src.exists():
                sys.exit(f"[{short}] ERROR: expected protein FASTA {faa_src} not found")

            pep_out = pep_dir / f"{short}.pep.fa"
            copy_with_header_fix(faa_src, pep_out)

            g = groups_by_proteome.get(f"{accession}.faa", {})
            gm = genome_map_by_proteome.get(f"{accession}.faa", {})
            manifest_records.append(build_record(
                source_url=f"(internal -- UCR HPCC, {faa_src})",
                source_release=SOURCE_RELEASE,
                license=LICENSE,
                local_path=pep_out,
                checksum=sha256_of(pep_out),
                derived_by=(
                    f"protein.faa copied from {faa_src} via "
                    "studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py: "
                    "headers already prefixed with their source genome id upstream "
                    "(Leila Shadmani's ingroup_prefixed/outgroup_prefixed, to avoid "
                    "bare-contig-id collisions across proteomes); any remaining "
                    "'<bin_id>^_<contig>_<gene>' prodigal header artifact rewritten "
                    "to '<bin_id>_gene_<contig>_<gene>' (no-op when absent)"
                ),
                extra={
                    "cluster_or_clade": g.get("cluster_or_clade", ""),
                    "source_genome": gm.get("source_genome", g.get("source_genome", "")),
                    "n_proteins": g.get("n_proteins", ""),
                    "note": g.get("note", ""),
                },
            ))

            gff_name = ""
            source_genome_id = gm.get("source_genome", "").strip()
            gff_src = gff_dir / f"{source_genome_id}.gff" if source_genome_id else None
            if gff_src is not None and gff_src.exists():
                gff_out_dir = study_dir / "data_dir" / "gff"
                gff_out_dir.mkdir(parents=True, exist_ok=True)
                gff_out = gff_out_dir / f"{short}.gff3"
                copy_gff_with_artifact_fix(gff_src, gff_out)
                gff_name = gff_out.name

                manifest_records.append(build_record(
                    source_url=f"(internal -- UCR HPCC, {gff_src})",
                    source_release=GFF_SOURCE_RELEASE,
                    license=LICENSE,
                    local_path=gff_out,
                    checksum=sha256_of(gff_out),
                    derived_by=(
                        f"GFF3 copied from {gff_src} via "
                        "studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py: "
                        "prodigal gene calls, same '^_' contig-id artifact fix as this "
                        "row's .pep.fa applied to the seqid column and seqhdr= comment "
                        "(no-op when absent)"
                    ),
                ))

            dna_name = ""
            ncbi_accession = row.get("NCBI_Accession", "").strip()
            if ncbi_accession:
                genome_dir = ncbi_cache / ncbi_accession / "extracted" / "ncbi_dataset" / "data" / ncbi_accession
                fna_candidates = sorted(genome_dir.glob("*.fna")) if genome_dir.exists() else []
                if not fna_candidates:
                    sys.exit(f"[{short}] ERROR: expected genome under {genome_dir} -- run studies/bacteria/UHM_Akkermansia/bin/fetch_akkermansia_outgroup_dna.sh first")
                dna_dir.mkdir(parents=True, exist_ok=True)
                dna_out = dna_dir / f"{short}.dna.fa"
                dna_out.write_bytes(fna_candidates[0].read_bytes())
                dna_name = dna_out.name

                for f in fna_candidates[:1]:
                    sc = f.with_suffix(f.suffix + ".provenance.yaml")
                    if sc.exists():
                        with open(sc) as fh_sc:
                            manifest_records.append(yaml.safe_load(fh_sc))

                manifest_records.append(build_record(
                    source_url=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{ncbi_accession}/",
                    source_release=f"NCBI genome assembly {ncbi_accession}",
                    license=NCBI_LICENSE,
                    local_path=dna_out,
                    checksum=sha256_of(dna_out),
                    derived_by=(
                        f"genomic.fna from the {ncbi_accession} NCBI Datasets genome package "
                        f"(resolved from proteome {accession}.faa's own WGS project prefix -- "
                        "see bin/fetch_akkermansia_outgroup_dna.sh), copied via "
                        "studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py"
                    ),
                ))
            else:
                nucleotide_fasta = gm.get("nucleotide_fasta", "").strip()
                if nucleotide_fasta and nucleotide_fasta != "MISSING":
                    dna_src = Path(nucleotide_fasta)
                    if not dna_src.exists():
                        sys.exit(f"[{short}] ERROR: genome_map.tsv nucleotide_fasta {dna_src} not found")
                    dna_dir.mkdir(parents=True, exist_ok=True)
                    dna_out = dna_dir / f"{short}.dna.fa"
                    dna_out.write_bytes(dna_src.read_bytes())
                    dna_name = dna_out.name

                    manifest_records.append(build_record(
                        source_url=f"(internal -- UCR HPCC, {dna_src})",
                        source_release=SOURCE_RELEASE,
                        license=LICENSE,
                        local_path=dna_out,
                        checksum=sha256_of(dna_out),
                        derived_by=(
                            f"genome FASTA copied from {dna_src} via "
                            "studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py "
                            f"(genome_map.tsv: proteome {accession}.faa -> source_genome "
                            f"{gm.get('source_genome', '')!r})"
                        ),
                    ))

            config_rows.append({
                "GROUP": row["Group"],
                "Species": row["Species"],
                "Strain": row["Strain"],
                "Protein": pep_out.name,
                "DNA": dna_name,
                "GFF3": gff_name,
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
        derived_by=f"studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py --study-dir {study_dir}",
    ))
    append_manifest(manifest_records, study_dir / "DATA_MANIFEST.yaml")

    print(f"\nWrote {config_csv} ({len(config_rows)} species)", file=sys.stderr)
    print(f"Wrote {pep_dir} (pass as --data_dir to nf_NovInvenio)", file=sys.stderr)
    print(f"Wrote {study_dir / 'DATA_MANIFEST.yaml'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
