#!/usr/bin/env python3
"""Build studies/bacteria/UHM_Akkermansia's config.csv + data_dir/ from its species.csv.

Source data: Leila Shadmani's ch3-chitin-evolution project (--staging-dir, default
/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novinvenio_set) -- 6 candidate-
genus ingroup MAGs + 8 Akkermansia-clade outgroup representative genomes, each a
metashot-binned assembly with prodigal gene calling (see --staging-dir's own
groups.tsv for the per-proteome cluster/clade and note columns this script folds into
DATA_MANIFEST.yaml). Protein-only: no assembled-genome DNA was staged for either group,
so config.csv's DNA column is empty for every row -- --run_tool pairwise TBLASTN
validate/loss-search workflows aren't runnable for this study until DNA is supplied.

Header fixes applied on copy:
  1. A subset of the ingroup proteomes (the ones whose prodigal headers start with a
     MAG bin id, e.g. "EHA02359_bin.2") carry a "^_" artifact between the bin id and
     the contig/gene coordinates (">EHA02359_bin.2^_257_4") instead of a real
     separator. Rewritten to "_gene_" (">EHA02359_bin.2_gene_257_4"). A no-op when
     the artifact is absent (every other proteome).
  2. Every header's id is then prefixed "<Short>__", uniformly across all 14
     proteomes -- bare megahit contig ids with no embedded genome identity
     (">k141_<contig>_<gene>", 3 ingroup + 3 outgroup proteomes) are exactly
     UHM_Koxytoca's colliding-MAG-contig problem waiting to happen (see that
     study's bin/prefix_uhm_bin_headers.py) once more proteomes are added, but the
     bin-id and NCBI WGS-accession headers get the same "<Short>__" prefix too, so
     every id in data_dir/pep/ is traceable to its proteome by a consistent,
     single naming convention rather than three different ones. Idempotent: a
     header already starting with its own "<Short>__" is left alone.

This is a study-specific script (lives under this study's own bin/, not NII's shared
bin/ -- see CLAUDE.md's "Where new code goes"), so it hardcodes NII_ROOT below rather
than deriving paths from __file__/cwd.

Input: <study-dir>/species.csv, columns:
    Short,Species,Strain,Group,TaxonGroup,Source,Accession
  Source is "ingroup" or "outgroup" (selects the --staging-dir subdirectory);
  Accession is the original proteome's filename stem (<Accession>.faa under that
  subdirectory).

Writes:
  - <study-dir>/data_dir/pep/<Short>.pep.fa (header-fixed copy of the source .faa)
  - <study-dir>/config.csv (GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup;
    DNA always empty -- see module docstring)
  - <study-dir>/DATA_MANIFEST.yaml (provenance: Leila Shadmani, ch3-chitin-evolution,
    metashot assembly + prodigal gene calling, frozen for this study 2026-09-01)
"""
import argparse
import csv
import re
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import append_manifest, build_record, sha256_of  # noqa: E402

LICENSE = "Internal / unpublished (Stajich lab, Leila Shadmani ch3-chitin-evolution project)"
SOURCE_RELEASE = (
    "Leila Shadmani ch3-chitin-evolution MAG/genome set (metashot binning + prodigal "
    "gene calling), frozen for this study 2026-09-01"
)

CARET_ARTIFACT_RE = re.compile(r"\^_")
HEADER_ID_RE = re.compile(r"^>(\S+)")


def fix_header_line(line: str, short: str) -> str:
    """Apply both header fixes (see module docstring) to one FASTA line: the
    '^_' artifact rewrite, then a uniform '<Short>__' id prefix (idempotent --
    a no-op if the id already starts with it). A no-op entirely on non-header
    lines."""
    if not line.startswith(">"):
        return line
    line = CARET_ARTIFACT_RE.sub("_gene_", line, count=1)
    if line.startswith(f">{short}__"):
        return line
    return HEADER_ID_RE.sub(rf">{short}__\1", line, count=1)


def copy_with_header_fix(src: Path, dst: Path, short: str) -> None:
    with open(src) as fh_in, open(dst, "w") as fh_out:
        for line in fh_in:
            fh_out.write(fix_header_line(line, short))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study-dir", default=str(NII_ROOT / "studies/bacteria/UHM_Akkermansia"))
    ap.add_argument("--staging-dir", default="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novinvenio_set")
    args = ap.parse_args()

    study_dir = Path(args.study_dir)
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        sys.exit(f"ERROR: {species_csv} not found")

    staging_dir = Path(args.staging_dir)
    groups_tsv = staging_dir / "groups.tsv"
    groups_by_proteome = {}
    if groups_tsv.exists():
        with open(groups_tsv, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                groups_by_proteome[row["proteome"]] = row

    pep_dir = study_dir / "data_dir" / "pep"
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

            faa_src = staging_dir / source / f"{accession}.faa"
            if not faa_src.exists():
                sys.exit(f"[{short}] ERROR: expected protein FASTA {faa_src} not found")

            pep_out = pep_dir / f"{short}.pep.fa"
            copy_with_header_fix(faa_src, pep_out, short)

            g = groups_by_proteome.get(f"{accession}.faa", {})
            manifest_records.append(build_record(
                source_url=f"(internal -- UCR HPCC, {faa_src})",
                source_release=SOURCE_RELEASE,
                license=LICENSE,
                local_path=pep_out,
                checksum=sha256_of(pep_out),
                derived_by=(
                    f"protein.faa copied from {faa_src} via "
                    "studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py: "
                    "any '<bin_id>^_<contig>_<gene>' prodigal header artifact rewritten "
                    "to '<bin_id>_gene_<contig>_<gene>' (no-op when absent), then every "
                    f"header id uniformly prefixed '{short}__'"
                ),
                extra={
                    "cluster_or_clade": g.get("cluster_or_clade", ""),
                    "source_genome": g.get("source_genome", ""),
                    "n_proteins": g.get("n_proteins", ""),
                    "note": g.get("note", ""),
                },
            ))

            config_rows.append({
                "GROUP": row["Group"],
                "Species": row["Species"],
                "Strain": row["Strain"],
                "Protein": pep_out.name,
                "DNA": "",
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
        derived_by=f"studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py --study-dir {study_dir}",
    ))
    append_manifest(manifest_records, study_dir / "DATA_MANIFEST.yaml")

    print(f"\nWrote {config_csv} ({len(config_rows)} species)", file=sys.stderr)
    print(f"Wrote {pep_dir} (pass as --data_dir to nf_NovInvenio)", file=sys.stderr)
    print(f"Wrote {study_dir / 'DATA_MANIFEST.yaml'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
