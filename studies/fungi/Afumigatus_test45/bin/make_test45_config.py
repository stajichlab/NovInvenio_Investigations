#!/usr/bin/env python3
"""Build studies/fungi/Afumigatus_test45/config.csv: the 45-genome
A. fumigatus test set for fast iteration on island generation and the island
interface.

Genome choice is the N=45, seed=3 subset of the island genome-count sweep
(Afumigatus_pangenome/analysis/island_genome_count, ISLAND_GENOME_COUNT.md):
45 representative ingroup strains of full_v070 (Af293 forced, 44 drawn with
numpy default_rng(3)), plus the 2 UniProt outgroup references. In the sweep,
this subset gave 75,194 co-occurring pairs and 3,152 islands, hit all 7
high-confidence Af293 Starships, and matched 10 of 13 scorable Starships by
presence pattern. Chosen over seeds 1 and 2 because it includes A1163 and had
the fewest co-occurring pairs (fastest). It has no Barber cluster 6 genome.

Rows (GROUP, TaxonGroup and file names) are copied unchanged from
Afumigatus_pangenome/config.csv. No data is copied: pangenome_runs.yaml points
data_dir at Afumigatus_pangenome/data_dir.
"""
import csv
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import build_record, append_manifest, sha256_of  # noqa: E402

SRC = NII_ROOT / "studies/fungi/Afumigatus_pangenome"
SUBSET = SRC / "analysis/island_genome_count/outputs/runs/n45_s3/subset_strains.txt"
STUDY = NII_ROOT / "studies/fungi/Afumigatus_test45"
EXPECTED = 47  # 45 IN + 2 OUT


def main() -> int:
    keep = SUBSET.read_text().split()
    with open(SRC / "config.csv", newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames
        cfg = {r["Short"]: r for r in reader}
    missing = [s for s in keep if s not in cfg]
    if missing or len(keep) != EXPECTED:
        raise SystemExit(f"subset mismatch: {len(keep)} strains, missing {missing}")

    out = STUDY / "config.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(cfg[s] for s in keep)
    rec = build_record(
        source_url=f"file://{SRC / 'config.csv'}",
        source_release=f"Afumigatus_pangenome/config.csv sha256 {sha256_of(SRC / 'config.csv')}",
        license="Subset of Afumigatus_pangenome inputs; terms per that study's DATA_MANIFEST.yaml",
        local_path=out.relative_to(NII_ROOT),
        checksum=sha256_of(out),
        derived_by="studies/fungi/Afumigatus_test45/bin/make_test45_config.py -- rows of "
                   "Afumigatus_pangenome/config.csv for the island genome-count sweep subset n45_s3 "
                   "(45 full_v070 representative IN strains incl. Af293 + 2 OUT references)",
        extra={"n_genomes": len(keep),
               "subset_source": str(SUBSET.relative_to(NII_ROOT)),
               "subset_sha256": sha256_of(SUBSET)},
    )
    append_manifest([rec], STUDY / "DATA_MANIFEST.yaml")
    print(f"{out}: {len(keep)} genomes", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
