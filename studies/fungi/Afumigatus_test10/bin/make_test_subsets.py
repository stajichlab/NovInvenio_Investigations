#!/usr/bin/env python3
"""Build the two small A. fumigatus pangenome test sets from Afumigatus_pangenome.

Set 1, studies/fungi/Afumigatus_test10 (ingroup only, 10 genomes, all IN):
  Asfu_Af293 and Asfu_A1163 (reference strains), the best assembly from each
  Barber cluster 1-7, and a second genome from Barber cluster 1 (the most
  distinct population: mean mash distance 0.0068-0.0075 to every other
  cluster, 0.0009 within).
Set 2, studies/fungi/Afumigatus_test_c1c7 (two populations as groups):
  the best 5 of Barber cluster 1 as IN, the best 5 of Barber cluster 7 as OUT
  (the most distant cluster pair, mean mash 0.0075).

"Best assembly" = fewest contigs among genomes with BUSCO Complete >= 98.9%
(busco_completeness_summary.tsv). BUSCO alone does not separate these
genomes (98.8-99.2%); contig count does (19 to ~1000), and fragmented
assemblies inflate singleton families.

No data is copied: each study's pangenome_runs.yaml points data_dir at
Afumigatus_pangenome/data_dir. Writes config.csv and DATA_MANIFEST.yaml in
each study.
"""
import csv
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import build_record, append_manifest, sha256_of  # noqa: E402

SRC = NII_ROOT / "studies/fungi/Afumigatus_pangenome"
MIN_BUSCO = 98.9
REFERENCES = ["Asfu_Af293", "Asfu_A1163"]


def ranked(cluster, cfg, busco):
    members = [s for s, r in cfg.items() if r["TaxonGroup"] == cluster
               and s in busco and float(busco[s]["Complete_pct"]) >= MIN_BUSCO]
    return sorted(members, key=lambda s: (int(busco[s]["Contigs"]), s))


def write_set(study: Path, rows: list[dict], fields: list[str], rule: str) -> None:
    out = study / "config.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    rec = build_record(
        source_url=f"file://{SRC / 'config.csv'}",
        source_release=f"Afumigatus_pangenome/config.csv sha256 {sha256_of(SRC / 'config.csv')}",
        license="Subset of Afumigatus_pangenome inputs; terms per that study's DATA_MANIFEST.yaml",
        local_path=out.relative_to(NII_ROOT),
        checksum=sha256_of(out),
        derived_by=f"studies/fungi/Afumigatus_test10/bin/make_test_subsets.py -- {rule}",
        extra={"n_genomes": len(rows),
               "busco_source": str((SRC / "busco_completeness_summary.tsv").relative_to(NII_ROOT)),
               "busco_sha256": sha256_of(SRC / "busco_completeness_summary.tsv")},
    )
    append_manifest([rec], study / "DATA_MANIFEST.yaml")
    print(f"{out}: {len(rows)} genomes", file=sys.stderr)


def main() -> int:
    with open(SRC / "config.csv", newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames
        cfg = {r["Short"]: r for r in reader}
    with open(SRC / "busco_completeness_summary.tsv", newline="") as fh:
        busco = {r["Short"]: r for r in csv.DictReader(fh, delimiter="\t")}
    clusters = [f"Barber_cluster{i}" for i in range(1, 8)]

    picks = list(REFERENCES)
    for c in clusters:
        picks.append(ranked(c, cfg, busco)[0])
    picks.append(ranked("Barber_cluster1", cfg, busco)[1])
    rows = [dict(cfg[s], GROUP="IN") for s in picks]
    write_set(NII_ROOT / "studies/fungi/Afumigatus_test10", rows, fields,
              "Af293 + A1163 + fewest-contig genome (BUSCO>=98.9) per Barber cluster 1-7 "
              "+ second of cluster 1; all GROUP=IN")

    rows = ([dict(cfg[s], GROUP="IN") for s in ranked("Barber_cluster1", cfg, busco)[:5]]
            + [dict(cfg[s], GROUP="OUT") for s in ranked("Barber_cluster7", cfg, busco)[:5]])
    write_set(NII_ROOT / "studies/fungi/Afumigatus_test_c1c7", rows, fields,
              "5 fewest-contig genomes (BUSCO>=98.9) of Barber cluster 1 as IN, "
              "5 of Barber cluster 7 as OUT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
