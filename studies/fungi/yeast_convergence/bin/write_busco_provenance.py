#!/usr/bin/env python3
"""Write analysis/busco/summary.tsv.provenance.yaml for the BUSCO table of the 16 contrasts.

Records the SHA-256 of summary.tsv and manifest_all.tsv, of every input proteome FASTA in the
manifest, the BUSCO version and lineage, the jobs, and the rows that were not made by
bin/busco_manifest.py. Run after bin/run_busco_proteomes.sh:
    python3.12 bin/write_busco_provenance.py
Absolute paths on purpose (study-specific script, NII CLAUDE.md).
"""
import csv, re, sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
B = NII_ROOT / "studies/fungi/yeast_convergence/analysis/busco"
LINEAGE = "/bigdata/stajichlab/shared/lib/funannotate_db/fungi_odb10"
HAND_ROWS = ["Microsporum_canis_BMU_08102", "Microsporum_canis_BMU_11032", "Microsporum_canis_CBS_113480",
             "Microsporum_audouinii_IHEM_2214", "Microsporum_ferrugineum_M1a", "Podosphaera_xanthii_YZU573"]
sys.path.insert(0, str(NII_ROOT))
import yaml  # noqa: E402
from lib.provenance import build_record, write_record, sha256_of  # noqa: E402


def main():
    manifest = list(csv.DictReader(open(B / "manifest_all.tsv"), delimiter="\t"))
    inputs = {}
    for r in manifest:
        f = Path(r["fasta"]) if r["fasta"] else None
        inputs[r["name"]] = {"fasta": r["fasta"], "sha256": sha256_of(f) if f and f.is_file() else None}
    versions = set()
    for f in sorted((B / "runs").glob("*.short_summary.txt")):
        m = re.search(r"BUSCO version is: (\S+)", f.read_text())
        if m:
            versions.add(m.group(1))
    rec = build_record(
        source_url="derived: BUSCO on the proteomes listed in manifest_all.tsv (no download)",
        source_release=f"BUSCO {', '.join(sorted(versions)) or 'version not recorded'}; lineage fungi_odb10 "
                       f"(creation 2019-12-13, 549 genomes, 758 BUSCOs) at {LINEAGE}",
        license="derived table; inputs carry their own licenses (UniProt CC BY 4.0; BFD lab annotations unpublished)",
        local_path=B / "summary.tsv",
        derived_by="bin/busco_manifest.py then bin/run_busco_proteomes.sh (busco -m proteins -l fungi_odb10 --offline, "
                   "8 CPUs); this record: bin/write_busco_provenance.py. SLURM jobs 30326454 to 30326456 (106 proteomes) "
                   "and 30335621, 30335622 (6 candidate genomes).",
        extra={"manifest_sha256": sha256_of(B / "manifest_all.tsv"),
               "input_proteomes": inputs,
               "rows_added_by_hand": {
                   "names": HAND_ROWS,
                   "note": "These six manifest rows were written by hand (the candidate genomes were chosen after the "
                           "first BUSCO pass), so their names are directory names and not the names busco_manifest.py "
                           "would write. Five are Microsporum candidates and one is Podosphaera xanthii YZU573."},
               "stale_rows": {
                   "note": "summary.tsv keeps rows of proteomes that species.csv no longer uses (for example "
                           "Micrcani = M. canis BMU_08359, 60.4%, replaced by BMU_08102; P. fusca; two other M. canis "
                           "candidates) so the replacement choices stay traceable. Match rows to species.csv by "
                           "accession (manifest_all.tsv), not by Short."}})
    write_record(rec, B / "summary.tsv.provenance.yaml")
    print(f"wrote {B / 'summary.tsv.provenance.yaml'}; {len(inputs)} inputs, missing FASTA: "
          f"{[k for k, v in inputs.items() if not v['sha256']]}")


if __name__ == "__main__":
    main()
