#!/usr/bin/env python3
"""Write PROVENANCE.yaml for a nr taxonomy screen (lichen_taxcheck.sh or candidate_taxcheck.sh).

Usage: write_screen_provenance.py <screen dir> <script> <input> [<input> ...]
Records the SHA-256 of hits.tsv.zst, queries.txt, the summary, the screen script, the summary
script and every input FASTA, plus the diamond version and database listing that the job wrote
to diamond_version.txt. hits.tsv.zst and queries.txt are not tracked (hits.tsv.zst is 2 MB per
screen); the tracked files are summary.tsv, best_hits.tsv.gz (candidate screen) and this record.
Absolute paths on purpose (study-specific script, NII CLAUDE.md).
"""
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT))
import yaml  # noqa: E402
from lib.provenance import build_record, write_record, sha256_of  # noqa: E402

SUMMARY_SCRIPT = NII_ROOT / "studies/fungi/yeast_convergence/bin/lichen_taxcheck_summary.py"


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    d, script, inputs = Path(sys.argv[1]), Path(sys.argv[2]), [Path(p) for p in sys.argv[3:]]
    for need in ("hits.tsv.zst", "queries.txt", "summary.tsv", "diamond_version.txt"):
        if not (d / need).is_file():
            raise SystemExit(f"ERROR: {d / need} missing")
    ver = (d / "diamond_version.txt").read_text().splitlines()
    rec = build_record(
        source_url="derived: diamond blastp against NCBI nr (clustered, taxonomy built in); no download",
        source_release="NCBI nr 20260128 (nr_cluster_seq.dmnd); taxonomy nodes.dmp 20260128 under "
                       "/srv/projects/db/ncbi/preformatted/20260128",
        license="derived table; inputs carry their own licenses",
        local_path=d / "summary.tsv",
        derived_by=f"sbatch {script.name} (diamond blastp --fast --max-target-seqs 5 --evalue 1e-5 -b 4 -c 1, 32 threads), "
                   f"then python3 {SUMMARY_SCRIPT.name} {d} (best hit by bitscore, group by NCBI lineage of the first taxon id)",
        extra={"screen_script": str(script), "screen_script_sha256": sha256_of(script),
               "summary_script_sha256": sha256_of(SUMMARY_SCRIPT),
               "diamond": ver[0] if ver else None, "database": ver[1:] if len(ver) > 1 else None,
               "hits_tsv_zst_sha256": sha256_of(d / "hits.tsv.zst"),
               "queries_txt_sha256": sha256_of(d / "queries.txt"),
               "inputs": {str(p): sha256_of(p) for p in inputs},
               "note": "A screen, not a taxonomic assignment: best hit only, --fast, E-value 1e-5."})
    write_record(rec, d / "PROVENANCE.yaml")
    print(f"wrote {d / 'PROVENANCE.yaml'}")


if __name__ == "__main__":
    main()
