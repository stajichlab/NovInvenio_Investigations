#!/usr/bin/env python3
"""List every distinct proteome used by the yeast_convergence contrasts.

Writes a TSV (name, source, accession, fasta) for BUSCO runs. local_faa rows use the
path in Protein_Accession. uniprot rows use the cached FASTA under data/uniprot/<UPID>/
(found by glob). ncbi rows use data_dir/pep/<Short>.pep.fa; rows not yet fetched are listed with an empty
fasta column; run bin/run_study.sh (or build_study_config.py) to fetch them, then rerun.
"""
import argparse, csv, glob, sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--uniprot-cache", default=str(NII_ROOT / "data/uniprot"))
    a = ap.parse_args()
    seen = {}
    for f in sorted(glob.glob(str(NII_ROOT / "studies/fungi/yeast_convergence/yc_*/species.csv"))):
        for r in csv.DictReader(open(f)):
            src, acc = r["Protein_Source"], r["Protein_Accession"]
            if not src:
                print(f"WARNING blank row {Path(f).parent.name} {r['Short']}", file=sys.stderr)
                continue
            key = (src, acc)
            if key in seen:
                continue
            fasta = ""
            if src == "local_faa" and Path(acc).exists():
                fasta = acc
            elif src == "uniprot":
                hits = sorted(glob.glob(f"{a.uniprot_cache}/{acc}/*.fasta.gz"))
                fasta = hits[0] if hits else ""
            elif src == "ncbi":
                hits = sorted(glob.glob(str(NII_ROOT / f"studies/fungi/yeast_convergence/yc_*/data_dir/pep/{r['Short']}.pep.fa")))
                fasta = hits[0] if hits else ""
            name = r["Short"] if src == "local_faa" else acc
            seen[key] = (f"{name}", src, acc, fasta)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "source", "accession", "fasta"])
        w.writerows(seen.values())
    have = sum(1 for v in seen.values() if v[3])
    print(f"{len(seen)} distinct proteomes; {have} with a local FASTA; {len(seen) - have} still to fetch")


if __name__ == "__main__":
    main()
