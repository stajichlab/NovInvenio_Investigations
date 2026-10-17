#!/usr/bin/env python3
"""Coordinates for Af293 proteins that the run's gene_positions table lacks.

full_v072 gene_positions has 8,686 Af293 rows, but tier1_cluster has 9,647 Af293
members; the Swiss-Prot ('sp|') entries are not placed. This script recovers a
coordinate for each unplaced member through the UniProt record's locus tag
(xrefs KEGG:afm:AFUA_xGxxxxx) and the gene feature (locus_tag=) of the study's
Af293 GFF3 (data_dir/gff3/UP000002530_330879.gff3).

Inputs: --cluster-tsv, --gene-positions, --uniprot-tsv (annotations/UP000002530.tsv.gz;
NOTE: no provenance record found for this file), --gff3. Output TSV columns:
Short, protein_id, contig, start, end, locus_tag (protein_id as in the cluster member).
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "bin"))
import starship_island_overlap as sio  # noqa: E402

SHORT = "Asfu_Af293"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cluster-tsv", required=True, type=Path)
    ap.add_argument("--gene-positions", required=True, type=Path)
    ap.add_argument("--uniprot-tsv", required=True, type=Path)
    ap.add_argument("--gff3", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    members = set()
    with sio.open_text(a.cluster_tsv) as fh:
        for line in fh:
            m = line.rstrip("\n").split("\t")[1]
            if m.startswith(SHORT + "|"):
                members.add(m[len(SHORT) + 1:])
    placed = {g[1] for g in sio.load_strain_genes(a.gene_positions, {SHORT})}
    missing = members - placed
    acc2locus = {}
    with gzip.open(a.uniprot_tsv, "rt") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            m = re.search(r"KEGG:afm:(AFUA_\w+)", p[10] if len(p) > 10 else "")
            if m:
                acc2locus[p[0]] = m.group(1)
    loc = {}
    with open(a.gff3) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 9 and p[2] == "gene":
                m = re.search(r"locus_tag=([^;]+)", p[8])
                if m:
                    loc[m.group(1)] = (p[0], int(p[3]), int(p[4]))
    n_ok = 0
    with open(a.output, "w") as out:
        out.write("Short\tprotein_id\tcontig\tstart\tend\tlocus_tag\n")
        for pid in sorted(missing):
            acc = pid.split("|")[1] if "|" in pid else pid
            lt = acc2locus.get(acc)
            if lt and lt in loc:
                c, s, e = loc[lt]
                out.write(f"{SHORT}\t{pid}\t{c}\t{s}\t{e}\t{lt}\n")
                n_ok += 1
    print(f"Af293 members {len(members)}, placed {len(placed)}, unplaced {len(missing)}, recovered {n_ok}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
