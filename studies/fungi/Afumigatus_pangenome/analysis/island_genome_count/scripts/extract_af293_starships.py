#!/usr/bin/env python3
"""Extract the Af293 Starship coordinates from Table S7 of Gluck-Thaler et al.
2025 (mBio, doi:10.1128/mbio.01092-25) and map the RefSeq chromosome IDs to the
GenBank IDs used by this study's Af293 assembly (UP000002530 / CM000169-176).

The RefSeq -> GenBank map was checked by sequence length on 2026-09-29: all 8
chromosome lengths match exactly (NCBI esummary for NC_007194.1-NC_007201.1 vs
data_dir/dna/UP000002530_330879.dna.fa).

Output columns: starship_id, name, haplotype, confidence, status, chrom_refseq,
contig, start, end.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import openpyxl

# Table S7 writes RefSeq IDs with a hyphen ("NC-007194.1").
REFSEQ_TO_GENBANK = {
    "NC_007194.1": "CM000169.1",
    "NC_007195.1": "CM000170.1",
    "NC_007196.1": "CM000171.1",
    "NC_007197.1": "CM000172.1",
    "NC_007198.1": "CM000173.1",
    "NC_007199.1": "CM000174.1",
    "NC_007200.1": "CM000175.1",
    "NC_007201.1": "CM000176.1",
}

EXPECTED_HEADER = ["Confidence", "StarshipID", "Navis", "Haplotype", "Status"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    ws = openpyxl.load_workbook(args.xlsx, read_only=True)["Table S7"]
    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[1])
    if header[:5] != EXPECTED_HEADER:
        raise SystemExit(f"Table S7 header changed: {header}")
    col = {name: i for i, name in enumerate(header)}

    out_rows = []
    for r in rows[2:]:
        sid = r[col["StarshipID"]]
        if not sid or not str(sid).startswith("AF293_"):
            continue
        refseq = str(r[col["ContigID"]]).replace("-", "_")
        if refseq not in REFSEQ_TO_GENBANK:
            raise SystemExit(f"{sid}: contig {refseq} not in RefSeq->GenBank map")
        out_rows.append({
            "starship_id": sid,
            "name": r[col["Navis"]],
            "haplotype": r[col["Haplotype"]],
            "confidence": r[col["Confidence"]],
            "status": r[col["Status"]],
            "chrom_refseq": refseq,
            "contig": REFSEQ_TO_GENBANK[refseq],
            "start": int(r[col["Start"]]),
            "end": int(r[col["End"]]),
        })

    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {len(out_rows)} Af293 Starships to {args.output}")


if __name__ == "__main__":
    main()
