#!/usr/bin/env python3
"""Genomic positions for rescue-pass (GENOME_ONLY) presence calls -- closes
a real gap found 2026-09-15: `bin/build_family_positions.py` only knows
positions for GFF3-ANNOTATED proteins (`bin/build_gene_positions.py`), so a
family present ONLY via the genome-level tblastn rescue (a real hit, but no
annotated gene model in that strain) has no resolvable position at all.
`bin/pair_classification.py`'s `classify_pair()` counts a strain as
"co-carrying with a resolvable position" only via `family_positions.tsv`,
so those strains silently couldn't count toward physical-linkage
classification -- verified against the real rescued-matrix pair
classification run: `insufficient_data` jumped from 6.0% (pre-rescue) to
96.2% of FDR-significant pairs (post-rescue), exactly the size of gap this
closes.

A tblastn hit already carries a genomic position directly (subject
start/end on a specific contig) -- no protein_id/gene-model crosswalk
needed, unlike the GFF3 path. This re-parses the SAME per-strain tblastn
output files rescue_pass.py already consumed, keeping only hits for
(family, strain) pairs the final rescued matrix actually calls
GENOME_ONLY (not every qualifying hit -- a family already PRESENT at the
protein level in a strain already has a real GFF3-based position and
doesn't need one from here). When a GENOME_ONLY (family, strain) pair has
multiple qualifying hits (e.g. real multi-copy duplication), only the
single highest-bitscore hit's position is kept -- enough to make the pair
countable for `classify_pair`'s co-carrying-strain-count and
`linkage_fraction`'s adjacency check; representing every genomic copy of a
rescued family is the same "rank-window only, not full copy-aware
synteny" limitation `pair_classification.py`'s own module docstring
already documents as a known partial-implementation gap, not a new one
introduced here.

Usage:
  extract_rescue_positions.py --matrix presence_matrix.rescued.tsv \\
      --tblastn_tsv results/rescue_pass/per_strain_chunks/*.tblastn.tsv.zst \\
      --output rescue_positions.tsv
  # Short<TAB>family<TAB>contig<TAB>start -- feed into
  # build_family_positions.py's --rescue_positions
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import GENOME_ONLY, PresenceMatrix  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402


def parse_tblastn_best_hit_positions(
    lines: list[str], min_pident: float = 90.0, min_qcov: float = 80.0,
) -> dict[tuple[str, str], tuple[str, int, float]]:
    """Parse tblastn outfmt6+qcovs lines into
    {(family, strain): (contig, start, bitscore)} for the single
    highest-bitscore qualifying hit per (family, strain) pair. `start` is
    min(sstart, send) -- tblastn reports subject coordinates in hit
    orientation, so a minus-strand hit has sstart > send; normalizing to
    the smaller value matches `build_gene_positions.py`'s GFF3-derived
    convention (min(start), max(end)) so ranks computed later are ordered
    consistently regardless of strand."""
    best: dict[tuple[str, str], tuple[str, int, float]] = {}
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        try:
            family, subject = parts[0], parts[1]
            pident, qcovs = float(parts[2]), float(parts[-1])
            sstart, send, bitscore = int(parts[8]), int(parts[9]), float(parts[11])
        except (ValueError, IndexError):
            continue
        if pident < min_pident or qcovs < min_qcov:
            continue
        strain, _, contig = subject.partition("|")
        if not contig:
            continue
        key = (family, strain)
        start = min(sstart, send)
        existing = best.get(key)
        if existing is None or bitscore > existing[2]:
            best[key] = (contig, start, bitscore)
    return best


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument(
        "--tblastn_tsv", required=True, action="append",
        help="plain, .gz, or .zst tblastn outfmt6+qcovs file; repeat for "
        "several per-strain/chunked files",
    )
    ap.add_argument("--min_pident", type=float, default=90.0)
    ap.add_argument("--min_qcov", type=float, default=80.0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)

    best: dict[tuple[str, str], tuple[str, int, float]] = {}
    for tblastn_path in args.tblastn_tsv:
        with open_maybe_compressed(tblastn_path) as fh:
            hits = parse_tblastn_best_hit_positions(fh.readlines(), args.min_pident, args.min_qcov)
        for key, (contig, start, bitscore) in hits.items():
            existing = best.get(key)
            if existing is None or bitscore > existing[2]:
                best[key] = (contig, start, bitscore)

    n_written, n_not_genome_only = 0, 0
    with open(args.output, "w") as out:
        out.write("Short\tfamily\tcontig\tstart\n")
        for (family, strain), (contig, start, _bitscore) in sorted(best.items()):
            if strain not in matrix.strains or family not in matrix.families:
                continue
            if matrix.call(family, strain) != GENOME_ONLY:
                n_not_genome_only += 1
                continue
            out.write(f"{strain}\t{family}\t{contig}\t{start}\n")
            n_written += 1

    print(
        f"extract_rescue_positions: {n_written} GENOME_ONLY positions written, "
        f"{n_not_genome_only} qualifying hits skipped (matrix call was not "
        "GENOME_ONLY -- e.g. already PRESENT at the protein level)", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
