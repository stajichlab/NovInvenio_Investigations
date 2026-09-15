#!/usr/bin/env python3
"""Split a multi-FASTA file into N roughly balanced chunks, for scattering
a large BLAST/tblastn/hmmsearch query set across parallel SLURM tasks.

Round-robin assignment (sequence i -> chunk i % n_chunks), not contiguous
blocks: family representative sequences vary meaningfully in length, and a
contiguous split can accidentally cluster several long sequences into one
chunk while another gets mostly short ones, unbalancing per-chunk runtime
even though sequence COUNTS are equal. Round-robin distributes that length
variance evenly without needing a real bin-packing pass.

Usage:
  split_fasta_chunks.py --input tier1_rep_seq.fasta --n_chunks 10 \\
      --out_prefix chunks/tblastn_query_chunk
  # writes chunks/tblastn_query_chunk_00.fa ... chunk_09.fa
"""
from __future__ import annotations

import argparse
from pathlib import Path


def split_fasta_chunks(input_path: str, n_chunks: int) -> list[list[str]]:
    """Return `n_chunks` lists of FASTA records (each a list of lines,
    header included), round-robin assigned in input order."""
    chunks: list[list[str]] = [[] for _ in range(n_chunks)]
    current_record: list[str] = []
    seq_index = 0
    with open(input_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_record:
                    chunks[seq_index % n_chunks].extend(current_record)
                    seq_index += 1
                current_record = [line]
            else:
                current_record.append(line)
        if current_record:
            chunks[seq_index % n_chunks].extend(current_record)
            seq_index += 1
    return chunks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--n_chunks", type=int, required=True)
    ap.add_argument("--out_prefix", required=True)
    args = ap.parse_args()

    chunks = split_fasta_chunks(args.input, args.n_chunks)
    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
    for i, lines in enumerate(chunks):
        out_path = f"{args.out_prefix}_{i:02d}.fa"
        with open(out_path, "w") as fh:
            fh.writelines(lines)
        n_seqs = sum(1 for line in lines if line.startswith(">"))
        print(f"{out_path}: {n_seqs} sequences")


if __name__ == "__main__":
    main()
