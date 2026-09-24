#!/usr/bin/env python3
"""For every contig in every strain, flag whether each end carries a telomeric
repeat array: >=4 tandem TTAGGG or CCCTAA within the terminal 500 bp.
Output: Short, contig, len, left_telo, right_telo (gzip TSV)."""
import gzip, re, sys
from pathlib import Path
import pandas as pd

STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/hrmA_2026-09-24"
RX = re.compile(r"(?:TTAGGG){4}|(?:CCCTAA){4}")
W = 500
ss = pd.read_csv(P / "samplesheet.with_clades.csv")
with gzip.open(A / "contig_telomere_ends.tsv.gz", "wt") as out:
    out.write("Short\tcontig\tlen\tleft_telo\tright_telo\n")
    for s in ss.Short:
        name, chunks = None, []
        def flush():
            if name is None:
                return
            seq = "".join(chunks).upper()
            out.write(f"{s}\t{name}\t{len(seq)}\t{int(bool(RX.search(seq[:W])))}\t{int(bool(RX.search(seq[-W:])))}\n")
        with open(STUDY / f"data_dir/dna/{s}.dna.fa") as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    name, chunks = line[1:].split()[0], []
                else:
                    chunks.append(line.strip())
            flush()
print("done", file=sys.stderr)
