#!/usr/bin/env python3
"""Write the full-length sequences of every PF11001 --cut_ga hit protein to one
FASTA (IDs '<Short>|<protein_id>'), plus a per-protein length table.

Usage: extract_hit_proteins.py <hmmsearch_dir> <out.fa>
"""
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

hd = Path(sys.argv[1])
out = Path(sys.argv[2])
man = {}
with open(hd / "manifest.tsv", newline="") as fh:
    for i, r in enumerate(csv.DictReader(fh, delimiter="\t"), 1):
        man[str(i)] = r
want = defaultdict(set)
txt = subprocess.run(["zstd", "-dc", str(hd / "all_proteomes.PF11001.domtbl.idx.zst")],
                     check=True, capture_output=True, text=True).stdout
for line in txt.splitlines():
    p = line.split()
    want[p[0]].add(p[1])
seen = set()
with open(out, "w") as fo, open(str(out) + ".lengths.tsv", "w") as fl:
    fl.write("id\tlength\n")
    for idx, pids in want.items():
        r = man[idx]
        sid = None
        seq = []

        def flush():
            if sid is not None and sid.split("|", 1)[1] in pids:
                key = sid
                if key in seen:
                    print(f"duplicate id {key}", file=sys.stderr)
                    return
                seen.add(key)
                s = "".join(seq)
                fo.write(f">{key}\n")
                for k in range(0, len(s), 60):
                    fo.write(s[k:k + 60] + "\n")
                fl.write(f"{key}\t{len(s.rstrip('*'))}\n")
        with open(r["pep_path"]) as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    sid = f"{r['short']}|{line[1:].split(None, 1)[0]}"
                    seq = []
                else:
                    seq.append(line.strip())
            flush()
print(f"wrote {len(seen)} sequences", file=sys.stderr)
