#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 1 -c 32 --mem 128G -t 06:00:00 --constraint=milan
#SBATCH --job-name yc-cand-tax
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/yc_cand_tax_%j.out
#
# Taxonomic screen of the novelty and loss CANDIDATE proteins of one contrast: DIAMOND blastp
# against NCBI nr (clustered, taxonomy built in), best hit per protein. Companion of
# lichen_taxcheck.sh (whole proteomes). A non-fungal best hit among candidates points to
# contamination or horizontal transfer. Query names: <direction>:<Short>|<protein id>, so the
# summary has one row per direction and species.
# Usage: sbatch --export=ALL,CONTRAST=yc_pseudo_lichen bin/candidate_taxcheck.sh
# Summary: bin/lichen_taxcheck_summary.py <contrast>/analysis/candidate_taxcheck
# Absolute paths; no BASH_SOURCE (SLURM).
set -euo pipefail
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/diamond/2.1.12/bin:$PATH
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
C=${CONTRAST:?set CONTRAST, for example yc_pseudo_lichen}
S=$NII/studies/fungi/yeast_convergence/$C
OUT=$S/analysis/candidate_taxcheck
DB=/srv/projects/db/ncbi/diamond/current/nr_cluster_seq.dmnd
W=${SCRATCH:?}
mkdir -p "$OUT"

# protein id -> Short, from the contrast's own proteomes (config.csv: Protein, Short)
$NII/.pixi/envs/default/bin/python3 - "$S" "$NII/results/$C" "$W/q.faa" <<'PY'
import csv, sys
from pathlib import Path
s, res, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
idmap = {}
for r in csv.DictReader(open(s / "config.csv")):
    for line in open(s / "data_dir/pep" / r["Protein"]):
        if line.startswith(">"):
            idmap[line[1:].split()[0]] = r["Short"]
n = {"novelty": 0, "loss": 0}
with open(out, "w") as fo:
    for direction, fname in (("novelty", "candidates.fa"), ("loss", "loss_candidates.fa")):
        keep = False
        for line in open(res / fname):
            if line.startswith(">"):
                pid = line[1:].split()[0]
                short = idmap.get(pid)
                keep = short is not None
                if keep:
                    fo.write(f">{direction}:{short}|{pid}\n"); n[direction] += 1
                else:
                    print("WARNING: no species for", pid, file=sys.stderr)
            elif keep:
                fo.write(line)
print("queries", n, file=sys.stderr)
PY
grep '^>' "$W/q.faa" | sed 's/^>//; s/ .*//' > "$OUT/queries.txt"

diamond blastp -d "$DB" -q "$W/q.faa" -o "$W/hits.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore staxids skingdoms sphylums sscinames \
    --max-target-seqs 5 --evalue 1e-5 --fast -b 4 -c 1 --threads "${SLURM_CPUS_PER_TASK:?set by -c}" --tmpdir "$W" \
    2> "$OUT/diamond.log"
zstd -q -T0 -o "$OUT/hits.tsv.zst" "$W/hits.tsv"
diamond version > "$OUT/diamond_version.txt"
ls -l "$DB" >> "$OUT/diamond_version.txt"
