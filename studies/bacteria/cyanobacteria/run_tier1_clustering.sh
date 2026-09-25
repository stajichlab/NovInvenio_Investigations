#!/usr/bin/bash
#SBATCH -p stajichlab -c 16 --mem 64gb --time=6:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/cyano_tier1_cluster.log

# Tier-1 mmseqs clustering "quick look" for studies/bacteria/cyanobacteria
# (Nostocales IN + Melainabacteria OUT, 544 genomes, 3,036,094 proteins).
# First attempt run interactively inside a 16GB/4-CPU shell job
# was OOM-killed at the alignment step (mmseqs --threads 256 oversubscribed
# to 4 real cores; RSS grew past the 16GB cgroup cap) -- same failure class
# already documented for this repo's cooccurrence.py run. Real SLURM batch
# job with proper resources instead, per that same precedent.

set -euo pipefail

STUDY="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/bacteria/cyanobacteria"
RUN="$STUDY/results/quick_look"
cd "$RUN"

export NOVINVENIO_ROOT=/bigdata/stajichlab/jstajich/projects/NovInvenio
rm -rf tmp_mmseqs tier1_cluster.tsv tier1_rep_seq.fasta tier1_all_seqs.fasta

pixi run --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/pixi.toml \
  python3 /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/bin/cluster_backend.py mmseqs-tier1 \
  --fasta all_strains.fa --out_prefix tier1

echo "Done: $(grep -c '' tier1_cluster.tsv 2>/dev/null || echo 0) cluster-membership lines"
