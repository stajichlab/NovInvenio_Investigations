#!/bin/bash
# Reproduce the island genome-count sweep: one SLURM job per seed, each job
# looping over all subset sizes, plus one job for the full representative set
# (N=121, deterministic, so one seed). Then run scripts/summarize.py once all
# jobs finish.
set -euo pipefail
AN=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/island_genome_count
NS="10 15 20 30 45 60 80 100"
SEEDS="1 2 3"
cd "$AN"
mkdir -p logs
for s in $SEEDS; do
    sbatch --output="logs/sweep_s${s}_%j.log" scripts/run_chain.sh "$s" $NS
done
sbatch --output="logs/full_n121_%j.log" scripts/run_chain.sh 0 121
