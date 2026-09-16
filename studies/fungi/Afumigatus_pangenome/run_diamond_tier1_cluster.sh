#!/usr/bin/bash
#SBATCH -p stajichlab -c 16 --mem 64gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/logs/diamond_tier1_cluster.log --job-name diamond_tier1

# Diamond-backend tier-1 clustering of the same full 295-strain (293 IN + 2
# OUT), 2,788,402-protein input FASTA already used for the real mmseqs-tier1
# run (results/full_293run/all_strains.fa, Short-prefixed) -- this study has
# so far ONLY ever run mmseqs2 for its real 295-strain data (mmseqs-tier1
# finished in ~14 min wall-clock, 47,983 families, 2026-09-13). No diamond
# equivalent has been run at this scale before. This job produces one so the
# benchmark scorecard (component 6 of
# notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md) can
# actually compare mmseqs vs. diamond against real ground truth, instead of
# scoring mmseqs alone.
#
# Same two-tier identity/coverage regime as the mmseqs run
# (cluster_backend.py's diamond-tier1: --approx-id 90 --member-cover 80),
# via bin/cluster_backend.py so the output TSV shape (rep\tmember) matches
# what build_presence_matrix.py / benchmark_scorecard.py already expect.
#
# Runtime is genuinely unknown ahead of time -- mmseqs's 14 min doesn't
# predict diamond's cost at the same scale, hence the generous 48h wall-clock
# ask. If this doesn't finish before the rest of the benchmark work is done,
# the benchmark scorecard is reported for mmseqs only (see
# PANGENOME_CLUSTER_PROFILE_NOTES.md's benchmark-suite section) and this job
# either finishes later for a follow-up comparison or is cancelled.

set -euo pipefail
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
mkdir -p studies/fungi/Afumigatus_pangenome/logs

pixi run python studies/fungi/Afumigatus_pangenome/bin/cluster_backend.py diamond-tier1 \
    --fasta studies/fungi/Afumigatus_pangenome/results/full_293run/all_strains.fa \
    --out_prefix studies/fungi/Afumigatus_pangenome/results/full_293run/tier1_diamond \
    --threads "${SLURM_CPUS_PER_TASK:-16}"
