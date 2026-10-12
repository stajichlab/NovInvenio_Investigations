#!/usr/bin/bash
#SBATCH -p stajichlab
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 3-00:00:00
#SBATCH --job-name nf-zoholo
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/nf_zoholo_%j.out
#SBATCH -e /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/nf_zoholo_%j.err

# Pairwise (diamond) run of zoosporic_holozoa_vs_dikarya. Builds config.csv + data_dir from species.csv
# (UniProt proteomes + NCBI genomes, with DATA_MANIFEST.yaml provenance) if absent,
# then runs the pipeline at NovInvenio main d03f5a2 (worktree). Submit with:
#   sbatch studies/eukaryotes/zoosporic_holozoa_vs_dikarya/run.sh
set -euo pipefail
mkdir -p /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm
source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/rhome/jstajich/projects/NovInvenio-worktrees/main-d03f5a2"
export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"
# NII pixi env: ncbi-datasets-cli, python deps for bin/build_study_config.py
export PATH="$NII_ROOT/.pixi/envs/default/bin:$PATH"

"$NII_ROOT/bin/run_study.sh" eukaryotes/zoosporic_holozoa_vs_dikarya \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml" \
    -resume
