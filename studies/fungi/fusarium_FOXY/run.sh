#!/usr/bin/bash
#SBATCH -J nf-fusarium_FOXY
#SBATCH -p exfab -A exfab -c 1 --mem 8gb --time=7-00:00:00
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/fusarium_FOXY.%j.log

# Pairwise DIAMOND novelty run for studies/fungi/fusarium_FOXY/ (PI, 2026-09-25):
# what is generally novel in the F. oxysporum species complex. 19 FOXY genomes
# (one per forma specialis) vs 4 other Fusarium lineages. Params come from this
# study's run_params.txt (read by bin/run_study.sh); site target from publish.yaml
# (fusarium_FOXY / pairwise-diamond).
#
# config.csv/data_dir were built by bin/build_study_config.py from species.csv
# (19 rows from the local 1KFG NCBI mirror, 4 fetched from NCBI; see
# DATA_MANIFEST.yaml).
#
# Pipeline: a worktree pinned at NovInvenio a72b4bc (origin/main on 2026-09-25),
# not the shared checkout, so -resume stays valid and the pipeline commit is
# recorded. Do not move or recreate the worktree while this run is in progress.
# Paths are absolute (BASH_SOURCE does not work under SLURM).

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/foxy-run-a72b4bc"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"
export NOVINVENIO_ROOT

echo "== pipeline commit: $(git -C "$NOVINVENIO_ROOT" rev-parse HEAD) =="

"$NII_ROOT/bin/run_study.sh" fungi/fusarium_FOXY \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    -c "$NII_ROOT/studies/fungi/fusarium_FOXY/foxy_resources.config" \
    -resume
