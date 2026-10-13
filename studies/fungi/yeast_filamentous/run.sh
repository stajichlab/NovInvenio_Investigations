#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/yeast_filamentous.log

# Rerun (2026-10-13) at NovInvenio main d03f5a2 (worktree), with results/yeast_filamentous deleted first
# (search_cache is a storeDir: old outputs would be reused otherwise).
# Real run of yeast_filamentous (studies/fungi/yeast_filamentous/), built the same way as
# pezizo_set1 -- see that study's run.sh for the full rationale (NII_PIPELINE
# local-checkout requirement, why --pfam_hmm/--swissprot_dmnd are omitted in
# favor of the UniProt-derived annotation merge, --modelorgs_config kept as a
# cheap no-conflict lookup). --run_tool diamond comes from this study's own
# run_params.txt.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/rhome/jstajich/projects/NovInvenio-worktrees/main-d03f5a2"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" fungi/yeast_filamentous \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml" \
    -resume
