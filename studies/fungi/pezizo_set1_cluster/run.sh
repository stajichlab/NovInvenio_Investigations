#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/pezizo_set1_cluster.log

# --cluster_tool mmseqs comparison run against the same config.csv/data_dir as
# studies/fungi/pezizo_set1/ (symlinked into this study dir), so the two runs'
# results/ directories (pezizo_set1 vs pezizo_set1_cluster) are directly
# comparable -- same species, same annotations, only the presence-matrix
# producer and the hmm_presence_cov/min_residues cutoffs differ.
#
# #85 (the "cluster member not found" crash this study first surfaced) is
# merged to main as of 2026-09-09 and validated against this exact study --
# back to the primary NovInvenio checkout.
#
# --run_tool diamond, --cluster_tool mmseqs, --hmm_presence_cov 0.3,
# --hmm_presence_min_residues 100 all come from this study's own
# run_params.txt (read automatically by bin/run_study.sh), not repeated here.
#
# Deliberately NOT passing --pfam_hmm/--swissprot_dmnd, same rationale as
# pezizo_set1/run.sh: candidates already carry real UniProt DR annotation,
# merged in afterward by bin/sync_reports.sh -> bin/merge_uniprot_annotations.py.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

# -w work_fix85: this study's original work/ directory has stale task-hash
# matches against a run pointed at a since-deleted checkout (NovInvenio_
# prod-main) -- Nextflow's task hash doesn't cover externally-called bin/
# script content, so -resume against the old work/ would replay a
# now-invalid .command.run rather than regenerate it (see the #85 fix
# validation run, 2026-09-08/09, for the full story). work_fix85 holds the
# first genuinely clean run's cache -- keep using -resume against it (not
# plain work/) for any future rerun of this study.
"$NII_ROOT/bin/run_study.sh" fungi/pezizo_set1_cluster \
    -profile slurm \
    -resume \
    -w "$NII_ROOT/.nf_launch/fungi/pezizo_set1_cluster/work_fix85" \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml"
