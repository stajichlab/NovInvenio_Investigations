#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/agaricomycetes_mmseqs.log

# --cluster_tool mmseqs comparison run against the same config.csv/data_dir as
# studies/fungi/agaricomycetes_pairwise/ (symlinked into this study dir), so the
# two runs' results/ directories are directly comparable -- same species, same
# annotations, only the presence-matrix producer and hmm_presence_* cutoffs
# differ. --run_tool/--cluster_tool/--hmm_presence_* all come from this study's
# own run_params.txt (read automatically by bin/run_study.sh), not repeated here.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" fungi/agaricomycetes_mmseqs \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --pfam_hmm "$NOVINVENIO_ROOT/db/pfam/Pfam-A.hmm" \
    --swissprot_dmnd "$NOVINVENIO_ROOT/db/uniprot/uniprot_sprot.fasta.dmnd" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml" \
    -resume
