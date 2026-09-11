#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/agaricomycetes_pairwise.log

# Agaricomycetes novelty-candidate study: Scom/Ccin/Agbi/Lbic (IN, Agaricales) vs
# Cneo/Umay/Rtor (OUT, Tremellomycetes/Ustilaginomycotina/Pucciniomycotina -- true
# outgroups outside Agaricomycotina; Cryptococcus/Tremellomycetes is itself inside
# Agaricomycotina, so it can only test Agaricales-level restriction, not
# Agaricomycotina-wide -- see configs/controls/Agaricales.controls.csv in the main
# NovInvenio repo). All proteomes NCBI RefSeq-sourced (built the same way as
# pezizo_set1 -- see that study's run.sh for the full NII_PIPELINE rationale).
#
# --cluster_tool pairwise (this study's default direction). Same config.csv/
# data_dir (symlinked) is reused by agaricomycetes_mmseqs and
# agaricomycetes_novelty_discovery for a direct three-way comparison, mirroring
# the pezizo_set1 / pezizo_set1_cluster pattern.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" fungi/agaricomycetes_pairwise \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --pfam_hmm "$NOVINVENIO_ROOT/db/pfam/Pfam-A.hmm" \
    --swissprot_dmnd "$NOVINVENIO_ROOT/db/uniprot/uniprot_sprot.fasta.dmnd" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml" \
    -resume
