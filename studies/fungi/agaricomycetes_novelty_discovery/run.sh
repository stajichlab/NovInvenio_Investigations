#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/agaricomycetes_novelty_discovery.log

# --cluster_tool novelty_discovery comparison run, third leg of the
# agaricomycetes_pairwise / agaricomycetes_mmseqs / agaricomycetes_novelty_discovery
# three-way comparison. Own config.csv (DISCOVERY_TARGET/DISCOVERY_OUT/
# NEAR_INGROUP groups, not IN/OUT) but the same data_dir (symlinked) as the
# other two studies. --run_tool/--cluster_tool/--hmm_presence_domain_evalue come
# from this study's own run_params.txt.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" fungi/agaricomycetes_novelty_discovery \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --pfam_hmm "$NOVINVENIO_ROOT/db/pfam/Pfam-A.hmm" \
    --swissprot_dmnd "$NOVINVENIO_ROOT/db/uniprot/uniprot_sprot.fasta.dmnd" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml" \
    -resume
