#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/UHM_Koxytoca.log

# Real run of UHM_Koxytoca (studies/bacteria/UHM_Koxytoca/): 16 NCBI RefSeq
# Klebsiella oxytoca outgroup references vs. 18 UHM metashot/prodigal MAG ingroup
# genomes, project UHM_Koxytoca.
#
# config.csv/data_dir/DATA_MANIFEST.yaml were NOT built by the generic
# bin/build_study_config.py (this study's proteomes are pre-existing local .faa,
# not a UniProt/GCA-fetch-per-species pull) -- they come from:
#   studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh          (outgroup genome DNA, for TBLASTN)
#   studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py --study-dir studies/bacteria/UHM_Koxytoca
# CAUTION: bin/run_study.sh's own "config.csv/data_dir missing -> rebuild" fallback
# calls the *generic* bin/build_study_config.py, which does NOT understand this
# study's species.csv schema (Source/Accession columns, not
# UniProt_Proteome_ID/GCA_Accession) -- if config.csv/data_dir/ ever get deleted,
# re-run the two commands above yourself; do not let bin/run_study.sh try to
# rebuild them.
#
# NII_PIPELINE points at the local NovInvenio checkout (not a git-fetch of
# stajichlab/NovInvenio) so relative paths inside its own nextflow.config (db/,
# conf/ symlinks) resolve against a real, already-populated checkout -- see
# NovInvenio/CLAUDE.md's own note on why BASH_SOURCE-based path resolution breaks
# on SLURM; every path below is hardcoded absolute for the same reason.
#
# --run_tool diamond, --cluster_tool pairwise, --pfam_hmm, --swissprot_dmnd, and
# --modelorgs_config all come from this study's own run_params.txt (read
# automatically by bin/run_study.sh), not repeated here.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" bacteria/UHM_Koxytoca \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    -resume
