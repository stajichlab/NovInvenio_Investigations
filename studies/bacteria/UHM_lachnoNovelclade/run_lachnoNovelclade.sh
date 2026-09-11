#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/UHM_lachnoNovelclade.log

# Real run of UHM_lachnoNovelclade (studies/bacteria/UHM_lachnoNovelclade/):
# 4 unclassified UHM MAG-bin outgroup genomes vs. 5 ingroup (1 more
# unclassified MAG bin + 4 named Lachnospiraceae/Lacrimispora/Enterocloster
# reference genomes), project UHM_lachnoNovelclade.
#
# config.csv/data_dir/DATA_MANIFEST.yaml are built by the unified
# bin/build_study_config.py (species.csv uses the shared
# Protein_Source/Genome_Source/GFF3_Source schema):
#   bin/build_study_config.py --study-dir studies/bacteria/UHM_lachnoNovelclade
# bin/run_study.sh's own "config.csv/data_dir missing -> rebuild" fallback calls
# this same generic builder -- no special-case rebuild needed.
#
# --run_tool diamond, --cluster_tool pairwise, --pfam_hmm, --swissprot_dmnd all
# come from this study's own run_params.txt (read automatically by
# bin/run_study.sh), not repeated here.
#
# NII_PIPELINE points at the local NovInvenio checkout (not a git-fetch of
# stajichlab/NovInvenio) so relative paths inside its own nextflow.config (db/,
# conf/ symlinks) resolve against a real, already-populated checkout -- see
# NovInvenio/CLAUDE.md's own note on why BASH_SOURCE-based path resolution breaks
# on SLURM; every path below is hardcoded absolute for the same reason.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" bacteria/UHM_lachnoNovelclade \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    -resume
