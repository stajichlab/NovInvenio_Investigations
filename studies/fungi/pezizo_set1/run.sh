#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/pezizo_set1.log

# Real run of pezizo_set1 (studies/fungi/pezizo_set1/): the first NII study,
# UniProt-reference-proteome-sourced replacement for NovInvenio's BFD-based
# pezizo5 config. Data already pulled (species.csv -> config.csv + data_dir via
# bin/build_study_config.py, see DATA_MANIFEST.yaml for full provenance).
#
# NII_PIPELINE points at the local NovInvenio checkout (not a git-fetch of
# stajichlab/NovInvenio) so relative paths inside its own nextflow.config
# (db/, configs/ symlinks) resolve against a real, already-populated checkout --
# see NovInvenio/CLAUDE.md's own note on why BASH_SOURCE-based path resolution
# breaks on SLURM; every path below is hardcoded absolute for the same reason.
#
# --run_tool diamond comes from this study's own run_params.txt (read
# automatically by bin/run_study.sh), not repeated here.
#
# Deliberately NOT passing --pfam_hmm/--swissprot_dmnd: candidates here already
# carry real UniProt DR Pfam/InterPro/GO cross-references (parsed by
# bin/extract_dat_annotations.py), so a fresh hmmscan/diamond-vs-swissprot pass
# would be redundant and could disagree with UniProt's own calls (different
# Pfam-A version, different thresholds). NII merges the UniProt-derived
# annotation into the matrix itself afterward (bin/merge_uniprot_annotations.py)
# instead of asking nf_NovInvenio's ANNOTATE_MATRIX to recompute it.
# --modelorgs_config is kept: it's a cheap lookup (no Pfam-A/hmmscan involved)
# and orthogonal to the Pfam/SwissProt redundancy concern above.

module load nextflow

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

export NII_PIPELINE="$NOVINVENIO_ROOT/main.nf"

"$NII_ROOT/bin/run_study.sh" fungi/pezizo_set1 \
    -profile slurm \
    -c "$NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config" \
    --modelorgs_config "$NOVINVENIO_ROOT/configs/modelorgs.yaml"
