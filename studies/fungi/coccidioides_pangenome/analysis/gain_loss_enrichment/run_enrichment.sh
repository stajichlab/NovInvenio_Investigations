#!/usr/bin/bash
# Per-family Pfam table + gain/loss Pfam enrichment for one pangenome run (login node, < 5 min).
# Usage: run_enrichment.sh <study dir (absolute)> <run> <pangenome dir (absolute)>
# Needs <study>/analysis/gain_loss_enrichment/<run>/pfam_hmmsearch/chunk_*.domtbl.zst
# (pfam_reps_chunk.sbatch) and <study>/analysis/gain_loss/<run>/family_gain_loss.tsv.gz.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
PY=$NII/.pixi/envs/default/bin/python3
CHECK=$NII/studies/fungi/coccidioides_pangenome/analysis/gain_loss_enrichment/check_hypergeom.py
STUDY=$1; RUN=$2; PDIR=$3
OUT=$STUDY/analysis/gain_loss_enrichment/$RUN
GL=$STUDY/analysis/gain_loss/$RUN/family_gain_loss.tsv.gz
ANN=$OUT/family_pfam.tsv.gz
export PATH=$NII/.pixi/envs/default/bin:"$PATH"

$PY $NII/bin/hmm_domtbl_to_annotations.py --domtbl $OUT/pfam_hmmsearch/chunk_*.domtbl.zst \
    --fasta $PDIR/cluster/tier1_rep_seq.fasta --out $ANN \
    --pfam-hmm /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm

for spec in shell_soft_core:shell,soft_core all_accessory:soft_core,shell,cloud,singleton \
            cloud_singleton:cloud,singleton; do
    name=${spec%%:*}; cls=${spec#*:}
    $PY $NII/bin/pangenome_gain_loss_enrichment.py --gain-loss-table $GL --family-annotation $ANN \
        --out-dir $OUT/$name --classes $cls --label "$RUN"
    $PY $CHECK --gain-loss-table $GL --family-annotation $ANN \
        --enrichment $OUT/$name/enrichment_domains.tsv.gz --classes $cls --n 5 --seed 1 \
        > $OUT/$name/check_hypergeom.tsv
done
