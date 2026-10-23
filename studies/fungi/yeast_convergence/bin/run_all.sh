#!/usr/bin/bash
# Submit one SLURM job per yeast_convergence contrast. Each job runs
# bin/run_study.sh fungi/yeast_convergence/<id>, which builds config.csv +
# data_dir from that contrast's species.csv (if absent) and runs the pipeline.
# Same pipeline worktree and profile as studies/fungi/yeast_filamentous/run.sh,
# so the contrasts stay comparable with that pooled run.
#
# Before the first submit, every species.csv must be resolved:
#   pixi run python bin/ni resolve --study-dir studies/fungi/yeast_convergence/<id>
# (rows with blank Protein_Source). This script refuses to submit a contrast
# that still has blank source columns.
#
# NII_ANNOTATE_PFAM=1 is set for every job: several contrasts resolve to NCBI
# proteomes with no UniProt annotation, and the cross-contrast comparison needs
# Pfam terms for every species (CLAUDE.md "Annotating non-UniProt proteomes").
# Those hmmsearch jobs run on their own; wait for them before
# run_enrichment_all.sh.
#
# Usage:
#   studies/fungi/yeast_convergence/bin/run_all.sh                 # every contrast in study_index.tsv
#   studies/fungi/yeast_convergence/bin/run_all.sh yc_saccharomycotina yc_rev_eremothecium
#   DRY_RUN=1 studies/fungi/yeast_convergence/bin/run_all.sh       # print, do not submit
#
# Absolute paths, no BASH_SOURCE (see ~/.claude/CLAUDE.md and NII CLAUDE.md).
set -euo pipefail

NII_ROOT="${NII_ROOT:-/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations}"
NOVINVENIO_ROOT="${NOVINVENIO_ROOT:-/rhome/jstajich/projects/NovInvenio-worktrees/main-d03f5a2}"
FOLDER="$NII_ROOT/studies/fungi/yeast_convergence"
LOGDIR="$NII_ROOT/logs/slurm"
mkdir -p "$LOGDIR"

if [ "$#" -gt 0 ]; then
    IDS=("$@")
else
    mapfile -t IDS < <(tail -n +2 "$FOLDER/study_index.tsv" | cut -f1)
fi

for id in "${IDS[@]}"; do
    species_csv="$FOLDER/$id/species.csv"
    if [ ! -f "$species_csv" ]; then
        echo "ERROR: $species_csv not found (run bin/yc_batches.py render)" >&2
        exit 1
    fi
    # Protein_Source is column 6; a blank one means bin/ni resolve has not filled it.
    unresolved=$(tail -n +2 "$species_csv" | awk -F, '$6==""{print $2}' | paste -sd ';' -)
    if [ -n "$unresolved" ]; then
        echo "SKIP $id: unresolved species ($unresolved); run bin/ni resolve first" >&2
        continue
    fi
    cmd=(sbatch -p stajichlab -N 1 -n 2 --mem 8G -t 3-00:00:00
         --job-name "nf-$id"
         -o "$LOGDIR/nf_${id}_%j.out" -e "$LOGDIR/nf_${id}_%j.err"
         --wrap "source /etc/profile.d/modules.sh 2>/dev/null || true; module load nextflow; \
export NII_PIPELINE=$NOVINVENIO_ROOT/main.nf; \
export PATH=$NII_ROOT/.pixi/envs/default/bin:\$PATH; \
export NII_ANNOTATE_PFAM=1; \
$NII_ROOT/bin/run_study.sh fungi/yeast_convergence/$id \
  -profile slurm -c $NOVINVENIO_ROOT/conf/ucr_hpcc_slurm.config \
  --modelorgs_config $NOVINVENIO_ROOT/configs/modelorgs.yaml -resume")
    if [ "${DRY_RUN:-0}" = 1 ]; then
        printf '%q ' "${cmd[@]}"; echo
    else
        "${cmd[@]}"
    fi
done
