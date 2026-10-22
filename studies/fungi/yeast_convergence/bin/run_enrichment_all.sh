#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 4 --mem 16G -t 1-00:00:00
#SBATCH --job-name yc-enrich
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/yc_enrich_%j.out
#
# Novelty and loss enrichment for every finished yeast_convergence contrast,
# with the same script and settings as yeast_filamentous and zoosporic_dikarya
# (enrich_novelties_v2.py: proteome universe, mmseqs-clustered background,
# ingroup-min-frac 0.75, loss other-max-frac 0.0), so all tables are
# comparable. Output per contrast:
#   <id>/analysis/enrichment/outputs/novelty/proteome_clustered/
#   <id>/analysis/enrichment/outputs/loss/proteome_clustered/
# These are ephemeral analysis outputs (class 1, regenerable): do not commit
# them. bin/convergence_compare.py reads their enrichment_*.tsv.gz.
#
# A contrast is skipped if results/<id> does not exist yet, or if its
# outputs already exist (delete the directory to redo it).
#
# Usage: sbatch studies/fungi/yeast_convergence/bin/run_enrichment_all.sh [id ...]
set -euo pipefail

NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
FOLDER=$NII/studies/fungi/yeast_convergence
ENRICH=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY=$NII/.pixi/envs/default/bin/python3
PFAM_DAT=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
IPR_ENTRIES=/srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json

if [ "$#" -gt 0 ]; then
    IDS=("$@")
else
    mapfile -t IDS < <(tail -n +2 "$FOLDER/study_index.tsv" | cut -f1)
fi

for id in "${IDS[@]}"; do
    S=$FOLDER/$id
    R=$NII/results/$id
    O=$S/analysis/enrichment/outputs
    if [ ! -d "$R" ]; then
        echo "SKIP $id: no pipeline results at $R" >&2
        continue
    fi
    COMMON=(--results "$R" --config "$S/config.csv" --annotations "$S/annotations"
            --study "$id" --cache-dir "$O/cache" --threads "${SLURM_CPUS_ON_NODE:-4}"
            --pfam-dat "$PFAM_DAT" --interpro-entries "$IPR_ENTRIES"
            --controls "$ENRICH/actin_controls.tsv"
            --ingroup-min-frac 0.75 --universe proteome --cluster-background)
    for direction in novelty loss; do
        out=$O/$direction/proteome_clustered
        if [ -d "$out" ]; then
            echo "SKIP $id $direction: $out exists" >&2
            continue
        fi
        mkdir -p "$O/$direction"
        extra=()
        [ "$direction" = loss ] && extra=(--other-max-frac 0.0)
        echo "== $id $direction ==" >&2
        "$PY" "$ENRICH/enrich_novelties_v2.py" "${COMMON[@]}" --direction "$direction" "${extra[@]}" \
            --out-dir "$out" > "$O/$direction/run_proteome_clustered.log" 2>&1 \
            || echo "WARNING: $id $direction failed; see $O/$direction/run_proteome_clustered.log" >&2
    done
done
