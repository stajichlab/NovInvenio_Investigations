#!/usr/bin/bash
# Runs enrich_novelties_v2.py for zoosporic_dikarya: primary run plus two sensitivity runs,
# then the summary figure. Absolute paths (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
S=$NII/studies/fungi/zoosporic_dikarya
E=$S/analysis/enrichment
O=$E/outputs
PY=$NII/.pixi/envs/default/bin/python3
COMMON=(--results "$NII/results/zoosporic_dikarya" --config "$S/config.csv" --annotations "$S/annotations"
        --study zoosporic_dikarya --cache-dir "$O/cache" --threads "${SLURM_CPUS_ON_NODE:-2}"
        --pfam-dat /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
        --interpro-entries /srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json
        --controls "$E/actin_controls.tsv")
mkdir -p "$O"
# Primary: whole IN proteomes as universe, background clustered (mmseqs 0.3/0.8/cov-mode 0).
$PY "$E/enrich_novelties_v2.py" "${COMMON[@]}" --universe proteome --cluster-background \
    --out-dir "$O/proteome_clustered" > "$O/run_proteome_clustered.log" 2>&1
# Sensitivity 1: universe = IN-conserved proteins (pass the ingroup rule), clustered.
$PY "$E/enrich_novelties_v2.py" "${COMMON[@]}" --universe in_conserved --cluster-background \
    --out-dir "$O/sensitivity_in_conserved_clustered" > "$O/run_in_conserved.log" 2>&1
# Sensitivity 2: whole proteome, background NOT clustered (prototype behaviour).
$PY "$E/enrich_novelties_v2.py" "${COMMON[@]}" --universe proteome \
    --out-dir "$O/sensitivity_proteome_unclustered" > "$O/run_unclustered.log" 2>&1
$PY "$E/make_summary_figure.py" --enrich "$O/proteome_clustered" --config "$S/config.csv" \
    --title "Zoosporic fungi (4 IN) vs Dikarya (5 OUT): enrichment of novelty candidates (diamond pairwise, NovInvenio d03f5a2)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the 4 IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Labels: k = candidate units with term / all units with term." \
    --out "$E/figures/zoosporic_dikarya_summary"
