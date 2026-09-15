#!/usr/bin/bash
#SBATCH -p stajichlab -c 1 --mem 64gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/cooccurrence_rescued.log

# Re-run of run_cooccurrence_full293.sh against the RESCUED presence matrix
# (bin/rescue_pass.py folded in the chunked tblastn results, 2026-09-15:
# 100,952 ABSENT -> GENOME_ONLY calls). Eligible shell+cloud family count
# nearly doubled as a result (13,181 -> 25,993), so candidate pairs scale
# ~(25993/13181)^2 =~ 3.9x (~338M vs. ~87M) -- memory and time bumped up
# from run_cooccurrence_full293.sh's 32gb/24h accordingly, not because the
# per-pair algorithm changed.

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RUN="$STUDY/results/full_293run"

export NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

cd "$NII_ROOT"
pixi run python3 "$STUDY/bin/cooccurrence.py" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --frequency_table "$RUN/frequency_table.rescued.tsv" \
    --config "$STUDY/config.csv" \
    --inventory "$RUN/strain_inventory.tsv" \
    --n_perms 200 \
    --output "$RUN/cooccurring_pairs.rescued.tsv"

echo "Done. Output: $RUN/cooccurring_pairs.rescued.tsv"
