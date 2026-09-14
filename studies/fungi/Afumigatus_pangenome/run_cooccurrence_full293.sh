#!/usr/bin/bash
#SBATCH -p stajichlab -c 1 --mem 32gb --time=1-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/cooccurrence_full293.log

# Co-occurrence (co-loss/co-gain) analysis across the real 295-strain
# presence matrix and frequency table -- design spec component 4, scaled per
# component 10 step 7's mitigations (bitset presence vectors, an analytic
# prefilter before the exact Fisher test, and BH-FDR against the true total
# pair count rather than materializing every pair; see cooccurrence.py's own
# module-level notes).
#
# Single-threaded on purpose: the screening/Fisher/permutation loop is not
# parallelized (no -p/--cpu flag exists for this script) -- cooccurrence.py
# is pure Python + scipy, not a multi-process tool. Requesting -c 1 keeps
# the SLURM allocation honest about what will actually be used.
#
# --mem 32gb, NOT the interactive session's default: a real run against this
# study's data (2026-09-14) was OOM-killed by a tight 8GB job memory cgroup
# after ~10 hours -- confirmed via `dmesg` (oom-kill:constraint=CONSTRAINT_
# MEMCG, not a host-wide OOM; the node had hundreds of GB free). 32gb is a
# generous multiple of the ~5-6GB the process was already using at kill time,
# not a tight fit -- cooccurrence.py's `survivors` list (millions of
# (family_a, family_b, p) tuples during screening) is the real memory driver
# at this family count (295 strains -> tens of thousands of eligible shell/
# cloud families -> tens of millions of candidate pairs).
#
# --time 1-00:00:00: the same real run logged 20.6M/25.8M pairs screened
# (80%) at 9h47m elapsed, still climbing -- a full run (screening +
# permutation testing on whatever clears FDR) plausibly needs most of a day.
# Progress is logged to stderr throughout (both the screening loop and,
# as of this script, the permutation-testing loop too) -- tail this script's
# --out log to check on a running job rather than guessing from silence.

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RUN="$STUDY/results/full_293run"

export NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"

cd "$NII_ROOT"
pixi run python3 "$STUDY/bin/cooccurrence.py" \
    --matrix "$RUN/presence_matrix.tsv" \
    --frequency_table "$RUN/frequency_table.tsv" \
    --config "$STUDY/config.csv" \
    --inventory "$RUN/strain_inventory.tsv" \
    --n_perms 200 \
    --output "$RUN/cooccurring_pairs.tsv"

echo "Done. Output: $RUN/cooccurring_pairs.tsv"
