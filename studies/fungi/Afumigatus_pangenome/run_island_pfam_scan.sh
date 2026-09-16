#!/usr/bin/bash
#SBATCH -p stajichlab -c 16 --mem 16gb --time=4:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/island_pfam_scan.log

# Full Pfam-A domain annotation for the 9,224 unique family representative
# sequences that appear as members of a significant accessory island
# (results/accessory_islands/significant_islands.tsv) -- answers "what
# functional groups (Pfam domains) are represented in these islands", the
# natural follow-up to the captain-gene/PKS/NRPS-only screen already done.
#
# hmmscan (query=protein sequences, db=Pfam-A.hmm), not hmmsearch, since
# this is the "annotate my modest sequence set against the whole Pfam
# database" direction -- Pfam-A.hmm is already pressed (.h3f/.h3i/.h3m/.h3p
# present) so hmmscan's heuristic filters apply. 9,224 sequences x 30,134
# Pfam profiles; -c 16 for real wall-clock parallelism (only 4 cores
# available in an interactive session, not enough for this at a
# reasonable time budget).

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RUN="$STUDY/results/accessory_islands"
PFAM="/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm"

cd "$NII_ROOT"
pixi run hmmscan \
    --domtblout "$RUN/island_family_reps_vs_pfam.domtblout" \
    -E 1e-3 --cpu 16 \
    "$PFAM" \
    "$RUN/island_family_reps.fa" \
    > "$RUN/island_family_reps_vs_pfam.hmmscan.log"

echo "Done. domtblout: $RUN/island_family_reps_vs_pfam.domtblout"
wc -l "$RUN/island_family_reps_vs_pfam.domtblout"
