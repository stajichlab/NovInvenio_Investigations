#!/usr/bin/bash
#SBATCH -p stajichlab -c 8 --mem 16gb --time=4:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/background_pfam_scan.log

# Pfam-A domain annotation for the 2,620 shell+cloud families NOT already
# covered by run_island_pfam_scan.sh's 9,224 island-member sequences --
# together these two scans cover all 11,052 families eligible for
# co-occurrence testing (cooccurrence.py's own shell+cloud selection),
# giving a real, correctly-scoped background for a domain-enrichment test
# (island members vs. the population that was actually eligible to be
# tested, NOT the whole genome -- core genes were never eligible, so a
# whole-genome background would spuriously enrich for "accessory-typical"
# domains regardless of which island is being tested).

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RUN="$STUDY/results/accessory_islands"
PFAM="/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm"

cd "$NII_ROOT"
pixi run hmmscan \
    --domtblout "$RUN/extra_background_vs_pfam.domtblout" \
    -E 1e-3 --cpu 8 \
    "$PFAM" \
    "$RUN/extra_background_reps.fa" \
    > "$RUN/extra_background_vs_pfam.hmmscan.log"

echo "Done. domtblout: $RUN/extra_background_vs_pfam.domtblout"
wc -l "$RUN/extra_background_vs_pfam.domtblout"
