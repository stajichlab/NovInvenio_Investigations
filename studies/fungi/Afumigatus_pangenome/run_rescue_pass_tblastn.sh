#!/usr/bin/bash
#SBATCH -p stajichlab -c 32 --mem 64gb --time=4-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/rescue_pass_tblastn.log

# Genome-level tblastn rescue pass (design spec component 1b) -- every
# tier-1 family representative against a single combined, Short-prefixed
# genome database covering all 295 strains, to catch coverage-failed
# protein-model absences (fragmented/split gene models, draft-assembly
# artifacts). Feeds bin/rescue_pass.py's parse_tblastn_hits/apply_rescue.
#
# One indexed tblastn run (47,983 queries vs. one ~8GB combined genome
# database), not 295 separate per-strain searches -- BLAST's own indexing
# makes "many queries against one big database" the efficient shape; 80.3%
# of the real presence matrix's 14.15M (family, strain) cells are already
# "absent" (2026-09-14 measurement), so restricting the query set to only
# currently-absent pairs would barely reduce the actual search work anyway
# (BLAST's cost scales with database size x query count regardless of which
# specific cells we already know are absent).
#
# Subject headers are Short-prefixed ("<Short>|<contig>") per
# rescue_pass.py's parse_tblastn_hits docstring, matching the same
# convention used throughout this study (build_presence_matrix.py's protein
# headers, the HAC screen's genome/proteome concatenation).
#
# -max_target_seqs 5: rescue_pass.py only needs to know WHETHER a
# qualifying hit exists per (family, strain), not every alignment -- capping
# hits per query keeps output size sane without changing the presence/
# absence call itself (min_pident/min_qcov gating happens in
# rescue_pass.py, not here).

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RESC="$STUDY/results/rescue_pass"
RUN="$STUDY/results/full_293run"

cd "$NII_ROOT"

echo "Building blast db from $RESC/all_genomes.fa ..."
pixi run makeblastdb -in "$RESC/all_genomes.fa" -dbtype nucl -out "$RESC/all_genomes_db"

echo "Running tblastn: $RUN/tier1_rep_seq.fasta vs $RESC/all_genomes_db ..."
pixi run tblastn \
    -query "$RUN/tier1_rep_seq.fasta" \
    -db "$RESC/all_genomes_db" \
    -outfmt "6 std qcovs" \
    -evalue 1e-10 \
    -max_target_seqs 5 \
    -num_threads 32 \
    -out "$RESC/tblastn_all_vs_all.tsv"

echo "Done. tblastn output: $RESC/tblastn_all_vs_all.tsv"
wc -l "$RESC/tblastn_all_vs_all.tsv"
