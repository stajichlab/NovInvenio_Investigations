#!/usr/bin/env bash
# SANS k-mer splits for the 529 Coccidioides pangenome strains.
# Submit from this directory:  sbatch run_sans.sh dna   |   sbatch run_sans.sh aa
#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem 240G
#SBATCH -t 2-00:00:00
#SBATCH --job-name sans-cocci
#SBATCH -o logs/sans_%x_%j.out
set -euo pipefail
MODE=${1:?dna or aa}
SANS=/bigdata/stajichlab/jstajich/software/sans/SANS
D=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/sans
cd "$D"
if [ "$MODE" = dna ]; then IN=genomes.kmt; EXTRA="-k 31"; else IN=proteins.kmt; EXTRA="-a -k 10"; fi
P=out/$MODE
/usr/bin/time -v "$SANS" -i $IN $EXTRA -T ${SLURM_CPUS_PER_TASK:-32} -t 10n -o $P.all.splits -v 2> logs/$MODE.time
"$SANS" -i $IN -s $P.all.splits -f strict -N $P.strict.nwk -o $P.strict.splits
"$SANS" -i $IN -s $P.all.splits -f weakly -o $P.weakly.splits -X $P.weakly.nex
gzip -f $P.all.splits
