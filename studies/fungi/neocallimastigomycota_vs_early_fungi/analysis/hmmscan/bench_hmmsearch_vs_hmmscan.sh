#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 4 --mem 8G -t 01:30:00 --job-name hmm-bench --constraint=milan
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/hmmbench_%j.out
# Benchmark: hmmscan vs hmmsearch on the first 500 proteins of PiroE2sp, same Pfam-A, --cut_ga.
# Compares run time and the set of (protein, Pfam accession) domain hits.
set -euo pipefail
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY="$NII_ROOT/studies/fungi/neocallimastigomycota_vs_early_fungi"
PFAM="/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
export PATH="$NII_ROOT/.pixi/envs/default/bin:$PATH"
cd "${SCRATCH:?}"
awk '/^>/{n++} n<=500' "$STUDY/data_dir/pep/PiroE2sp.pep.fa" > in.fa
echo "proteins: $(grep -c '^>' in.fa)"
hmmscan -h | sed -n 2p
echo "mpi build: $(hmmscan -h | grep -c -i mpi)"
s=$(date +%s)
hmmscan --cpu 4 --cut_ga --noali -o /dev/null --domtblout scan.domtbl "$PFAM" in.fa
echo "hmmscan seconds: $(( $(date +%s) - s ))"
s=$(date +%s)
hmmsearch --cpu 4 --cut_ga --noali -o /dev/null --domtblout search.domtbl "$PFAM" in.fa
echo "hmmsearch seconds: $(( $(date +%s) - s ))"
# domtblout: hmmscan col1 = domain name, col2 = domain acc, col4 = query (protein); hmmsearch col1 = target (protein), col4 = query name, col5 = query acc
grep -v '^#' scan.domtbl | awk '{print $4"\t"$2}' | sed 's/\.[0-9]*$//' | sort -u > scan.pairs
grep -v '^#' search.domtbl | awk '{print $1"\t"$5}' | sed 's/\.[0-9]*$//' | sort -u > search.pairs
echo "scan pairs: $(wc -l < scan.pairs)  search pairs: $(wc -l < search.pairs)  shared: $(comm -12 scan.pairs search.pairs | wc -l)"
