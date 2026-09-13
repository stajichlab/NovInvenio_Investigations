#!/usr/bin/bash
#SBATCH -p highclock -c 8 --mem 32gb --time=2-00:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/parsnp_ingroup293.log

# ParSNP core-genome SNP extraction across the 293 A. fumigatus IN strains
# (OUT-group Aslen_ref/Neofi_ref deliberately excluded -- ParSNP is a
# within-species core-genome aligner, and the Mash+PCoA prototype on this
# same config showed a single outgroup strain dominates a mixed-species
# distance structure; see assign_clades.py's docstring / the 2026-09-13
# session that discovered this).
#
# Reference: Asfu_Af293 (the community reference assembly, already in this
# study's data_dir -- UP000002530_330879).
#
# Prototype run on a 17-strain subset (2026-09-13, results/test_smoke_17strain/
# parsnp_run) showed the MUM-search/LCB-alignment core (`parsnp_core`) runs
# SINGLE-THREADED regardless of `-p`: parsnp only parallelizes across
# partitions (--max-concurrent-partitions), and with <50 genomes it collapses
# to one partition ("Too few genomes to run partitions of size >50. Running
# all genomes at once."). `-p 8` is requested here anyway for the later
# muscle/fasttree stages that DO thread, and because 293 genomes exceeds the
# ~50-genome partition threshold, so this run should actually get real
# partition-level parallelism where the 17-strain prototype could not.
# Confirm actual wall-clock/CPU behavior from this run's own log before
# assuming 8 threads bought a real speedup.
#
# --min-partition-size is NOT set here: leaving it at parsnp's default keeps
# one joint core-genome alignment across all 293 strains (what we actually
# want for a shared SNP matrix), rather than several independent per-subset
# alignments. If 48h is not enough, revisit --min-partition-size (a real
# quality/tradeoff, not a free speedup -- see assign_clades.py conversation
# notes) before reaching for more wall-clock.

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"

: "${SCRATCH:?SCRATCH is not set -- run this as a SLURM job, not interactively}"
WORK="$SCRATCH/parsnp_ingroup293"
GENOMES="$WORK/genomes"
mkdir -p "$GENOMES"

# Short-named symlinks for every IN strain's DNA FASTA, excluding OUT.
awk -F, 'NR>1 && $1=="IN" {print $7","$5}' "$STUDY/config.csv" | \
while IFS=, read -r short dna; do
    ln -sf "$STUDY/data_dir/dna/$dna" "$GENOMES/${short}.fa"
done

n_genomes=$(ls "$GENOMES" | wc -l)
echo "Prepared $n_genomes ingroup genomes in $GENOMES"

cd "$NII_ROOT"
pixi run parsnp \
    -r "$GENOMES/Asfu_Af293.fa" \
    -d "$GENOMES" \
    -o "$WORK/parsnp_run" \
    -c --vcf --use-fasttree -p 8

# Copy final results back to shared storage (results/ is gitignored -- never
# committed, per NII's DESIGN.md Sec 4/8 -- but kept for downstream use).
OUT="$STUDY/results/parsnp_ingroup293"
mkdir -p "$OUT"
cp -r "$WORK/parsnp_run"/. "$OUT/"
echo "Copied ParSNP results to $OUT"
