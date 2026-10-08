#!/bin/bash
#SBATCH -J bdmash -p epyc -c 8 --mem 16G --time 4:00:00
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/bdmash_%j.out
# Mash sketch + all-vs-all distances of the 344 B. dendrobatidis genomes (k=21, s=10000).
M=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/.tools/mashenv/bin/mash
OUT=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Bd_pangenome/analysis/mash
cd ${SCRATCH:?}
ls /bigdata/stajichlab/jstajich/projects/Bd/Pangenome/dna/Batrachochytrium_dendrobatidis_*.scaffolds.fa.gz > list.txt
wc -l list.txt
$M sketch -p 8 -k 21 -s 10000 -l list.txt -o bd
$M triangle -p 8 bd.msh > $OUT/bd_mash_triangle.tsv
echo done
