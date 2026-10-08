#!/bin/bash
#SBATCH -J bdcfg -p epyc -c 2 --mem 8G --time 6:00:00
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/bdcfg_%j.out
export PATH="$HOME/.pixi/bin:$PATH"
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
eval "$(pixi shell-hook -s bash --frozen --manifest-path pixi.toml)"
which datasets
python bin/build_study_config.py --study-dir studies/fungi/Bd_pangenome
