#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-120:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH -o ./%j.out.txt #stdout file
#SBATCH -e ./%j.err.txt #stderr file


mkdir -p runlogs

snakemake --unlock


snakemake \
--keep-going \
--rerun-incomplete \
--jobs 400 \
--rerun-triggers mtime \
--cluster 'sbatch -p park -A park_contrib -c {threads} --mem={resources.mem_mb} -t {resources.runtime} -o ./runlogs/slurm-%A.log' # \


# --max-status-checks-per-second 0.01
