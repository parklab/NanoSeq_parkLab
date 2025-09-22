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
--max-status-checks-per-second 0.01 \
--jobs 400 \
--cores 400 \
--latency-wait 60 \
--rerun-triggers mtime \
--cluster 'sbatch -p short -A park -c {threads} --mem={resources.mem_mb} -t {resources.runtime} -o ./runlogs/slurm-%A.log' # \


