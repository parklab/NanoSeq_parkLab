#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 0-120:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH -o ./%j.out.txt #stdout file
#SBATCH -e ./%j.err.txt #stderr file


####SBATCH -p park

config=$1

if [[ $config == "" ]]; then
    echo -e "config was not specified. Defaulting to hg19"
    config=config/grch37.yaml
fi

echo -e "Using config file $config"

mkdir -p runlogs

snakemake --unlock -s Snakefile.from_bams.smk --configfile $config

snakemake \
--keep-going \
--rerun-incomplete \
--jobs 400 \
--cores 400 \
--rerun-triggers mtime \
-s Snakefile.from_bams.smk \
--configfile $config \
--retries 3 \
--cluster 'sbatch -p park -A park_contrib -c {threads} --mem={resources.mem_mb} -t {resources.runtime} -o ./runlogs/slurm-%A.log' # \


# --max-status-checks-per-second 0.01
