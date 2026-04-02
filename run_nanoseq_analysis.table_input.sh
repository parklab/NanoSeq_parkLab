#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-120:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH -o ./%j.out.txt #stdout file
#SBATCH -e ./%j.err.txt #stderr file

CONFIG=$1

if [[ $CONFIG == "" ]]; then
    echo -e "config was not specified. Defaulting to hg19"
    CONFIG=config/grch37.yaml
fi

echo -e "Using config file $CONFIG"

SNAKEFILE=Snakefile.table_input_nanoSeq.smk

mkdir -p runlogs

# check that input.tsv exists and is not empty
if [[ ! -s input.tsv ]]; then
    echo -e "Error: input.tsv is missing or empty. Please provide or link a valid input.tsv file."
    exit 1
fi

echo -e "Unlocking Snakemake workflow $SNAKEFILE (if previously locked)..."

snakemake --unlock --configfile $CONFIG -s $SNAKEFILE

echo "dry run to check for errors and get an estimate of the runtime..."
snakemake -np --configfile $CONFIG -s $SNAKEFILE --cores 400 --rerun-incomplete --rerun-triggers mtime > runlogs/dry_run.txt

echo "submitting job..."
snakemake \
--keep-going \
--rerun-incomplete \
--jobs 400 \
--cores 400 \
--rerun-triggers mtime \
--configfile $CONFIG \
-s $SNAKEFILE \
--cluster 'sbatch -p park -A park_contrib -c {threads} --mem={resources.mem_mb} -t {resources.runtime} -o ./runlogs/slurm-%A.log' # \


# --max-status-checks-per-second 0.01
