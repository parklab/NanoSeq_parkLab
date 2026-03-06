#!/usr/bin/sh
#SBATCH --partition=bch-compute
#SBATCH -A rc-dst-bioinfo
#SBATCH --qos=priority
#SBATCH --job-name=runnanoseqs_snakemake
#SBATCH --time=24:00:00
#SBATCH -o Logs/runnanoseqs_snakemake_%j.log
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80

echo Running on `hostname`
echo Started: `date`
source /programs/biogrids.shrc

snakemake \
  --snakefile Snakefile \
  --configfile config/grch37_E3.yaml \
  --cores 80 \
  -p

echo Finished: `date`
