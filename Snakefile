##GATK: v4.1.0.0
#https://groups.google.com/forum/#!topic/snakemake/e0XNmXqL7Bg
from snakemake.utils import min_version
min_version("6.0")

# packages
# make sure the conda environment `nanoseq` is loaded
shell.prefix("module load gcc/14.2.0 ; \
module load bcftools/1.21; \
export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

print("WARNING: make sure you have loaded the micromamba environment `nanoseq`")

configfile: "config/grch37.yaml"
# config: "config/grch38.yaml"

module nanoseq:
    snakefile:
        "workflows/Snakefile.testSample_fastq.ctrlSample_bam.smk"
    config: config
use rule * from nanoseq as run_nanoseq_*


rule all:
    input:
        rules.run_nanoseq_all.input,
    default_target: True
