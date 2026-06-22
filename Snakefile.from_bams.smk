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

# configfile: "config/grch37.yaml"
configfile: "config/grch38.yaml"

print(config)

module nanoseq:
    snakefile:
        # "workflows/Snakefile.testSample_bam.ctrlSample_bam.smk"
        # "Snakefile.testSample_bam.ctrlSample_bam.table_input.smk"
        "Snakefile.testSample_bam.ctrlSample_bam.table_input.by_donor.smk"
    config: config
use rule * from nanoseq as run_nanoseq_*


rule all:
    input:
        rules.run_nanoseq_all.input,
    default_target: True