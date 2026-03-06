##GATK: v4.1.0.0
#https://groups.google.com/forum/#!topic/snakemake/e0XNmXqL7Bg
from snakemake.utils import min_version
min_version("6.0")

shell.executable("bash")

shell.prefix(
    "source /programs/biogrids.shrc; "
    "export PATH=../NanoSeq_pipeline_install_new/bin:../conda_environment/nanoseq/bin:$PATH; "
    "export LD_LIBRARY_PATH=../conda_environment/nanoseq/lib:$LD_LIBRARY_PATH; "
)

configfile: "config/grch37_E3.yaml"

module nanoseq:
    snakefile:
        "workflows/Snakefile.testSample_bam.ctrlSample_bam_E3.smk"
    config: config

use rule * from nanoseq as run_nanoseq_*

rule all:
    input:
        rules.run_nanoseq_all.input,
    default_target: True

