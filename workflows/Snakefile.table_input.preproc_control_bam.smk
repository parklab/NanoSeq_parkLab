# BAM-input variant of Snakefile.table_input.align_control.smk.
# Assumes a pre-aligned control BAM exists at input_control_bam/{control}.bam
# (control IDs are drawn from the ControlFastqID column of input.tsv).
# If the BAM is specified as undiluted NanoSeq, then we'll proceed with name-sort -> RcMcOd -> optical-duplicate marking -> read-bundle tagging before we dilute the normal
# else, we just jump straight to dilute-normal
# Skips extract_tags + bwa alignment in any case.

def get_control_input_bam(wildcards):
    control_bam = f"input_control_bam/{wildcards.control}.bam"
    if not os.path.exists(control_bam):
        raise ValueError(f"Control BAM file not found for control {wildcards.control}. Expected: {control_bam}")
    return control_bam


# # modify rules to proceed straight to dilution
# rule dilute_normal:
#     input:
#         get_control_input_bam,
#     output:
#         diluted_control="control_duplex/{control}.diluted.ctrl.bam",
#         diluted_control_bai="control_duplex/{control}.diluted.ctrl.bam.bai",
#     benchmark:
#         "benchmarks/dilute_normal/{control}.txt"
#     log:
#         "logs/dilute_normal/{control}.log"
#     params:
#         dilution_factor=0.1,
#     resources:
#         mem_mb=10000,
#         runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=240, increment=120, label="dilute_normal"),
#     threads: 1
#     conda:
#         workflow.source_path("../envs/nanoseq_snakemake.yml")
#     group:
#         "dilute_normal"
#     shell:
#         """
#         PATH=$PATH:$PWD/bin/
#         randomreadinbundle -I {input} -O {output.diluted_control} && \
#         samtools index {output.diluted_control}
#         """


rule name_sort_control_bam:
    input:
        get_control_input_bam,
    output:
        tmpdir=temp(directory("{control}_ctrl_tmp/")),
        outbam=temp("control_duplex/{control}.name_sorted.bam"),
    benchmark:
        "benchmarks/name_sort_control_bam/{control}.txt"
    log:
        "logs/name_sort_control_bam/{control}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=60*8, increment=60*4, label="name_sort_control_bam"),
    threads: 12
    group:
        "name_sort_control_bam"
    shell:
        """
        mkdir -p {output.tmpdir}
        samtools sort -n -@ {threads} -T {output.tmpdir}/nsort -o {output.outbam} {input}
        """


rule prepare_RcMcOd_tags_control:
    input:
        inbam=rules.name_sort_control_bam.output.outbam,
        tmpdir=rules.name_sort_control_bam.output.tmpdir,
    output:
        outbam=temp("control_duplex/{control}.od.bam"),
    benchmark:
        "benchmarks/prepare_RcMcOd_tags_control/{control}.txt"
    log:
        "logs/prepare_RcMcOd_tags_control/{control}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=60*8, increment=60*4, label="prepare_RcMcOd_tags_control"),
    threads: 12
    group:
        "prepare_RcMcOd_tags_control"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        mkdir -p {input.tmpdir}
        samtools view -h {input.inbam} | \
        bamsormadup inputformat=sam rcsupport=1 threads={threads} tmpfile={input.tmpdir}/bsmd \
        > {output.outbam} 2> {log}
        """


rule mark_optical_duplicates_control:
    input:
        rules.prepare_RcMcOd_tags_control.output.outbam,
    output:
        outbam=temp("control_duplex/{control}.marked_od.bam"),
        outbai=temp("control_duplex/{control}.marked_od.bam.bai"),
        metrics="control_duplex/{control}.mark_od_metrics.txt",
    benchmark:
        "benchmarks/mark_optical_duplicates_control/{control}.txt"
    log:
        "logs/mark_optical_duplicates_control/{control}.log"
    params:
        fasta=config["fasta"],
        tmpdir=temp(directory("{control}_ctrl_tmp2/")),
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=480, increment=60*4, label="mark_optical_duplicates_control"),
    threads: 12
    group:
        "mark_optical_duplicates_control"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={params.tmpdir}
        mkdir -p $tmpdir
        echo -e "{wildcards.control}" >> {log}
        bammarkduplicatesopt \
        I={input} \
        O={output.outbam} \
        M={output.metrics} \
        tmpfile=$tmpdir/{wildcards.control} \
        index=1 \
        optminpixeldif=2500 \
        inputformat=bam \
        outputformat=bam \
        inputthreads={threads} \
        outputthreads={threads}
        """


rule mark_read_bundles_control:
    input:
        rules.mark_optical_duplicates_control.output.outbam,
    output:
        outbam="control_duplex/{control}.filtered.bam",
        outbamInd="control_duplex/{control}.filtered.bam.bai",
    benchmark:
        "benchmarks/mark_read_bundles_control/{control}.txt"
    log:
        "logs/mark_read_bundles_control/{control}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=480, increment=60*4, label="mark_read_bundles_control"),
    threads: 36
    group:
        "mark_read_bundles_control"
    shell:
        """
        bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
        samtools index {output.outbam} || exit 1
        """


rule dilute_normal:
    input:
        ctrl_bam=rules.mark_read_bundles_control.output.outbam,
    output:
        diluted_control="control_duplex/{control}.diluted.ctrl.bam",
        diluted_control_bai="control_duplex/{control}.diluted.ctrl.bam.bai",
    benchmark:
        "benchmarks/dilute_normal/{control}.txt"
    log:
        "logs/dilute_normal/{control}.log"
    params:
        dilution_factor=0.1,
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=240, increment=120, label="dilute_normal"),
    threads: 1
    conda:
        workflow.source_path("../envs/nanoseq_snakemake.yml")
    group:
        "dilute_normal"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        randomreadinbundle -I {input.ctrl_bam} -O {output.diluted_control} && \
        samtools index {output.diluted_control}
        """
