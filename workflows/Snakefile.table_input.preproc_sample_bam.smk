# BAM-input variant of Snakefile.table_input.align_sample.smk.
# Assumes a pre-aligned NanoSeq test BAM exists at input_bam/{sample}.bam
# (sample IDs are drawn from the TestFastqID column of input.tsv).
# Skips extract_tags + bwa alignment; picks up at name-sort -> RcMcOd ->
# optical-duplicate marking -> read-bundle tagging.

def get_sample_bam(wildcards):
    sample_bam = f"input_bam/{wildcards.sample}.bam"
    if not os.path.exists(sample_bam):
        raise ValueError(f"Sample BAM file not found for sample {wildcards.sample}. Expected: {sample_bam}")
    return sample_bam


rule name_sort_sample_bam:
    input:
        get_sample_bam,
    output:
        # tmpdir=temp(directory("{sample}_tmp/")),
        outbam=temp("rcMcOd_duplex/{sample}.name_sorted.bam"),
    benchmark:
        "benchmarks/name_sort_sample_bam/{sample}.txt"
    log:
        "logs/name_sort_sample_bam/{sample}.log"
    params:
        fasta=config["fasta"],
        # tmpdir=temp(directory("{sample}_tmp/")),
        tmpdir = lambda wildcards: f"{wildcards.sample}_tmp/",
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=60*8, increment=60*4, label="name_sort_sample_bam"),
    threads: 12
    group:
        "name_sort_sample_bam"
    shell:
        """
        mkdir -p {params.tmpdir}
        samtools sort -n -@ {threads} -T {params.tmpdir}/nsort -o {output.outbam} {input}
        """
        # mkdir -p {params.tmpdir}


rule prepare_RcMcOd_tags:
    input:
        inbam=rules.name_sort_sample_bam.output.outbam,
    output:
        outbam=temp("rcMcOd_duplex/{sample}.od.bam"),
    benchmark:
        "benchmarks/prepare_RcMcOd_tags/{sample}.txt"
    log:
        "logs/prepare_RcMcOd_tags/{sample}.log"
    params:
        fasta=config["fasta"],
        tmpdir = lambda wildcards: f"{wildcards.sample}_tmp/",
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=60*8, increment=60*4, label="prepare_RcMcOd_tags"),
    threads: 12
    group:
        "prepare_RcMcOd_tags"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        mkdir -p {params.tmpdir}
        samtools view -h {input.inbam} | \
        bamsormadup inputformat=sam rcsupport=1 threads={threads} tmpfile={params.tmpdir}/bsmd \
        > {output.outbam} 2> {log}
        """
        # tmpdir=rules.name_sort_sample_bam.params.tmpdir,
        # mkdir -p {params.tmpdir}

rule mark_optical_duplicates:
    input:
        rules.prepare_RcMcOd_tags.output.outbam,
    output:
        outbam=temp("rcMcOd_duplex/{sample}.marked_od.bam"),
        outbai=temp("rcMcOd_duplex/{sample}.marked_od.bam.bai"),
        metrics="rcMcOd_duplex/{sample}.mark_od_metrics.txt",
    benchmark:
        "benchmarks/mark_optical_duplicates/{sample}.txt"
    log:
        "logs/mark_optical_duplicates/{sample}.log"
    params:
        fasta=config["fasta"],
        outduplicates="rcMcOd_duplex/{sample}.duplicates.bam",
        tmpdir=temp(directory("{sample}_tmp2/")),
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=480, increment=60*4, label="mark_optical_duplicates"),
    threads: 12
    group:
        "mark_optical_duplicates"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={params.tmpdir}
        mkdir -p $tmpdir
        echo -e "{wildcards.sample}" >> {log}
        bammarkduplicatesopt \
        I={input} \
        O={output.outbam} \
        M={output.metrics} \
        tmpfile=$tmpdir/{wildcards.sample} \
        index=1 \
        optminpixeldif=2500 \
        inputformat=bam \
        outputformat=bam \
        inputthreads={threads} \
        outputthreads={threads}
        """


rule mark_read_bundles:
    input:
        rules.mark_optical_duplicates.output.outbam,
    output:
        outbam="readBundle_duplex/{sample}.filtered.bam",
        outbamInd="readBundle_duplex/{sample}.filtered.bam.bai",
    benchmark:
        "benchmarks/mark_read_bundles/{sample}.txt"
    log:
        "logs/mark_read_bundles/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=480, increment=60*4, label="mark_read_bundles"),
    threads: 36
    group:
        "mark_read_bundles"
    shell:
        """
        bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
        samtools index {output.outbam} || exit 1
        """
