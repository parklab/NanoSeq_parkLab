# shell.prefix("export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990; \
# export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.03; \
# export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.03/bin; \
# module load gcc/14.2.0 ; \
# module load bcftools/1.21; \
# export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

rule extract_tags:
    input:
        expand("input_duplex_fq/{sample}.{pe}.fastq.gz",pe=PAIRED_ENDS,allow_missing=True)
        # "input_duplex_fq/{sample}.R1.fastq.gz",
        # "input_duplex_fq/{sample}.R2.fastq.gz",
    output:
        expand("extracted_tags_duplex_fq/{sample}.{pe}.fastq.gz",pe=PAIRED_ENDS,allow_missing=True)
        # "extracted_tags_duplex_fq/{sample}.R1.fastq.gz",
        # "extracted_tags_duplex_fq/{sample}.R2.fastq.gz",
    benchmark:
        "benchmarks/extract_tags/{sample}.txt"
    log:
        "logs/extract_tags/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        bases_trim=3,
        bases_skip=2,
        read_len=150,
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 1
    shell:
        """
        python ./python/extract_tags.py \
        -a {input[0]} \
        -b {input[1]} \
        -c {output[0]} \
        -d {output[1]} \
        -m {params.bases_trim} \
        -s {params.bases_skip} \
        -l {params.read_len}
        """

rule align:
    input:
        expand(rules.extract_tags.output,allow_missing=True),
    output:
        outsam="aligned_duplex/{sample}.sam"
    benchmark:
        "benchmarks/align/{sample}.txt"
    log:
        "logs/align/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=120,
    threads: 36
    shell:
        """
        sentieon \
        bwa mem \
        -t {threads} \
        -K 10000000 \
        -C \
        {config[fasta]} \
        {input} > {output.outsam} || exit 1
        """

rule bam_index_aligned:
    input:
        rules.align.output.outsam
    output:
        outbam="aligned_duplex/{sample}.bam",
        outbamind="aligned_duplex/{sample}.bam.bai",
    benchmark:
        "benchmarks/align/{sample}.txt"
    log:
        "logs/align/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=120,
    threads: 36
    shell:
        """
        samtools view -b -o {output.outbam} {input}
        samtools index {output.outbam}
        """

rule prepare_RcMcOd_tags:
    input:
        rules.align.output.outsam
    output:
        tmpdir=temp(directory("{sample}_tmp/")),
        outbam="rcMcOd_duplex/{sample}.od.bam"
    benchmark:
        "benchmarks/prepare_RcMcOd_tags/{sample}.txt"
    log:
        "logs/prepare_RcMcOd_tags/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads: 12
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={wildcards.sample}_tmp
        mkdir -p $tmpdir
        echo -e "{wildcards.sample}" >> {log}
        bamsormadup inputformat=sam rcsupport=1 threads=1 tmpfile=$tmpdir/{wildcards.sample} < {input} > {output.outbam}
        """
        # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 

rule mark_optical_duplicates:
    input:
        rules.prepare_RcMcOd_tags.output.outbam
    output:
        outbam="rcMcOd_duplex/{sample}.marked_od.bam",
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
        runtime=480,
    threads: 12
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
        # D={params.outduplicates} \
        # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 


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
        runtime=480,
    threads: 36
    shell:
        """
        bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
        samtools index {output.outbam} || exit 1
        """