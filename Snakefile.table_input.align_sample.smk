# function to get the sample FASTQ files for a given sample and paired-end read. this will be used in the alignment rule to specify the input FASTQ files for each sample. we will check for the existence of the FASTQ files in the extract_tags rule, so we can assume that the FASTQ files exist when we get to the align_samples rule. this function will allow us to easily modify the naming convention for the FASTQ files if needed without having to modify the entire workflow.
def get_sample_fastq(wildcards):
    sample_fastq = [f"fastq/{wildcards.sample}.{pe}.fq.gz" for pe in PAIRED_ENDS] # this will be used to check for the existence of the sample FASTQ files, but the actual alignment rule will use the control BAM file for the downstream analysis, so we don't need to worry about the sample FASTQ files for now. we can just check for their existence and then use the control BAM file for the downstream analysis.
    for i in sample_fastq:
        if not os.path.exists(i):
            raise ValueError(f"Sample FASTQ file not found for sample {wildcards.sample} and paired-end {PAIRED_ENDS}. Expected FASTQ: {i}")
    return sample_fastq

# ALIGNMENT AND READ BUNDLING FOR SAMPLES
# alignment rule:

rule extract_sample_tags:
    input:
        get_sample_fastq
    output:
        temp(expand("extracted_sample_tags_duplex_fq/{sample}.{pe}.fastq.gz",pe=PAIRED_ENDS,allow_missing=True))
    # when: check_if_fastq_exists # SNAKEMAKE V8+
    benchmark:
        "benchmarks/extract_sample_tags/{sample}.txt"
    log:
        "logs/extract_sample_tags/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        bases_trim=3,
        bases_skip=2,
        read_len=150,
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 1
    group:
        "extract_sample_tags"
    shell:
        """
        # TAG EXTRACTION
        python ./python/extract_tags.py \
        -a {input[0]} \
        -b {input[1]} \
        -c {output[0]} \
        -d {output[1]} \
        -m {params.bases_trim} \
        -s {params.bases_skip} \
        -l {params.read_len}
        """

rule align_samples:
    input:
        expand(rules.extract_sample_tags.output,pe=PAIRED_ENDS,allow_missing=True),
    output: 
        tempSam=temp("aligned_bam/{sample}.unsorted.sam"), # formerly .bam
        # sortedBam="aligned_bam/{sample}.bam",
        # sortedBamIndex="aligned_bam/{sample}.bam.bai"
    log:
        "logs/align_samples/{sample}.log"
    params:
        fasta=config["fasta"],
        fastq_dir=FASTQ_DIR,
    resources:
        mem_mb=30000,
        runtime=60*24,
    threads: 12
    group:
        "align_samples"
    shell:        
        """
        ulimit -c unlimited
        ulimit -n 8192
        # SENTIEON        
        sentieon \
        bwa mem \
        -t {threads} \
        -K 10000000 \
        -C \
        {config[fasta]} \
        {input} > {output.tempSam} || exit 1
        """

        # sentieon \
        # bwa mem \
        # -t {threads} \
        # -K 10000000 \
        # -C \
        # {config[fasta]} \
        # {input} | \
        # samtools view -@ {threads} -bS > {output.tempBam} || exit 1
        
        # samtools sort -@ {threads} -o {output.sortedBam} {output.tempBam}
        # samtools index {output.sortedBam}



        # control_bam=control_bam/${{control}}.diluted.ctrl.bam
        # control_bam_index=control_bam/${{control}}.diluted.ctrl.bai

rule prepare_RcMcOd_tags:
    input:
        # "aligned_bam/{sample}.bam",
        rules.align_samples.output.tempSam,
    output:
        tmpdir=temp(directory("{sample}_tmp/")),
        outbam=temp("rcMcOd_duplex/{sample}.od.bam")
    benchmark:
        "benchmarks/prepare_RcMcOd_tags/{sample}.txt"
    log:
        "logs/prepare_RcMcOd_tags/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=60*24,
    threads: 12
    group:
        "prepare_RcMcOd_tags"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={wildcards.sample}_tmp
        mkdir -p $tmpdir
        echo -e "{wildcards.sample}" >> {log}
        bamsormadup inputformat=sam rcsupport=1 threads={threads} tmpfile=$tmpdir/{wildcards.sample} < {input} > {output.outbam} 
        # threads=1 to avoid multithreading issues with bamsormadup
        """
        # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 

        # ulimit -c unlimited
        # ulimit -n 8192
        # tmpdir={wildcards.sample}_tmp
        # mkdir -p $tmpdir
        # echo -e "{wildcards.sample}" >> {log}
        # bamsormadup inputformat=bam rcsupport=1 threads=1 tmpfile=$tmpdir/{wildcards.sample} < {input} > {output.outbam}

rule mark_optical_duplicates:
    input:
        rules.prepare_RcMcOd_tags.output.outbam
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
        runtime=480,
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
    group:
        "mark_read_bundles"
    shell:
        """
        bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
        samtools index {output.outbam} || exit 1
        """