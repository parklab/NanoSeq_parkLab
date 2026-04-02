# function to get the control FASTQ files for a given sample. this will be used in the extract_tags rule to specify the input FASTQ files for each sample. we will check for the existence of the FASTQ files in the extract_tags rule, so we can assume that the FASTQ files exist when we get to the align_samples rule. this function will allow us to easily modify the naming convention for the FASTQ files if needed without having to modify the entire workflow.
# def get_control_fastq(wildcards,pe):
#     ext={".fastq.gz",".fq.gz",".fastq",".fq"} # possible extensions for the FASTQ files
#     pe={"R1","R2"} # paired-end read identifiers
#     inTable = pd.read_csv("input.tsv",sep="\t",header=0)
#     control = wildcards.control
#     control_fastq = f"fastq/{control}.{pe}.{ext}" # this will be used to check for the existence of the control FASTQ files, but the actual alignment rule will use the control BAM file for the downstream analysis, so we don't need to worry about the control FASTQ files for now. we can just check for their existence and then use the control BAM file for the downstream analysis.
#     if os.path.exists(control_fastq):
#         return control_fastq
#     else:
#         raise ValueError(f"Control FASTQ file not found for control {control} and paired-end {pe}. Expected FASTQ: {control_fastq}")


def get_control_fastq(wildcards):
    control_fastq = [f"fastq/{wildcards.control}.{pe}.fq.gz" for pe in PAIRED_ENDS] # this will be used to check for the existence of the control FASTQ files, but the actual alignment rule will use the control BAM file for the downstream analysis, so we don't need to worry about the control FASTQ files for now. we can just check for their existence and then use the control BAM file for the downstream analysis.
    for i in control_fastq:
        if not os.path.exists(i):
            raise ValueError(f"Control FASTQ file not found for control {wildcards.control} and paired-end {PAIRED_ENDS}. Expected FASTQ: {i}")
    return control_fastq


# ALIGNMENT AND READ BUNDLING FOR CONTROLS
# alignment rule:
rule extract_control_tags:
    input:
        # expand(get_control_fastq,pe=PAIRED_ENDS,allow_missing=True)
        get_control_fastq
    output:
        expand("extracted_control_tags_duplex_fq/{control}.{pe}.fastq.gz",pe=PAIRED_ENDS,allow_missing=True)
    # when: check_if_fastq_exists
    benchmark:
        "benchmarks/extract_control_tags/{control}.txt"
    log:
        "logs/extract_control_tags/{control}.log"
    params: # a preset for the ultrashear or covaris libraries
        bases_trim=3,
        bases_skip=2,
        read_len=150,
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 1
    group:
        "extract_control_tags"
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

rule align_controls:
    input:
        expand(rules.extract_control_tags.output,pe=PAIRED_ENDS,allow_missing=True),
    output: 
        tempSam=temp("control_aligned_bam/{control}.unsorted.sam"), # formerly .bam
        # sortedBam="control_aligned_bam/{control}.bam",
        # sortedBamIndex="control_aligned_bam/{control}.bam.bai"
    log:
        "logs/align_controls/{control}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=60*24,
    threads: 12
    group:
        "align_controls"
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

        # control={wildcards.control}
        # control=$(awk 'NR==2{{print $2}}' {input})
        # control_fastq1=fastq/${{control}}.R1.fastq.gz
        # control_fastq2=fastq/${{control}}.R2.fastq.gz

        # bwa mem -t {threads} {params.fasta} $control_fastq1 $control_fastq2 | samtools view -@ {threads} -bS - > {output[0]}
        # samtools sort -@ {threads} -o {output[0]} {output[0]}
        # samtools index {output[0]} {output[1]}

        # ulimit -c unlimited
        # ulimit -n 8192
        # echo -e "Control: {wildcards.control}\nControl FASTQ 1: {input[0]}\nControl FASTQ 2: {input[1]}" >> {log}
        # # align the control FASTQ files to the reference genome and produce a BAM file that can be used for downstream analysis. we can use the control BAM file as a reference for the alignment if needed, or we can just align the control FASTQ files independently and then use the control BAM file for the downstream analysis. for now, let's just align the control FASTQ files independently and then use the control BAM file for the downstream analysis. this will allow us to easily modify the alignment parameters if needed without having to modify the entire workflow.
        # # sentieon bwa mem
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


rule prepare_RcMcOd_tags_control:
    input:
        # "control_aligned_bam/{control}.bam",
        rules.align_controls.output.tempSam,
    output:
        tmpdir=temp(directory("{control}_tmp/")),
        outbam=temp("control_duplex/{control}.od.bam"),
    benchmark:
        "benchmarks/control_duplex/{control}.txt"
    log:
        "logs/control_duplex/{control}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=60*24,
    threads: 12
    group:
        "prepare_RcMcOd_tags_control"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={wildcards.control}_tmp
        mkdir -p $tmpdir
        echo -e "{wildcards.control}" >> {log}
        bamsormadup inputformat=sam rcsupport=1 threads={threads} tmpfile=$tmpdir/{wildcards.control} < {input} > {output.outbam} 
        # threads=1 to avoid multithreading issues with bamsormadup
        """
        # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 

        # ulimit -c unlimited
        # ulimit -n 8192
        # tmpdir={wildcards.control}_tmp
        # mkdir -p $tmpdir
        # echo -e "{wildcards.control}" >> {log}
        # bamsormadup inputformat=bam rcsupport=1 threads=1 tmpfile=$tmpdir/{wildcards.control} < {input} > {output.outbam}


rule mark_optical_duplicates_control:
    input:
        rules.prepare_RcMcOd_tags_control.output.outbam
    output:
        outbam=temp("control_duplex/{control}.marked_od.bam"),
        metrics="control_duplex/{control}.mark_od_metrics.txt",
    benchmark:
        "benchmarks/mark_optical_duplicates_control/{control}.txt"
    log:
        "logs/mark_optical_duplicates_control/{control}.log"
    params:
        fasta=config["fasta"],
        # outduplicates="control_duplex/{control}.duplicates.bam",
        tmpdir=temp(directory("{control}_tmp2/")),
    resources:
        mem_mb=30000,
        runtime=480,
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
        # D={params.outduplicates} \
        # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 

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
        runtime=480,
    threads: 36
    group:
        "mark_read_bundles_control"
    shell:
        """
        bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
        samtools index {output.outbam} || exit 1
        """

# if config["dilute_ctrl_bam"]:
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
    params: # a preset for the ultrashear or covaris libraries
        dilution_factor=0.1, # 0.1 means 10%, 0.01 means 1%
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 1
    conda:
        "envs/nanoseq_snakemake.yml"
    group:
        "dilute_normal"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        randomreadinbundle -I {input.ctrl_bam} -O {output.diluted_control} && \
        samtools index {output.diluted_control}
        """


