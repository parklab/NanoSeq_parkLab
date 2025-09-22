
rule verifyBamID:
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

rule efficiency_check: # requires installing File::Which --> maybe rewrite in python? 
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
