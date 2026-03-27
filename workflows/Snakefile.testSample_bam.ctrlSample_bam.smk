# THIS WORKFLOW WILL ASSUME (1) A PAIR OF FASTQS FOR THE ACTUAL NANOSEQ SAMPLE AND (2) A BAM FILE FOR THE CONTROL, TYPICALLY FROM BULK ILLUMINA WGS
import os,sys,re
import socket
import numpy as np
import pandas as pd

shell.prefix("export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990; \
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.03; \
export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.03/bin; \
module load gcc/14.2.0 ; \
module load bcftools/1.21; \
export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

# SAMPLE, = glob_wildcards("input_duplex_fq/{sample}.R1.fastq.gz")
SAMPLE, = glob_wildcards("input_bam/{sample}.bam")

INTERVALS, = glob_wildcards("intervals/{interval}.intervals.list")
n_jobs_partitioned = len(INTERVALS)
jobs_partitioned = list(range(n_jobs_partitioned+1))[1::]

# chroms
CHROMS = [i+1 for i in range(24)]

# NOISE MASK AND SNPS
NOISE = config["noise_wgns"]
SNP = config["snp_wgns"]
if config["restriction_enzyme_nanoseq"] is not None:
    # don't use restriction enzyme nanoseq
    print("ANALYZING RESTRICTION-ENZYME NANOSEQ")
    NOISE = config["noise_rens"]
    SNP = config["snp_rens"]
else:
    print("ANALYZING WHOLE-GENOME NANOSEQ")

print("NOISE mask: %s"%NOISE)
print("SNP for filtering: %s"%SNP)

#
rule all:
    input:
        #
        # JOB INITIATION
        expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",sample=SAMPLE,chroms=CHROMS),
        expand("{sample}.runNanoSeq/tmpNanoSeq/part/args.json",sample=SAMPLE),
        expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",sample=SAMPLE,job=jobs_partitioned),
        #
        # expand("control_bam/{sample}.diluted.ctrl.bam",sample=SAMPLE),
        # expand("control_bam/{sample}.diluted.ctrl.bam.bai",sample=SAMPLE),
        # expand("efficiency/{sample}.tsv",sample=SAMPLE),
        # #
        # # DSA
        expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.dsa.bed.gz",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # # SNVs
        expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.var",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # # INDELS
        expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.filtered.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.done",sample=SAMPLE,job=jobs_partitioned),

rule prepare_RcMcOd_tags:
    input:
        "input_bam/{sample}.bam",
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
        runtime=60*24,
    threads: 12
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={wildcards.sample}_tmp
        mkdir -p $tmpdir
        echo -e "{wildcards.sample}" >> {log}
        bamsormadup inputformat=bam rcsupport=1 threads=1 tmpfile=$tmpdir/{wildcards.sample} < {input} > {output.outbam}
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

include: "Snakefile.a4s2_bam_stats.smk"    

# is an undiluted nanoseq used as control?
if config["control_is_undiluted_nanoseq"] is None:
    print("Nanoseq will be analyzed against a standard WGS bulk (i.e. not undiluted Nanoseq)")
    include: "Snakefile.analyze_nanoseq.against_standard_bulk.parallelized.smk"
else:
    print("Nanoseq will be analyzed against an undiluted Nanoseq library")
    include: "Snakefile.analyze_nanoseq.parallelized.smk"