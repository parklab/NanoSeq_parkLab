# THIS WORKFLOW WILL ASSUME (1) A PAIR OF FASTQS FOR THE ACTUAL NANOSEQ SAMPLE AND (2) A BAM FILE FOR THE CONTROL, TYPICALLY FROM BULK ILLUMINA WGS
import os,sys,re
import socket
import numpy as np
import pandas as pd
from itertools import product

shell.prefix("VERSION=202503.01; \
export LD_LIBRARY_PATH=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-$VERSION/lib:$LD_LIBRARY_PATH; \
export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990; \
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-$VERESION; \
export PATH=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-$VERSION/bin:$PATH; \
module load gcc/14.2.0 ; \
module load bcftools/1.21; \
module load R/4.4.2; \
module load java; \
module load perl; \
export PATH=$PATH:/n/data1/hms/dbmi/park/yujie/biobambam2/src; \
export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")


# Load config early so the FASTQ input directory can be configured by the user
configfile: "/n/data1/hms/dbmi/park/yujie/Nanoseq/NanoSeq_indelRealign_pipeline/grch38.yaml"
FASTQ_DIR = config.get("fastq_dir")

SAMPLE, = glob_wildcards(
    f"{FASTQ_DIR}/{{sample}}_1.fq.gz"
)
PAIRED_ENDS = ["1","2"]

n_jobs = config["njobs"]
dsa_threads = [i+1 for i in range(n_jobs)]

# this is for intervals with callable bases in nanoseq
INTERVALS, = glob_wildcards("intervals/{interval}.intervals.list")
n_jobs_partitioned = len(INTERVALS)
jobs_partitioned = list(range(n_jobs_partitioned+1))[1::]

# # hard coding for hg38 !! caution needed !!
# jobs_partitioned = list(range(43))[1::]

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

# WE WANT TO LIMIT TO THOSE SAMPLES THAT HAVE A MATCHED NORMAL -- WE CAN ADD THIS LATER...

#
rule all:
    input:
        # ALIGNMENTS
        expand("extracted_tags_duplex_fq/{sample}_{pe}.fq.gz", sample=SAMPLE, pe=PAIRED_ENDS),
        expand("readBundle_duplex/{sample}.filtered.bam",sample=SAMPLE),
        #
        # # JOB INITIATION
        expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",sample=SAMPLE,chroms=CHROMS),
        expand("{sample}.runNanoSeq/tmpNanoSeq/part/args.json",sample=SAMPLE),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",sample=SAMPLE,job=jobs_partitioned),
       
        # CONTROL BAM
        expand("control_bam/{sample}.ctrl.bam",sample=SAMPLE),
        expand("control_bam/{sample}.ctrl.bam.bai",sample=SAMPLE),
        # EFFICIENCY
        expand("efficiency/{sample}.tsv",sample=SAMPLE),
        #
        # # DSA
        expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.dsa.bed.gz",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # # SNVs
        expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.var",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/var", sample=SAMPLE),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel", sample=SAMPLE),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/post", sample=SAMPLE),
        # # INDELS
        expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.filtered.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # SUMMARY
        expand("{sample}.runNanoSeq/tmpNanoSeq/post/results.muts.vcf.gz",sample=SAMPLE),
        expand("{sample}.runNanoSeq/tmpNanoSeq/post/1.done",sample=SAMPLE),
        # # how do I include the analyze_nanoseq rule when I don't know the output?!!!
        # expand("{sample}.runNanoSeq/", sample=SAMPLE),
        # rules.all_analyze_nanoseq.input

# include: "Snakefile.preprocess_duplex_fastq_hg38_yujieMOD.smk" # will include workflows for extracting barcodes, mapping reads, and the tag steps (add rc and mc, mark ODs, and create read bundle tags)
include: "Snakefile.preprocess_duplex_fastq_hg38_standard.smk" # test out why indel realignment generate messy bam files after add-read-bundles and empty dsa files
## will skip a step that prepares the normal FASTQs; we have a BAM and will skip straight to analyuzing the nanoseq
# include: "Snakefile.analyze_efficiency.smk" # this is an important step to analyze the efficiency of the aligned and deduplicated nanoseq BAM: efficiency_nanoseq.pl 
## Since these nanoseq files are aligned to hg38, we need to align the ctrl bam files to hg38 as well
include: "Snakefile.ctrl_fastqs_to_38.smk"
# include: "Snakefile.analyze_nanoseq.smk" # this step will compare control BAMs vs nanoseq BAMs
include: "Snakefile.analyze_nanoseq.parallelized.smk" # this step will compare control BAMs vs nanoseq BAMs


# #(trim 3 bases, skip 4 bases, add rb & mb tags), for reads of read length 151 bps
# python extract-tags.py -a R1.fastq -b R2.fastq -c extrR1.fastq -d extrR2.fastq -m 3 -s 4 -l 151

# #(align with bwa appending rb & mb tags)
# bwa mem -C reference_genome.fa extrR1.fastq extrR2.fastq > mapped.sam
