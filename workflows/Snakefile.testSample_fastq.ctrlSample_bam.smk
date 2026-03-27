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
module load perl/5.40.1; \
export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

SAMPLE, = glob_wildcards("input_duplex_fq/{sample}.R1.fastq.gz")
PAIRED_ENDS = ["R1","R2"] 
# n_jobs = config["njobs"]
# dsa_threads = [i+1 for i in range(n_jobs)]

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

# Are we using duplex or standard-bulk WGS as the control?
undiulted_nanoseq_bam_ctrl = config["undiulted_nanoseq_bam_ctrl"]
if undiulted_nanoseq_bam_ctrl is not None and undiulted_nanoseq_bam_ctrl == 1:
    print("USING UNDIULTED NANOSEQ BAM AS CONTROL")
    undiulted_nanoseq_bam_ctrl=True
else:    
    print("ASSUMING STANDARD ILLUMINA WGS AS CONTROL")
    undiulted_nanoseq_bam_ctrl=False


# WE WANT TO LIMIT TO THOSE SAMPLES THAT HAVE A MATCHED NORMAL -- WE CAN ADD THIS LATER...

#
rule all:
    input:
        # ALIGNMENTS
        expand("readBundle_duplex/{sample}.filtered.bam",sample=SAMPLE),
        expand("duplication_rates/duplication_rates_duplex.tsv"),
        # EFFICIENCY
        expand("efficiency/{sample}.tsv",sample=SAMPLE),
        expand("efficiency/{sample}.RBs",sample=SAMPLE),
        expand("efficiency/{sample}.RBs.pdf",sample=SAMPLE),
        expand("efficiency/{sample}.RBs.GC_inserts.tsv",sample=SAMPLE),
        # A4S2 BAM STATS
        expand("a4s2_bundles/{sample}.sort_by_byRb.bam",sample=SAMPLE),
        expand("a4s2_bundles/{sample}.a4s2.bam",sample=SAMPLE),
        expand("a4s2_bundles/{sample}.a4s2.bam.bai",sample=SAMPLE),
        expand("a4s2_bundles/{sample}.a4s2.covStats.tsv",sample=SAMPLE),
        expand("a4s2_bundles/{sample}.all.covStats.tsv",sample=SAMPLE),
        expand("a4s2_bundles/{sample}.a4s2.avg_depth.tsv",sample=SAMPLE),
        # JOB INITIATION
        expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",sample=SAMPLE,chroms=CHROMS),
        expand("{sample}.runNanoSeq/tmpNanoSeq/part/args.json",sample=SAMPLE),
        expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",sample=SAMPLE,job=jobs_partitioned),
        #
        # # CONTROL BAM
        expand("control_bam/{sample}.diluted.ctrl.bam",sample=SAMPLE),
        expand("control_bam/{sample}.diluted.ctrl.bam.bai",sample=SAMPLE),
        # # EFFICIENCY
        expand("efficiency/{sample}.tsv",sample=SAMPLE),
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
        # # # SSINDELS
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.w_ssIndel.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.w_ssIndel.filtered.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.w_ssIndel.done",sample=SAMPLE,job=jobs_partitioned),
        # # SUMMARY
        expand("{sample}.runNanoSeq/tmpNanoSeq/post/results.muts.vcf.gz",sample=SAMPLE),
        expand("{sample}.runNanoSeq/tmpNanoSeq/post/1.done",sample=SAMPLE),


include: "Snakefile.preprocess_duplex_fastq.smk" # will include workflows for extracting barcodes, mapping reads, and the tag steps (add rc and mc, mark ODs, and create read bundle tags)
## will skip a step that prepares the normal FASTQs; we have a BAM and will skip straight to analyuzing the nanoseq
# include: "Snakefile.analyze_efficiency.smk" # this is an important step to analyze the efficiency of the aligned and deduplicated nanoseq BAM: efficiency_nanoseq.pl 
# include: "Snakefile.analyze_nanoseq.smk" # this step will compare control BAMs vs nanoseq BAMs

include: "Snakefile.a4s2_bam_stats.smk"

if undiulted_nanoseq_bam_ctrl:
    include: "Snakefile.analyze_nanoseq_parallelized.smk" # this step will compare control BAMs vs nanoseq BAMs
else:
    include: "Snakefile.analyze_nanoseq.against_standard_bulk.parallelized.smk" 

# include: "Snakefile.analyze_nanoseq.parallelized.smk" # this step will compare control BAMs vs nanoseq BAMs


# #(trim 3 bases, skip 4 bases, add rb & mb tags), for reads of read length 151 bps
# python extract-tags.py -a R1.fastq -b R2.fastq -c extrR1.fastq -d extrR2.fastq -m 3 -s 4 -l 151

# #(align with bwa appending rb & mb tags)
# bwa mem -C reference_genome.fa extrR1.fastq extrR2.fastq > mapped.sam