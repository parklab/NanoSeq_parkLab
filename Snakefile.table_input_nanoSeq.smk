# THIS WORKFLOW WILL ASSUME (1) A PAIR OF FASTQS FOR THE ACTUAL NANOSEQ SAMPLE AND (2) A BAM FILE FOR THE CONTROL, TYPICALLY FROM BULK ILLUMINA WGS
import os,sys,re
import socket
import numpy as np
import pandas as pd

shell.prefix("export SENTIEON_LICENSE=license01.rc.hms.harvard.edu:8990; \
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503; \
export PATH=$PATH:$SENTIEON_INSTALL_DIR/bin; \
module load gcc/14.2.0 ; \
module load bcftools/1.21; \
export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

# inputFile: "input.tsv"
if not os.path.exists("input.tsv"):
    raise ValueError("input.tsv file not found. Please provide or link a valid input.tsv file.")

inTable = pd.read_csv("input.tsv",sep="\t",header=0)
print("Input table:")
print(inTable)
# DONOR = set(inTable["Donor"].values)
SAMPLE = set(inTable["TestFastqID"].values)
CONTROL = set(inTable["ControlFastqID"].values)

# print(DONOR)
# print(SAMPLE)
# print(CONTROL)

INTERVALS, = glob_wildcards("intervals/{interval}.intervals.list")
n_jobs_partitioned = len(INTERVALS)
jobs_partitioned = list(range(n_jobs_partitioned+1))[1::]

# chroms
CHROMS = [i+1 for i in range(24)]

configfile: "config/grch38.yaml"

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


# DEFAULTS
FASTQ_DIR = "fastq/"
ALIGNED_BAM_DIR = "aligned_bam/"
CONTROL_ALIGNED_BAM_DIR = "control_aligned_bam/"
RCMCD_OD_DIR = "rcMcOd_duplex/"
MARKED_OD_DIR = RCMCD_OD_DIR
READBUNDLE_DIR = "readBundle_duplex/"
CONTROL_RCMCD_OD_DIR = "control_duplex/"
CONTROL_MARKED_OD_DIR = CONTROL_RCMCD_OD_DIR
CONTROL_READBUNDLE_DIR = CONTROL_RCMCD_OD_DIR
postPath="{sample}.runNanoSeq/tmpNanoSeq/post/"
PAIRED_ENDS = ["R1","R2"]


# # FUNCTIONS TO DETERMINE WHICH SAMPLES TO RUN, BASED ON EXISTENCE OF CONTROL FASTQS
# def check_if_sample_fastq_exists(wildcards):
#     ext=[".fastq.gz",".fq.gz",".fastq",".fq"] # possible extensions for the FASTQ files
#     inTable = pd.read_csv("input.tsv",sep="\t",header=0)
#     control = inTable.loc[inTable["TestFastqID"]==wildcards.sample,"ControlFastqID"].values[0]
#     # control exists?
#     control_fastq1 = f"{FASTQ_DIR}/{control}.R1.{ext}" # this will be used to check for the existence of the control FASTQ files, but the actual alignment rule will use the control BAM file for the downstream analysis, so we don't need to worry about the control FASTQ files for now. we can just check for their existence and then use the control BAM file for the downstream analysis. this will allow us to easily modify the alignment parameters if needed without having to modify the entire workflow.
#     control_fastq2 = f"{FASTQ_DIR}/{control}.R2.{ext}"
#     control_exists = os.path.exists(control_fastq1) and os.path.exists(control_fastq2)
#     # sample exists?
#     sample_fastq1 = f"{FASTQ_DIR}/{wildcards.sample}.R1.{ext}"
#     sample_fastq2 = f"{FASTQ_DIR}/{wildcards.sample}.R2.{ext}"
#     sample_exists = os.path.exists(sample_fastq1) and os.path.exists(sample_fastq2)
#     return control_exists and sample_exists

# def specify_final_sample_txt(wildcards):
#     files = [f"SAMPLES/{wildcards.sample}.txt"]
#     if check_if_sample_fastq_exists(wildcards):
#         return expand(files,sample=wildcards.sample)

# def specify_final_nanoseqOutput_samples(wildcards):
#     files = [f"readBundle_duplex/{wildcards.sample}.filtered.bam"]
#     if check_if_sample_fastq_exists(wildcards):
#         return expand(files,sample=wildcards.sample)
#     # else:
#     #     raise ValueError(f"FASTQ files not found for sample {wildcards.sample}. Expected FASTQ files: {FASTQ_DIR}/{wildcards.sample}.R1.{ext} and {FASTQ_DIR}/{wildcards.sample}.R2.{ext}")

# def specify_efficiencyOutput_samples(wildcards):
#     files = [f"efficiency/{wildcards.sample}.tsv"]
#     if check_if_sample_fastq_exists(wildcards):
#         return expand(files,sample=wildcards.sample)

# def specify_finalBAM_samples(wildcards):
#     # files = [f"readBundle_duplex/{wildcards.sample}.filtered.bam",f"readBundle_duplex/{wildcards.sample}.filtered.bam.bai"]
#     files = [f"readBundle_duplex/{wildcards.sample}.filtered.bam"]
#     if check_if_sample_fastq_exists(wildcards):
#         return expand(files,sample=wildcards.sample)

# # FUNCTIONS TO DETERMINE WHICH CONTROLS TO RUN, BASED ON EXISTENCE OF CONTROL FASTQS
# def check_if_control_fastq_exists(wildcards):
#     ext=[".fastq.gz",".fq.gz",".fastq",".fq"] # possible extensions for the FASTQ files
#     inTable = pd.read_csv("input.tsv",sep="\t",header=0)
#     control = inTable.loc[inTable["ControlFastqID"]==wildcards.control,"ControlFastqID"].values[0]
#     # control exists?
#     control_fastq1 = f"{FASTQ_DIR}/{control}.R1.{ext}" # this will be used to check for the existence of the control FASTQ files, but the actual alignment rule will use the control BAM file for the downstream analysis, so we don't need to worry about the control FASTQ files for now. we can just check for their existence and then use the control BAM file for the downstream analysis. this will allow us to easily modify the alignment parameters if needed without having to modify the entire workflow.
#     control_fastq2 = f"{FASTQ_DIR}/{control}.R2.{ext}"
#     control_exists = os.path.exists(control_fastq1) and os.path.exists(control_fastq2)
#     return control_exists

# def specify_final_control_txt(wildcards):
#     files = [f"CONTROLS/{wildcards.control}.txt"]
#     if check_if_control_fastq_exists(wildcards):
#         return expand(files,control=wildcards.control)

# def specify_finalBAM_controls(wildcards):
#     # files = [f"control_aligned_bam/{wildcards.control}.bam",f"control_aligned_bam/{wildcards.control}.bam.bai"]
#     files = [f"control_duplex/{wildcards.control}.diluted.ctrl.bam"]
#     if check_if_control_fastq_exists(wildcards):
#         return expand(files,control=wildcards.control)

#
rule all:
    input:
        # specify_final_sample_txt,
        # specify_final_control_txt,
        # specify_finalBAM_samples,
        # specify_finalBAM_controls,
        # specify_final_nanoseqOutput_samples,
        # specify_efficiencyOutput_samples,
        expand("SAMPLES/{sample}.txt",sample=SAMPLE),
        expand("CONTROLS/{control}.txt",control=CONTROL),
        # # alignments
        expand("aligned_bam/{sample}.bam",sample=SAMPLE),
        expand("aligned_bam/{sample}.bam.bai",sample=SAMPLE),
        expand("control_aligned_bam/{control}.bam",control=CONTROL),
        expand("control_aligned_bam/{control}.bam.bai",control=CONTROL),
        # samples' read-bundled BAMs
        expand("readBundle_duplex/{sample}.filtered.bam",sample=SAMPLE),
        expand("readBundle_duplex/{sample}.filtered.bam.bai",sample=SAMPLE),
        # controls' read-bundled BAMs
        expand("control_duplex/{control}.filtered.bam",control=CONTROL),
        expand("control_duplex/{control}.filtered.bam.bai",control=CONTROL),
        # controls' diluted BAMs
        expand("control_duplex/{control}.diluted.ctrl.bam",control=CONTROL),
        expand("control_duplex/{control}.diluted.ctrl.bam.bai",control=CONTROL),
        # nanoseq efficiency
        expand("efficiency/{sample}.tsv",sample=SAMPLE),
        # rb-isolated BAM
        expand("a4s2_bundles/{sample}.a4s2.bam",sample=SAMPLE),
        expand("a4s2_bundles/{sample}.a4s2.bam.bai",sample=SAMPLE),
        # nanoseq post
        expand(postPath+"results.muts.vcf.gz",sample=SAMPLE),
        # files to transfer
        # expand("{sample}.final_files_to_transfer.txt",sample=SAMPLE),


        #
        # # JOB INITIATION
        # expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",sample=SAMPLE,chroms=CHROMS),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/part/args.json",sample=SAMPLE),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",sample=SAMPLE,job=jobs_partitioned),
        # #
        # # expand("control_bam/{sample}.diluted.ctrl.bam",sample=SAMPLE),
        # # expand("control_bam/{sample}.diluted.ctrl.bam.bai",sample=SAMPLE),
        # # expand("efficiency/{sample}.tsv",sample=SAMPLE),
        # # #
        # # # DSA
        # expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.dsa.bed.gz",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # # # SNVs
        # expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.var",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.done",sample=SAMPLE,job=jobs_partitioned),
        # # # INDELS
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.indel.filtered.vcf.gz",sample=SAMPLE,job=jobs_partitioned),
        # expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.done",sample=SAMPLE,job=jobs_partitioned),


rule associate_sample_control:
    # modify this rule to fish out the FASTQ files for the control sample and associate it with the test sample
    # can revisit this file to get the name of the control sample and then use that to fish out the diluted control BAM file for the downstream analysis
    input:
        "input.tsv"
    output:
        "SAMPLES/{sample}.txt",
    log:
        "logs/associate_sample_control/{sample}.log",
    run:
        inTable = pd.read_csv("input.tsv",sep="\t",header=0)
        donor = inTable.loc[inTable["TestFastqID"]==wildcards.sample,"Donor"].values[0]
        control = inTable.loc[inTable["TestFastqID"]==wildcards.sample,"ControlFastqID"].values[0]
        with open(output[0],"w") as f:
            f.write("Donor\t" + donor + "\n")
            f.write("Control\t" + control + "\n")

rule associate_control_donor:
    input:
        "input.tsv"
    output:
        "CONTROLS/{control}.txt",
    log:
        "logs/associate_control_donor/{control}.log",
    run:
        inTable = pd.read_csv("input.tsv",sep="\t",header=0)
        donor = inTable.loc[inTable["ControlFastqID"]==wildcards.control,"Donor"].values[0]
        control = inTable.loc[inTable["ControlFastqID"]==wildcards.control,"ControlFastqID"].values[0]
        with open(output[0],"w") as f:
            f.write("Donor\t" + donor + "\n")
            f.write("Control\t" + control + "\n")

include: "Snakefile.table_input.align_sample.smk"
include: "Snakefile.table_input.align_control.smk"
include: "Snakefile.table_input.run_nanoseq.smk"
include: "Snakefile.table_input.a4s2_bam_stats.smk"
# include: "Snakefile.table_input.list_final_files_to_transfer.smk"