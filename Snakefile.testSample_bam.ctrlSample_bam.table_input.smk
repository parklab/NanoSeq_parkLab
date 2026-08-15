# Table-input entry point for the BAM-on-BAM NanoSeq workflow.
#
# Inputs:
#   - input.tsv with columns: TestBamID, ControlBamID, Donor
#       (column names retained for compatibility with the BAM table_input
#        workflow; values here are BAM IDs.)
#   - input_bam/{sample}.bam        : pre-aligned NanoSeq test BAM(s)
#   - input_control_bam/{control}.bam: pre-aligned control / bulk BAM(s)
#
# Pipeline:
#   1. Sample BAM:    name_sort -> RcMcOd -> mark optical dups -> read-bundle
#   2. Control BAM:   name_sort -> RcMcOd -> mark optical dups -> read-bundle -> dilute
#   3. NanoSeq:       efficiency, coverage, partition, dsa, var, indel, post,
#                     all parallelized over intervals/{interval}.intervals.list
#
# All long-running rules use lambda-based, retry-aware `runtime` so SLURM
# wall-clock grows on each retry attempt (--retries N).

import os, sys, re
import socket
import numpy as np
import pandas as pd

shell.prefix("export SENTIEON_LICENSE=license01.rc.hms.harvard.edu:8990; \
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202503; \
export PATH=$PATH:$SENTIEON_INSTALL_DIR/bin; \
module load gcc/14.2.0 ; \
module load bcftools/1.21; \
export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

# ---- input table -----------------------------------------------------------
if not os.path.exists("input.tsv"):
    raise ValueError("input.tsv file not found. Please provide or link a valid input.tsv file.")

inTable = pd.read_csv("input.tsv", sep="\t", header=0)
SAMPLE = set(inTable["TestBamID"].values)
CONTROL = set(inTable["ControlBamID"].values)

INTERVALS, = glob_wildcards("intervals/{interval}.intervals.list")
n_jobs_partitioned = len(INTERVALS)
jobs_partitioned = list(range(n_jobs_partitioned + 1))[1::]

CHROMS = [i + 1 for i in range(24)]

# ---- NOISE / SNP mask selection -------------------------------------------
NOISE = config["noise_wgns"]
SNP = config["snp_wgns"]
if config["restriction_enzyme_nanoseq"] is not None:
    print("ANALYZING RESTRICTION-ENZYME NANOSEQ AS DEFAULT IF NOT SPECIFIED OTHERWISE UNDER 'Fragmentation' IN input.tsv")
    NOISE = config["noise_rens"]
    SNP = config["snp_rens"]
else:
    print("ANALYZING WHOLE-GENOME NANOSEQ AS DEFAULT IF NOT SPECIFIED OTHERWISE UNDER 'Fragmentation' IN input.tsv")
print("Default NOISE mask: %s" % NOISE)
print("Default SNP for filtering: %s" % SNP)

# ---- shared layout constants (mirror Snakefile.table_input_nanoSeq.smk) ----
ALIGNED_BAM_DIR = "aligned_bam/"
CONTROL_ALIGNED_BAM_DIR = "control_aligned_bam/"
RCMCD_OD_DIR = "rcMcOd_duplex/"
MARKED_OD_DIR = RCMCD_OD_DIR
READBUNDLE_DIR = "readBundle_duplex/"
CONTROL_RCMCD_OD_DIR = "control_duplex/"
CONTROL_MARKED_OD_DIR = CONTROL_RCMCD_OD_DIR
CONTROL_READBUNDLE_DIR = CONTROL_RCMCD_OD_DIR
postPath = "{sample}.runNanoSeq/tmpNanoSeq/post/"
diagnosticPath = "{sample}.runNanoSeq/tmpNanoSeq/diagnostic/"
covPath = "{sample}.runNanoSeq/tmpNanoSeq/cov/"
SUMMARY_RESULT_FILES = ["burdens", "callvsqpos", "coverage", "pyrvsmask", "readbundles"]

# current working dir
cwd = os.getcwd()
print("Current working directory: %s" % cwd)

# Optional outputs
isolate_a4s2 = config.get("isolate_a4s2", False)
if isolate_a4s2:
    print("CONFIG SPECIFIES TO ISOLATE A4S2 READS. FINAL BAM FILES WILL BE GENERATED.")
else:
    print("CONFIG SPECIFIES NOT TO ISOLATE A4S2 READS. FINAL BAM FILES WILL NOT BE GENERATED.")


rule all:
    input:
        # sample <-> control associations
        expand("SAMPLES/{sample}.txt", sample=SAMPLE),
        expand("CONTROLS/{control}.txt", control=CONTROL),
        # read-bundled sample BAMs
        expand("readBundle_duplex/{sample}.filtered.bam", sample=SAMPLE),
        expand("readBundle_duplex/{sample}.filtered.bam.bai", sample=SAMPLE),
        # # read-bundled control BAMs
        # expand("control_duplex/{control}.filtered.bam", control=CONTROL),
        # expand("control_duplex/{control}.filtered.bam.bai", control=CONTROL),
        # diluted control BAMs
        expand("control_duplex/{control}.diluted.ctrl.bam", control=CONTROL),
        expand("control_duplex/{control}.diluted.ctrl.bam.bai", control=CONTROL),
        # efficiency
        expand("efficiency/{sample}.tsv", sample=SAMPLE),
        # coverage diagnostic threshold
        expand(diagnosticPath + "{sample}.bulk_coverage_diagnostic.logistic_pred.txt", sample=SAMPLE),
        # per-chromosome coverage
        expand(covPath + "{chroms}.done", sample=SAMPLE, chroms=CHROMS),
        # NanoSeq post results
        expand(postPath + "results.muts.vcf.gz", sample=SAMPLE),
        expand(postPath + "{summary_results}.annotated.tsv", sample=SAMPLE, summary_results=SUMMARY_RESULT_FILES),
        # optional a4s2-isolated BAMs
        expand("a4s2_bundles/{sample}.a4s2.bam", sample=SAMPLE) if isolate_a4s2 else [],
        expand("a4s2_bundles/{sample}.a4s2.bam.bai", sample=SAMPLE) if isolate_a4s2 else [],
        # ---- mitochondrial NanoSeq ------------------------------------------
        # Coverage tracks and Circos plots key on TestBamID: donors differing
        # only in their matched control share a test BAM, so keying on Donor
        # would recompute identical mitochondrial tracks.
        expand("mito_tracks/{sample}.per_base.tsv.gz", sample=SAMPLE),
        expand("mito_plots/{sample}.circos.pdf", sample=SAMPLE),
        "mito_plots/all_samples.circos.pdf",
        "mito_tracks/summary.tsv",
        # Homopolymer stratifier for indel calls. Listed explicitly because,
        # unlike the low-MAPQ mask (pulled in by dsa -D) and the NUMT track
        # (pulled in by the plots), nothing else depends on it.
        "mito_mask/homopolymers.chrM.bed.gz",
        # Donor-keyed mito outputs, wrapped in a function because the mito
        # modules are included AFTER this rule: DONOR and mito_call_targets
        # do not exist at definition time, only at DAG-build time.
        lambda wildcards: (
            expand("mito_a4s2_bundles/{donor}.a4s2.bam", donor=DONOR)
            + expand("mito_a4s2_bundles/{donor}.a4s2.avg_depth.tsv", donor=DONOR)
            + expand("mito_efficiency/{donor}.tsv", donor=DONOR)
            + mito_call_targets()
        ),


rule associate_sample_control:
    input:
        "input.tsv",
    output:
        "SAMPLES/{sample}.txt",
    log:
        "logs/associate_sample_control/{sample}.log",
    resources:
        mem_mb=2000,
        runtime=5,
    threads: 1
    group:
        "associate_sample_control"
    run:
        inTable = pd.read_csv("input.tsv", sep="\t", header=0)
        donor = inTable.loc[inTable["TestBamID"] == wildcards.sample, "Donor"].values[0]
        control = inTable.loc[inTable["TestBamID"] == wildcards.sample, "ControlBamID"].values[0]
        if 'Fragmentation' in inTable.columns:
            fragmentation = inTable.loc[inTable["TestBamID"] == wildcards.sample, "Fragmentation"].values[0]
            noise_mask = config.get("noise_"+fragmentation.lower(),NOISE)
            snp_mask = config.get("snp_"+fragmentation.lower(),SNP)
        else:
            noise_mask = NOISE
            snp_mask = SNP
        with open(output[0], "w") as f:
            f.write("Donor\t" + donor + "\n")
            f.write("Control\t" + control + "\n")
            f.write("Noise\t" + noise_mask + "\n")
            f.write("SNP\t" + snp_mask + "\n")

rule associate_control_donor:
    input:
        "input.tsv",
    output:
        "CONTROLS/{control}.txt",
    log:
        "logs/associate_control_donor/{control}.log",
    resources:
        mem_mb=2000,
        runtime=5,
    threads: 1
    group:
        "associate_control_donor"
    run:
        inTable = pd.read_csv("input.tsv", sep="\t", header=0)
        donor = inTable.loc[inTable["ControlBamID"] == wildcards.control, "Donor"].values[0]
        control = inTable.loc[inTable["ControlBamID"] == wildcards.control, "ControlBamID"].values[0]
        with open(output[0], "w") as f:
            f.write("Donor\t" + donor + "\n")
            f.write("Control\t" + control + "\n")


# Sub-snakefiles live at the project root; this entry point sits in workflows/.
include: "workflows/Snakefile.table_input.preproc_sample_bam.smk"
include: "workflows/Snakefile.table_input.preproc_control_bam.smk"
include: "workflows/Snakefile.table_input.run_nanoseq.from_bam.smk" # modify to take in the results from fragmentation
# AUGUST 14, 2026: ANALYZE MITOCHONDRIAL GENOME COVERAGE AND VARIANTS, ISOLATING A4S2 DUPLEXES
include: "workflows/Snakefile.table_input.analyze_mito.from_bam.smk" # modify to take in the results from fragmentation

if isolate_a4s2:
    include: "workflows/Snakefile.table_input.a4s2_bam_stats.smk"
