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
DONOR = set(inTable["Donor"].values)
SAMPLE = set(inTable["TestBamID"].values) # TEST
CONTROL = set(inTable["ControlBamID"].values)

# Controls split by ControlType: only "Undiluted NanoSeq" controls go through the
# read-bundle chain (control_duplex/{control}.filtered.bam). "Standard WGS" controls
# bypass it and are used directly as the matched normal. Absent column -> all NanoSeq.
if "ControlType" in inTable.columns:
    NANOSEQ_CONTROL = set(inTable.loc[inTable["ControlType"].str.strip() == "Undiluted NanoSeq", "ControlBamID"].values)
else:
    NANOSEQ_CONTROL = set(CONTROL)

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
postPath = "{donor}.runNanoSeq/tmpNanoSeq/post/"
diagnosticPath = "{donor}.runNanoSeq/tmpNanoSeq/diagnostic/"
covPath = "{donor}.runNanoSeq/tmpNanoSeq/cov/"
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


# a function that generates the expected sample BAM path for a given donor
def donor_expected_bams(donor):
    """Single source of truth for per-donor duplex/control BAM paths."""
    inTable = pd.read_csv("input.tsv", sep="\t", header=0)
    row = inTable.loc[inTable["Donor"] == donor]
    testbam = row["TestBamID"].values[0]
    control = row["ControlBamID"].values[0]
    return {
        "test": f"readBundle_duplex/{testbam}.filtered.bam",
        "control": f"control_duplex/{control}.diluted.ctrl.bam",
    }


rule all:
    input:
        # PER-BAM OUTPUTS
        # read-bundled sample BAMs
        expand("readBundle_duplex/{sample}.filtered.bam", sample=SAMPLE),
        expand("readBundle_duplex/{sample}.filtered.bam.bai", sample=SAMPLE),
        # # read-bundled control BAMs (only undiluted-NanoSeq controls; bulk WGS bypasses this)
        expand("control_duplex/{control}.filtered.bam", control=NANOSEQ_CONTROL),
        expand("control_duplex/{control}.filtered.bam.bai", control=NANOSEQ_CONTROL),
        # diluted control BAMs (matched normal consumed by variant calling, all controls)
        expand("control_duplex/{control}.diluted.ctrl.bam", control=CONTROL),
        expand("control_duplex/{control}.diluted.ctrl.bam.bai", control=CONTROL),
        # DONOR OUTPUTS
        # sample <-> control associations
        expand("DONOR/{donor}.txt",donor=DONOR), # each donor sheet will list the sample and control BAMs, as well as the intended analysis-ready BAMs
        # sample efficiency
        expand("efficiency/{donor}.tsv", donor=DONOR),
        # expand("SAMPLES/{sample}.txt", sample=SAMPLE),
        # expand("CONTROLS/{control}.txt", control=CONTROL),
        # verification
        expand("DONOR/{donor}.verified", donor=DONOR),
        # coverage diagnostic threshold
        expand(diagnosticPath + "{donor}.bulk_coverage_diagnostic.logistic_pred.txt", donor=DONOR),
        # per-chromosome coverage
        expand(covPath + "{chroms}.done", donor=DONOR, chroms=CHROMS),
        # NanoSeq post results
        expand(postPath + "results.muts.vcf.gz", donor=DONOR),
        expand(postPath + "{summary_results}.annotated.tsv", donor=DONOR, summary_results=SUMMARY_RESULT_FILES),
        # optional a4s2-isolated BAMs
        expand("a4s2_bundles/{donor}.a4s2.bam", donor=DONOR) if isolate_a4s2 else [],
        expand("a4s2_bundles/{donor}.a4s2.bam.bai", donor=DONOR) if isolate_a4s2 else [],


rule list_donor_info:
    input:
        "input.tsv",
    output:
        "DONOR/{donor}.txt",
    log:
        "logs/associate_sample_control/{donor}.log",
    resources:
        mem_mb=2000,
        runtime=5,
    threads: 1
    group:
        "associate_sample_control"
    run:
        inTable = pd.read_csv("input.tsv", sep="\t", header=0)
        testbam = inTable.loc[inTable["Donor"] == wildcards.donor, "TestBamID"].values[0]
        # donor = inTable.loc[inTable["TestBamID"] == wildcards.sample, "Donor"].values[0]
        # control = inTable.loc[inTable["TestBamID"] == wildcards.sample, "ControlBamID"].values[0]
        control = inTable.loc[inTable["Donor"] == wildcards.donor, "ControlBamID"].values[0]
        # list the expected final BAMs
        # expected_testbam = f"readBundle_duplex/{testbam}.filtered.bam"
        # expected_controlbam = f"control_duplex/{control}.filtered.bam"
        if 'Fragmentation' in inTable.columns:
            fragmentation = inTable.loc[inTable["Donor"] == wildcards.donor, "Fragmentation"].values[0]
            noise_mask = config.get("noise_"+fragmentation.lower(),NOISE)
            snp_mask = config.get("snp_"+fragmentation.lower(),SNP)
        else:
            noise_mask = NOISE
            snp_mask = SNP
        expected_bams = donor_expected_bams(wildcards.donor)
        with open(output[0], "w") as f:
            f.write("Type\tValue\n")
            f.write("Donor\t" + wildcards.donor + "\n")
            f.write("Control\t" + control + "\n")
            f.write("Noise\t" + noise_mask + "\n")
            f.write("SNP\t" + snp_mask + "\n")
            f.write("ExpectedTestBAM\t" + expected_bams["test"] + "\n")
            f.write("ExpectedControlBAM\t" + expected_bams["control"] + "\n")

# rule associate_sample_control:
#     input:
#         "input.tsv",
#     output:
#         "SAMPLES/{sample}.txt",
#     log:
#         "logs/associate_sample_control/{sample}.log",
#     resources:
#         mem_mb=2000,
#         runtime=5,
#     threads: 1
#     group:
#         "associate_sample_control"
#     run:
#         inTable = pd.read_csv("input.tsv", sep="\t", header=0)
#         donor = inTable.loc[inTable["TestBamID"] == wildcards.sample, "Donor"].values[0]
#         control = inTable.loc[inTable["TestBamID"] == wildcards.sample, "ControlBamID"].values[0]
#         if 'Fragmentation' in inTable.columns:
#             fragmentation = inTable.loc[inTable["TestBamID"] == wildcards.sample, "Fragmentation"].values[0]
#             noise_mask = config.get("noise_"+fragmentation.lower(),NOISE)
#             snp_mask = config.get("snp_"+fragmentation.lower(),SNP)
#         else:
#             noise_mask = NOISE
#             snp_mask = SNP
#         with open(output[0], "w") as f:
#             f.write("Donor\t" + donor + "\n")
#             f.write("Control\t" + control + "\n")
#             f.write("Noise\t" + noise_mask + "\n")
#             f.write("SNP\t" + snp_mask + "\n")

# rule associate_control_donor:
#     input:
#         "input.tsv",
#     output:
#         "CONTROLS/{control}.txt",
#     log:
#         "logs/associate_control_donor/{control}.log",
#     resources:
#         mem_mb=2000,
#         runtime=5,
#     threads: 1
#     group:
#         "associate_control_donor"
#     run:
#         inTable = pd.read_csv("input.tsv", sep="\t", header=0)
#         donor = inTable.loc[inTable["ControlBamID"] == wildcards.control, "Donor"].values[0]
#         control = inTable.loc[inTable["ControlBamID"] == wildcards.control, "ControlBamID"].values[0]
#         with open(output[0], "w") as f:
#             f.write("Donor\t" + donor + "\n")
#             f.write("Control\t" + control + "\n")


# Sub-snakefiles live at the project root; this entry point sits in workflows/.
include: "workflows/Snakefile.table_input.preproc_sample_bam.smk"
include: "workflows/Snakefile.table_input.preproc_control_bam.smk"
include: "workflows/Snakefile.table_input.run_nanoseq_by_donor.from_bam.smk" # modify to take in the results from fragmentation

if isolate_a4s2:
    include: "workflows/Snakefile.table_input.a4s2_bam_stats.smk"
