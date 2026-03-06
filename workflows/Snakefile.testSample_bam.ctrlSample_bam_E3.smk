import os,sys,re
import socket
import numpy as np
import pandas as pd

SAMPLE = config.get("sample", "sample")

COV_K  = int(config.get("cov_k", 4))
N_JOBS = int(config.get("njobs", 80))

cov_jobs = list(range(1, COV_K + 1))
jobs     = list(range(1, N_JOBS + 1))

RUN_DIR = f"{SAMPLE}.runNanoSeq"
TMP     = f"{RUN_DIR}/tmpNanoSeq"

include: "Snakefile.analyze_nanoseq.parallelized_E3.smk"

rule all:
    input:
        expand(f"{TMP}/cov/{{kjob}}.done", kjob=cov_jobs),
        f"{TMP}/cov/cov.all.done",

        # Partition
        rules.run_nanoseq_part.output.args,

        # DSA
        expand(f"{TMP}/dsa/{{job}}.dsa.bed.gz", job=jobs),
        expand(f"{TMP}/dsa/{{job}}.done", job=jobs),

        # VAR
        expand(f"{TMP}/var/{{job}}.var", job=jobs),
        expand(f"{TMP}/var/{{job}}.done", job=jobs),

        # INDEL
        expand(f"{TMP}/indel/{{job}}.indel.vcf.gz", job=jobs),
        expand(f"{TMP}/indel/{{job}}.indel.filtered.vcf.gz", job=jobs),
        expand(f"{TMP}/indel/{{job}}.done", job=jobs),

        # POST
        rules.run_nanoseq_post.output.vcf,
        rules.run_nanoseq_post.output.done
