import os

SAMPLE = config.get("sample", "sample")

DUPLEX = config["duplex_cram"]
NORMAL = config["normal_cram"]
FASTA  = config["fasta"]

RUNNANOSEQ  = config["runNanoSeq"]
NANOSEQ_BIN = config["nanoseq_bin"]
CONDA_ENV   = config.get("conda_env", None)

STRICT_BASH = bool(config.get("strict_bash", False))

COV_K  = int(config.get("cov_k", 4))
N_JOBS = int(config.get("njobs", 80))

SNP   = config["snp_mask"]
NOISE = config["noise_mask"]

RUN_DIR = f"{SAMPLE}.runNanoSeq"
TMP     = f"{RUN_DIR}/tmpNanoSeq"


def env_prefix(extra_path=""):
    lines = []
    if STRICT_BASH:
        lines.append("set -euo pipefail")

    lines += [
        "source /programs/biogrids.shrc",
    ]

    if CONDA_ENV:
        lines += [
            "source $(conda info --base)/etc/profile.d/conda.sh",
            f"conda activate {CONDA_ENV}",
        ]

    add_path = NANOSEQ_BIN
    if extra_path:
        add_path = f"{add_path}:{extra_path}"

    lines += [f'export PATH="$PATH:{add_path}"']

    if CONDA_ENV:
        lines += [f'export LD_LIBRARY_PATH="{CONDA_ENV}/lib:$LD_LIBRARY_PATH"']

    return "\n".join(lines)


# --------------------------
# COV
# --------------------------
rule run_nanoseq_cov:
    output:
        done=f"{TMP}/cov/{{kjob}}.done"
    params:
        RUN_DIR=RUN_DIR,
        TMP=TMP,
        RUNNANOSEQ=RUNNANOSEQ,
        DUPLEX=DUPLEX,
        NORMAL=NORMAL,
        FASTA=FASTA,
        k=COV_K,
        Q=str(config.get("cov", {}).get("Q", 0)),
        exclude=str(config.get("cov", {}).get("exclude", "MT,GL%,NC_%,hs37d5")),
        env=env_prefix(),
        done_rel=lambda wc, output: os.path.relpath(output.done, RUN_DIR),
    log:
        f"Logs/nanoseq_cov_{SAMPLE}" + "_{kjob}.log"
    resources:
        mem_mb=4000,
        runtime_min=6*60
    threads: 1
    shell:
        r"""
        mkdir -p {params.TMP}/cov Logs
        (
          {params.env}
          cd {params.RUN_DIR}

          python {params.RUNNANOSEQ} -k {params.k} -j {wildcards.kjob} \
            -A {params.NORMAL} \
            -B {params.DUPLEX} \
            -R {params.FASTA} \
            cov \
            -Q {params.Q} \
            --exclude "{params.exclude}"

          test -f {params.done_rel}
        ) &> {log}
        """


rule run_nanoseq_cov_all:
    input:
        expand(f"{TMP}/cov/{{kjob}}.done", kjob=list(range(1, COV_K + 1)))
    output:
        done=f"{TMP}/cov/cov.all.done"
    params:
        TMP=TMP
    shell:
        r"""
        mkdir -p {params.TMP}/cov
        touch {output.done}
        """


# --------------------------
# PARTITION
# --------------------------
rule run_nanoseq_part:
    input:
        rules.run_nanoseq_cov_all.output.done
    output:
        args=f"{TMP}/part/args.json"
    params:
        RUN_DIR=RUN_DIR,
        TMP=TMP,
        RUNNANOSEQ=RUNNANOSEQ,
        DUPLEX=DUPLEX,
        NORMAL=NORMAL,
        FASTA=FASTA,
        n=N_JOBS,
        env=env_prefix(),
        args_rel=lambda wc, output: os.path.relpath(output.args, RUN_DIR),
    log:
        f"Logs/nanoseq_part_{SAMPLE}.log"
    resources:
        mem_mb=10000,
        runtime_min=6*60
    threads: 1
    shell:
        r"""
        mkdir -p {params.TMP}/part Logs
        (
          {params.env}
          cd {params.RUN_DIR}

          python {params.RUNNANOSEQ} -t 1 \
            -A {params.NORMAL} \
            -B {params.DUPLEX} \
            -R {params.FASTA} \
            part \
            -n {params.n}

          test -f {params.args_rel}
        ) &> {log}
        """

# --------------------------
# DSA
# --------------------------
rule run_nanoseq_dsa:
    input:
        rules.run_nanoseq_part.output.args
    output:
        bed=f"{TMP}/dsa/{{job}}.dsa.bed.gz",
        done=f"{TMP}/dsa/{{job}}.done"
    params:
        RUN_DIR=RUN_DIR,
        TMP=TMP,
        RUNNANOSEQ=RUNNANOSEQ,
        DUPLEX=DUPLEX,
        NORMAL=NORMAL,
        FASTA=FASTA,
        k=N_JOBS,
        C=SNP,
        D=NOISE,
        d=str(config.get("dsa", {}).get("d", 2)),
        q=str(config.get("dsa", {}).get("q", 30)),
        env=env_prefix(),
        bed_rel=lambda wc, output: os.path.relpath(output.bed, RUN_DIR),
        done_rel=lambda wc, output: os.path.relpath(output.done, RUN_DIR),
    log:
        f"Logs/nanoseq_dsa_{SAMPLE}" + "_{job}.log"
    resources:
        mem_mb=12000,
        runtime_min=24*60
    threads: 1
    shell:
        r"""
        mkdir -p {params.TMP}/dsa Logs
        (
          {params.env}
          cd {params.RUN_DIR}

          python {params.RUNNANOSEQ} -k {params.k} -j {wildcards.job} \
            -A {params.NORMAL} \
            -B {params.DUPLEX} \
            -R {params.FASTA} \
            dsa \
            -C {params.C} \
            -D {params.D} \
            -d {params.d} \
            -q {params.q}

          test -f {params.bed_rel}
          test -f {params.done_rel}
        ) &> {log}
        """

JOB_IDS = range(1, N_JOBS + 1)

rule nanoseq_dsa_barrier:
    input:
        expand(f"{TMP}/dsa/{{job}}.done", job=JOB_IDS)
    output:
        f"{TMP}/dsa/ALL.done"
    shell:
        "touch {output}"


# --------------------------
# VAR
# --------------------------
rule run_nanoseq_var:
    input:
        barrier=f"{TMP}/dsa/ALL.done"
    output:
        var=f"{TMP}/var/{{job}}.var",
        done=f"{TMP}/var/{{job}}.done"
    params:
        RUN_DIR=RUN_DIR,
        TMP=TMP,
        RUNNANOSEQ=RUNNANOSEQ,
        DUPLEX=DUPLEX,
        NORMAL=NORMAL,
        FASTA=FASTA,
        k=N_JOBS,
        env=env_prefix(),
        var_rel=lambda wc, output: os.path.relpath(output.var, RUN_DIR),
        done_rel=lambda wc, output: os.path.relpath(output.done, RUN_DIR),
        **config.get("var", {})
    log:
        f"Logs/nanoseq_var_{SAMPLE}" + "_{job}.log"
    resources:
        mem_mb=12000,
        runtime_min=6*60
    threads: 1
    shell:
        r"""
        mkdir -p {params.TMP}/var Logs
        (
          {params.env}
          cd {params.RUN_DIR}

          python {params.RUNNANOSEQ} -k {params.k} -j {wildcards.job} \
            -A {params.NORMAL} \
            -B {params.DUPLEX} \
            -R {params.FASTA} \
            var \
            -a {params[a]} \
            -b {params[b]} \
            -c {params[c]} \
            -d {params[d]} \
            -f {params[f]} \
            -i {params[i]} \
            -m {params[m]} \
            -n {params[n]} \
            -p {params[p]} \
            -q {params[q]} \
            -r {params[r]} \
            -v {params[v]} \
            -x {params[x]} \
            -z {params[z]}

          test -f {params.var_rel}
          test -f {params.done_rel}
        ) &> {log}
        """


# --------------------------
# INDEL
# --------------------------
rule run_nanoseq_indel:
    input:
        rules.run_nanoseq_var.output.done
    output:
        vcf=f"{TMP}/indel/{{job}}.indel.vcf.gz",
        filtered=f"{TMP}/indel/{{job}}.indel.filtered.vcf.gz",
        done=f"{TMP}/indel/{{job}}.done"
    params:
        RUN_DIR=RUN_DIR,
        TMP=TMP,
        RUNNANOSEQ=RUNNANOSEQ,
        DUPLEX=DUPLEX,
        NORMAL=NORMAL,
        FASTA=FASTA,
        k=N_JOBS,
        env=env_prefix(extra_path=config.get("extra_path_for_indel", "")),
        vcf_rel=lambda wc, output: os.path.relpath(output.vcf, RUN_DIR),
        filtered_rel=lambda wc, output: os.path.relpath(output.filtered, RUN_DIR),
        done_rel=lambda wc, output: os.path.relpath(output.done, RUN_DIR),
        **config.get("indel", {})
    log:
        f"Logs/nanoseq_indel_{SAMPLE}" + "_{job}.log"
    resources:
        mem_mb=12000,
        runtime_min=6*60
    threads: 1
    shell:
        r"""
        mkdir -p {params.TMP}/indel Logs
        (
          {params.env}
          cd {params.RUN_DIR}

          python {params.RUNNANOSEQ} -k {params.k} -j {wildcards.job} \
            -A {params.NORMAL} \
            -B {params.DUPLEX} \
            -R {params.FASTA} \
            indel \
            -s {params[s]} \
            --rb {params[rb]} \
            --t3 {params[t3]} \
            --t5 {params[t5]} \
            -z {params[z]} \
            -a {params[a]} \
            -c {params[c]} \
            -v {params[v]}

          test -f {params.vcf_rel}
          test -f {params.filtered_rel}
          test -f {params.done_rel}
        ) &> {log}
        """


# --------------------------
# POST
# --------------------------
rule run_nanoseq_post:
    input:
        expand(f"{TMP}/indel/{{job}}.done", job=list(range(1, N_JOBS + 1)))
    output:
        vcf=f"{TMP}/post/results.muts.vcf.gz",
        done=f"{TMP}/post/1.done"
    params:
        RUN_DIR=RUN_DIR,
        TMP=TMP,
        RUNNANOSEQ=RUNNANOSEQ,
        DUPLEX=DUPLEX,
        NORMAL=NORMAL,
        FASTA=FASTA,
        env=env_prefix(),
        vcf_rel=lambda wc, output: os.path.relpath(output.vcf, RUN_DIR),
        done_rel=lambda wc, output: os.path.relpath(output.done, RUN_DIR),
    log:
        f"Logs/nanoseq_post_{SAMPLE}.log"
    resources:
        mem_mb=5000,
        runtime_min=6*60
    threads: 1
    shell:
        r"""
        mkdir -p {params.TMP}/post Logs
        (
          {params.env}
          cd {params.RUN_DIR}

          python {params.RUNNANOSEQ} -t 1 \
            -A {params.NORMAL} \
            -B {params.DUPLEX} \
            -R {params.FASTA} \
            post

          test -f {params.vcf_rel}
          test -f {params.done_rel}
        ) &> {log}
        """
