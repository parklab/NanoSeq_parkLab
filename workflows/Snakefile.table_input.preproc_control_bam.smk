# BAM-input variant of Snakefile.table_input.align_control.smk.
# Assumes a pre-aligned control BAM exists at input_control_bam/{control}.bam
# (control IDs are drawn from the ControlFastqID column of input.tsv).
# If the BAM is specified as undiluted NanoSeq, then we'll proceed with name-sort -> RcMcOd -> optical-duplicate marking -> read-bundle tagging before we dilute the normal
# else, we just jump straight to dilute-normal
# Skips extract_tags + bwa alignment in any case.

def get_control_input_bam(wildcards):
    control_bam = f"input_control_bam/{wildcards.control}.bam"
    if not os.path.exists(control_bam):
        # "Bring your own preprocessed control": when the finished diluted
        # control is already present -- linked in from a previous analysis --
        # the raw input is not needed, and demanding it would block DAG
        # construction for a file that never has to be rebuilt. Returning no
        # input leaves the existing output up to date under mtime triggers.
        final = f"control_duplex/{wildcards.control}.diluted.ctrl.bam"
        if os.path.exists(final):
            return []
        raise ValueError(
            f"Control BAM not found for {wildcards.control}. Provide either "
            f"{control_bam} (to preprocess from scratch) or {final} "
            f"(a previously diluted control, e.g. symlinked in).")
    return control_bam


# ControlType drives how the matched normal is prepared (see Sanger NanoSeq docs):
#   "Undiluted NanoSeq" -> the normal is itself a NanoSeq library, so it must be
#       read-bundled (bamaddreadbundles) and deduplicated to one read-pair per
#       bundle (randomreadinbundle) to produce a "neat" normal.
#   "Standard WGS"      -> the normal is an ordinary bulk WGS BAM. It has no
#       read-bundle structure (bamaddreadbundles would drop every read), so it is
#       used directly as the NanoSeq -A track; at most it needs coordinate-sort +
#       duplicate-marking + index, which production (e.g. Sentieon) BAMs already have.
CONTROL_TYPE_NANOSEQ = "Undiluted NanoSeq"
CONTROL_TYPE_WGS = "Standard WGS"
_VALID_CONTROL_TYPES = {CONTROL_TYPE_NANOSEQ, CONTROL_TYPE_WGS}

def get_control_type(control):
    inTable = pd.read_csv("input.tsv", sep="\t", header=0)
    if "ControlType" not in inTable.columns:
        # Back-compat: no column -> treat every control as an undiluted NanoSeq
        # library (the historical default of this workflow).
        return CONTROL_TYPE_NANOSEQ
    vals = inTable.loc[inTable["ControlBamID"] == control, "ControlType"].values
    if len(vals) == 0:
        raise ValueError(f"Control {control} not found in input.tsv ControlBamID column")
    ctype = str(vals[0]).strip()
    if ctype not in _VALID_CONTROL_TYPES:
        raise ValueError(
            f"Unrecognized ControlType {ctype!r} for control {control}. "
            f"Expected one of {sorted(_VALID_CONTROL_TYPES)}."
        )
    return ctype

def get_normal_prep_input(wildcards):
    # The file that feeds dilute_normal depends on the control type, so that the
    # read-bundle chain is only built for NanoSeq-library controls. Wildcards must
    # be substituted here (an input function returns concrete paths, not templates).
    if get_control_type(wildcards.control) == CONTROL_TYPE_WGS:
        return get_control_input_bam(wildcards)                       # raw bulk BAM
    return f"control_duplex/{wildcards.control}.filtered.bam"          # read-bundled BAM


# # modify rules to proceed straight to dilution
# rule dilute_normal:
#     input:
#         get_control_input_bam,
#     output:
#         diluted_control="control_duplex/{control}.diluted.ctrl.bam",
#         diluted_control_bai="control_duplex/{control}.diluted.ctrl.bam.bai",
#     benchmark:
#         "benchmarks/dilute_normal/{control}.txt"
#     log:
#         "logs/dilute_normal/{control}.log"
#     params:
#         dilution_factor=0.1,
#     resources:
#         mem_mb=10000,
#         runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=240, increment=120, label="dilute_normal"),
#     threads: 1
#     conda:
#         workflow.source_path("../envs/nanoseq_snakemake.yml")
#     group:
#         "dilute_normal"
#     shell:
#         """
#         PATH=$PATH:$PWD/bin/
#         randomreadinbundle -I {input} -O {output.diluted_control} && \
#         samtools index {output.diluted_control}
#         """


rule name_sort_control_bam:
    input:
        get_control_input_bam,
    output:
        outbam=temp("control_duplex/{control}.name_sorted.bam"),
        # tmpdir=temp(directory("{control}_ctrl_tmp/")),
    benchmark:
        "benchmarks/name_sort_control_bam/{control}.txt"
    log:
        "logs/name_sort_control_bam/{control}.log"
    params:
        fasta=config["fasta"],
        tmpdir = lambda wildcards: f"{wildcards.control}_ctrl_tmp/",
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=60*8, increment=60*4, label="name_sort_control_bam"),
    threads: 12
    group:
        "name_sort_control_bam"
    shell:
        """
        mkdir -p {params.tmpdir}
        samtools sort -n -@ {threads} -T {params.tmpdir}/nsort -o {output.outbam} {input}
        """

rule prepare_RcMcOd_tags_control:
    input:
        inbam=rules.name_sort_control_bam.output.outbam,
        # tmpdir=rules.name_sort_control_bam.output.tmpdir,
    output:
        outbam=temp("control_duplex/{control}.od.bam"),
    benchmark:
        "benchmarks/prepare_RcMcOd_tags_control/{control}.txt"
    log:
        "logs/prepare_RcMcOd_tags_control/{control}.log"
    params:
        fasta=config["fasta"],
        tmpdir = lambda wildcards: f"{wildcards.control}_ctrl_tmp/",
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=60*8, increment=60*4, label="prepare_RcMcOd_tags_control"),
    threads: 12
    group:
        "prepare_RcMcOd_tags_control"
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        mkdir -p {params.tmpdir}
        samtools view -h {input.inbam} | \
        bamsormadup inputformat=sam rcsupport=1 threads={threads} tmpfile={params.tmpdir}/bsmd \
        > {output.outbam} 2> {log}
        """
        # mkdir -p {input.tmpdir}


rule mark_optical_duplicates_control:
    input:
        rules.prepare_RcMcOd_tags_control.output.outbam,
    output:
        outbam=temp("control_duplex/{control}.marked_od.bam"),
        outbai=temp("control_duplex/{control}.marked_od.bam.bai"),
        metrics="control_duplex/{control}.mark_od_metrics.txt",
    benchmark:
        "benchmarks/mark_optical_duplicates_control/{control}.txt"
    log:
        "logs/mark_optical_duplicates_control/{control}.log"
    params:
        fasta=config["fasta"],
        tmpdir = lambda wildcards: f"{wildcards.control}_ctrl_tmp2/",
    resources:
        mem_mb=30000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=480, increment=60*4, label="mark_optical_duplicates_control"),
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
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=480, increment=60*4, label="mark_read_bundles_control"),
    threads: 36
    group:
        "mark_read_bundles_control"
    shell:
        """
        bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
        samtools index {output.outbam} || exit 1
        """


rule dilute_normal:
    # Produces the matched normal that variant calling consumes (-A track).
    #   Undiluted NanoSeq -> randomreadinbundle dedup of the read-bundled BAM.
    #   Standard WGS      -> use the bulk directly; coordinate-sort + markdup only
    #                        if it is not already sorted (production BAMs are).
    input:
        ctrl_bam=get_normal_prep_input,
    output:
        diluted_control="control_duplex/{control}.diluted.ctrl.bam",
        diluted_control_bai="control_duplex/{control}.diluted.ctrl.bam.bai",
    benchmark:
        "benchmarks/dilute_normal/{control}.txt"
    log:
        "logs/dilute_normal/{control}.log"
    params:
        dilution_factor=0.1,
        control_type=lambda wildcards: get_control_type(wildcards.control),
        wgs_label=CONTROL_TYPE_WGS,
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: _get_dynamic_runtime(attempt, basetime=240, increment=120, label="dilute_normal"),
    threads: 8
    conda:
        workflow.source_path("../envs/nanoseq_snakemake.yml")
    group:
        "dilute_normal"
    shell:
        """
        PATH=$PATH:$PWD/bin/

        if [ "{params.control_type}" = "{params.wgs_label}" ]; then
            # Standard bulk WGS normal: do NOT read-bundle. Use as the -A track
            # directly, ensuring it is coordinate-sorted + duplicate-marked + indexed.
            sortorder=$(samtools view -H {input.ctrl_bam} | awk -F 'SO:' '/^@HD/{{print $2}}' | cut -f1)
            if [ "$sortorder" = "coordinate" ]; then
                echo "Control {wildcards.control} is coordinate-sorted bulk WGS; using directly (assumed duplicate-marked)."
                ln -sf "$(readlink -f {input.ctrl_bam})" {output.diluted_control}
                samtools index {output.diluted_control} {output.diluted_control_bai}
            else
                echo "Control {wildcards.control} bulk WGS not coordinate-sorted; sorting + marking duplicates."
                samtools sort -@ {threads} -o {output.diluted_control}.tmp.sorted.bam {input.ctrl_bam}
                samtools markdup -@ {threads} {output.diluted_control}.tmp.sorted.bam {output.diluted_control}
                rm -f {output.diluted_control}.tmp.sorted.bam
                samtools index {output.diluted_control} {output.diluted_control_bai}
            fi
        else
            # Undiluted NanoSeq library normal: keep one read-pair per read bundle.
            randomreadinbundle -I {input.ctrl_bam} -O {output.diluted_control} && \
            samtools index {output.diluted_control} {output.diluted_control_bai}
        fi
        """
