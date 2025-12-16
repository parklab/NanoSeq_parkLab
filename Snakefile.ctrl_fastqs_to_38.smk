# shell.prefix("VERSION=202503.01; \
# export LD_LIBRARY_PATH=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-$VERSION/lib:$LD_LIBRARY_PATH; \
# export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990; \
# export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-$VERESION; \
# export PATH=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-$VERSION/bin:$PATH; \
# module load gcc/14.2.0 ; \
# module load bcftools/1.21; \
# module load java; \
# export PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab/bin;")

rule ctrl_extract_tags:
    input:
        expand("../20251201_doubleSidedSelectionTest_ctrl/S11239BA9NeuNpos_p002D05_CKDL250023474-1A_22WVNCLT4_L7_{pe}.fq.gz",pe=PAIRED_ENDS,allow_missing=True)

    output:
        expand("extracted_tags_ctrl/{sample}_{pe}.fq.gz",pe=PAIRED_ENDS,allow_missing=True)
    benchmark:
        "benchmarks/extract_tags_ctrl/{sample}.txt"
    log:
        "logs/extract_tags_ctrl/{sample}.log"
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
        python ../python/extract_tags.py \
        -a {input[0]} \
        -b {input[1]} \
        -c {output[0]} \
        -d {output[1]} \
        -m {params.bases_trim} \
        -s {params.bases_skip} \
        -l {params.read_len}
        """


rule realign_to_hg38:
    input:
        r1  = "extracted_tags_ctrl/{sample}_1.fq.gz",
        r2  = "extracted_tags_ctrl/{sample}_2.fq.gz",
        ref = config["fasta"] + ".bwt"  # ensures reference is indexed
    output:
        # outbam="control_bam/{sample}.ctrl.bam",  
        outsam="control_bam/{sample}.ctrl.sam",
        # outbai="control_bam/{sample}.ctrl.bam.bai",
    log:
        "logs/realign_control_bam/{sample}.log"
    resources:
        mem_mb=50000,
        runtime=240,
    threads: 20
    shell:
        """
        set -euo pipefail
        sentieon \
        bwa mem \
        -t {threads} \
        -K 10000000 \
        -C \
        {config[fasta]} \
        {input.r1} {input.r2} > {output.outsam} || exit 1
        """


rule ctrl_prepare_RcMcOd_tags:
    input:
        rules.realign_to_hg38.output.outsam
    output:
        tmpdir=temp(directory("{sample}_ctrl_tmp/")),
        outbam="control_bam/{sample}.od.bam"
    benchmark:
        "benchmarks/prepare_RcMcOd_tags/{sample}.ctrl.txt"
    log:
        "logs/prepare_RcMcOd_tags/{sample}.ctrl.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads: 12
    shell:
        """
        ulimit -c unlimited
        ulimit -n 8192
        tmpdir={wildcards.sample}_ctrl_tmp
        mkdir -p $tmpdir
        echo -e "{wildcards.sample}" >> {log}
        bamsormadup inputformat=sam rcsupport=1 threads=1 tmpfile=$tmpdir/{wildcards.sample} < {input} > {output.outbam}
        """

rule ctrl_mark_optical_duplicates:
    input:
        rules.ctrl_prepare_RcMcOd_tags.output.outbam
    output:
        outbam="control_bam/{sample}.marked_od.ctrl.bam",
        metrics="control_bam/{sample}.mark_od_metrics.txt",
    benchmark:
        "benchmarks/mark_optical_duplicates/{sample}.ctrl.txt"
    log:
        "logs/mark_optical_duplicates/{sample}.ctrl.log"
    params:
        fasta=config["fasta"],
        outduplicates="control_bam/{sample}.duplicates.ctrl.bam",
        tmpdir=temp(directory("{sample}_ctrl_tmp2/")),
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


rule control_bam_add_readbundle:
    input:
        rules.ctrl_mark_optical_duplicates.output.outbam,
    output:
        RG_control="control_bam/{sample}.ctrl.bam",
        RG_control_bai="control_bam/{sample}.ctrl.bam.bai",
    benchmark:
        "benchmarks/control_bam_add_RG/{sample}.txt"
    log:
        "logs/control_bam_add_RG/{sample}.log"
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 1
    shell: 
        """
        PATH=$PATH:$PWD/bin/
        bamaddreadbundles -I {input} -O {output.RG_control} || exit 1
        samtools index {output.RG_control} || exit 1
        """
