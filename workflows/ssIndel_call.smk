# WITH SINGLE-STRAND INDEL CALLING MOD
ssIndelPath_wJob=indelPath+"/{job}"
rule indelCall_per_partition_w_ssIndel:
    input:
        dsa                     =   rules.dsa_bed_per_partition.output.coverage,
        indir                   =   rules.coverage_histogram_controlBam.output.runNanoSeqDir,
        duplex_bam              =   rules.mark_read_bundles.output.outbam,
        ctrl_bam                =   check_ctrl_bam_exists,
        coverage                =   rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
        job                     =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",
        nfiles                  =   rules.start_indelCall.output.nfiles,
        argsJson                =   rules.start_indelCall.output.argsJson,
    output:
        indel_bed               =   ssIndelPath_wJob+".indel.w_ssIndel.bed.gz",
        indel_vcf               =   ssIndelPath_wJob+".indel.w_ssIndel.vcf.gz",
        indel_filtered_vcf      =   ssIndelPath_wJob+".indel.w_ssIndel.filtered.vcf.gz",
        indel_filtered_vcf_tbi  =   ssIndelPath_wJob+".indel.w_ssIndel.filtered.vcf.gz.tbi",
        doneFile                =   ssIndelPath_wJob+".w_ssIndel.done",
    benchmark:
        "benchmarks/indelCall_per_partition_w_ssIndel/{sample}.{job}.txt"
    log:
        "logs/indelCall_per_partition_w_ssIndel/{sample}.{job}.log"
    params: # a preset for the ultrashear or covaris libraries
        # DEFAULTS FROM NANOSEQ SOFTWARE `runNanoSeq.py`
        fasta                   =   config["fasta"],
        max_reads_bundle        =   2,
        min_as_xs               =   50,
        min_normal_coverage     =    15,
        max_bulk_vaf            =   0.01,
        ## FROM GILAD'S 2025 CEREBELLUM PAPER: https://www.biorxiv.org/content/10.1101/2025.09.29.679392v1.full.pdf 
        max_frac_clips          =   0,
        trim_from_3p_at_pos     =   135,
        trim_from_5p_at_pos     =   10,
    resources:
        mem_mb                  =   20000,
        # runtime=120,
        # runtime                 =   480,
        runtime                 =   60*15,
    threads: 1
    conda:
        workflow.source_path("../envs/nanoseq_snakemake.yml")
    shell:
        """
        # set path to file
        PATH=$PATH:$PWD/bin/
        PATH=$PATH:$PWD/perl/

        # navigate to output directory
        cd {input.indir}

        echo -e "starting indel call for job #{wildcards.job}"
        # indel call step 1
        # ./tmpNanoSeq/indel/{wildcards.job}.indel.bed.gz
        indelCaller_step1.withSingleStrandIndel.pl \
        -o ../{output.indel_bed} \
        -rb {params.max_reads_bundle} \
        -t3 {params.trim_from_3p_at_pos} \
        -t5 {params.trim_from_5p_at_pos} \
        -mc {params.min_normal_coverage} \
        -vaf {params.max_bulk_vaf} \
        -a {params.min_as_xs} \
        -c {params.max_frac_clips} \
        ../{input.dsa} &&
        # indel call step 2
        indelCaller_step2.pl \
        -t \
        -o ./tmpNanoSeq/indel/{wildcards.job}.indel \
        -r {params.fasta} \
        -b ../{input.duplex_bam} \
        ../{output.indel_bed} &&
        # indel call step 3
        indelCaller_step3.R \
        {params.fasta} \
        ../{output.indel_vcf} \
        ../{input.ctrl_bam} \
        {params.max_bulk_vaf} &&
        # touch
        touch ../{output.doneFile}
        echo -e "single-stranded indel calling job {wildcards.job} is done"
        """
