# function to check if control BAMs exist

CHROMS = [i+1 for i in range(24)]

# REWRITE THIS WROFKLOW TO SPECIFY THOUSANDS OF CHUNKS AND INSTEAD ASSESS 10KB INTERVALS AT A TIME. 

# I need to give this a run manually in a separate directory outside of Snakefile s thatI know what the workflow steps are like

# for this analysis, we will split BAMs by chromosome to conduct variant calling
INTERVALS, = glob_wildcards("intervals/{interval}.intervals.list")
n_jobs_partitioned = len(INTERVALS)
jobs_partitioned = list(range(n_jobs_partitioned+1))[1::]

# MAYBE A FUNCTON THAT RETURNS THE INDEX OF THE INTERVAL -- FEED BOTH THE INDEX AND THE INTERVAL INTO THE `dsa` JOB
# def return_output (wildcards,jobs):

# if exists("input_duplex_fq/{sample}.R1.fastq.gz")
def check_ctrl_bam_exists(wildcards):
    file_path = f"control_bam/{wildcards.sample}.ctrl.bam"
    if os.path.exists(file_path):
        return(file_path)
    else:
        return(None)

# step to downsample the BAM file such that we have just 1 read per read bundle

# # before this step make sure the following perl packages are installed
# # perl -MCPAN -e shell
# # install File::Which
# # install Capture::Tiny

if config["dilute_ctrl_bam"]:
    rule dilute_normal:
        input:
            ctrl_bam=check_ctrl_bam_exists,
        output:
            diluted_control="control_bam/{sample}.diluted.ctrl.bam",
            diluted_control_bai="control_bam/{sample}.diluted.ctrl.bam.bai",
        benchmark:
            "benchmarks/dilute_normal/{sample}.txt"
        log:
            "logs/dilute_normal/{sample}.log"
        params: # a preset for the ultrashear or covaris libraries
            dilution_factor=0.1, # 0.1 means 10%, 0.01 means 1%
        resources:
            mem_mb=10000,
            runtime=240,
        threads: 1
        conda:
            "../envs/nanoseq_snakemake.yml"
        shell:
            """
            PATH=$PATH:$PWD/bin/
            randomreadinbundle -I {input.ctrl_bam} -O {output.diluted_control} && \
            samtools index {output.diluted_control}
            """
else:
    rule dilute_normal:
        input:
            ctrl_bam=check_ctrl_bam_exists,
        output:
            diluted_control="control_bam/{sample}.diluted.ctrl.bam",
            diluted_control_bai="control_bam/{sample}.diluted.ctrl.bam.bai",
        benchmark:
            "benchmarks/dilute_normal/{sample}.txt"
        log:
            "logs/dilute_normal/{sample}.log"
        params: # a preset for the ultrashear or covaris libraries
            dilution_factor=0.1, # 0.1 means 10%, 0.01 means 1%
        resources:
            mem_mb=10000,
            runtime=240,
        threads: 1
        conda:
            "../envs/nanoseq_snakemake.yml"
        shell:
            """
            cd control_bam
            a=$(basename {input.ctrl_bam})
            b=$(basename {output.diluted_control})
            ln -s $a $b
            ln -s $a.bai $b.bai
            """

rule check_efficiency:
    input:
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=rules.dilute_normal.output.diluted_control,
    output:
        read_bundles="efficiency/{sample}.RBs",
        read_bundles_gc_inserts="efficiency/{sample}.RBs.GC_inserts.tsv",
        read_bundles_pdf="efficiency/{sample}.RBs.pdf",
        efficiency_stats="efficiency/{sample}.tsv",
    benchmark:
        "benchmarks/check_efficiency/{sample}.txt"
    log:
        "logs/check_efficiency/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
    resources:
        mem_mb=20000,
        runtime=240,
    threads: 20
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        PATH=$PATH:$PWD/bin/:$PWD/perl/
        perl $PWD/perl/efficiency_nanoseq.pl \
        -t {threads} \
        -duplex {input.duplex_bam} \
        -dedup {input.ctrl_bam} \
        -out efficiency/{wildcards.sample} \
        -ref {params.fasta} 
        """
        # 2> {logs} 1> {logs}


# efficiency_nanoseq.pl \
# -dedup ../20250630_analyis_parallelized/control_bam/S11239_BA9_NeuN_0_05fmol_CKDL250013506-1A_22T3K2LT4_L8.ctrl.bam \
# -duplex ../20250630_analyis_parallelized/readBundle_duplex/S11239_BA9_NeuN_0_05fmol_CKDL250013506-1A_22T3K2LT4_L8.filtered.bam \
# -ref /n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh37d5/human_g1k_v37_decoy.fasta \
# -out test_efficiency -t 12
        
# specify chromosomes that must finish first...go from 1 to 24
rule coverage_histogram_controlBam:
    input:
        duplex_bam=rules.check_efficiency.input.duplex_bam,
        ctrl_bam=rules.dilute_normal.output.diluted_control,
    output:
        # coverage="cov/gIntervals.dat"
        runNanoSeqDir=directory("{sample}.runNanoSeq"),
        coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",chroms=CHROMS,allow_missing=True),
    benchmark:
        "benchmarks/coverage_histogram_controlBam/{sample}.txt"
    log:
        "logs/coverage_histogram_controlBam/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 10
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        
        module load gcc/9.2.0 samtools/1.14

        # newDir="{wildcards.sample}.runNanoSeq/"

        mkdir -p {output.runNanoSeqDir}

        cd {output.runNanoSeqDir}

        runNanoSeq.py \
        -t {threads} \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -R {params.fasta} \
        cov \
        -Q 0 \
        --exclude "MT,GL%,NC_%,hs37d5"

        """

rule partition_coverage:
    input:
        indir=rules.coverage_histogram_controlBam.output.runNanoSeqDir,
        duplex_bam=rules.check_efficiency.input.duplex_bam,
        ctrl_bam=rules.dilute_normal.output.diluted_control,
    output:
        coverage="{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
    benchmark:
        "benchmarks/partition_coverage/{sample}.txt"
    log:
        "logs/partition_coverage/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        n_partitions=60,
        jobs=n_jobs_partitioned,
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 10
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        PATH=$PATH:$PWD/bin/

        cd {input.indir}

        runNanoSeq.py \
        -t 1 \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -R {params.fasta} \
        part \
        -n {params.jobs}

        """

rule list_intervals:
    input:
        intvl=expand("intervals/{interval}.intervals.list",interval=INTERVALS)
    output:
        intvl_list=temp(".intervals.txt")
    benchmark:
        "benchmarks/list_intervals/all.txt"
    log:
        "logs/list_intervals/all.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        snp=SNP,
        noise=NOISE,
        # jobs=n_jobs,
    resources:
        mem_mb=5000,
        runtime=20,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        ls -1 {input.intvl} | sed "s/_/\t/g" | sort -k1,1 -k2,2n | sed "s/\t/_/g" > {output.intvl_list}
        """

rule start_dsa:
    #  next time, modify to add `nfiles` files
    # next time, add the args.json here...
    input:
        intvl_list=rules.list_intervals.output.intvl_list,
    output:
        jobDone=expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",job=jobs_partitioned,allow_missing=True), # job=n_jobs_partitioned,allow_missing=True),
        jobToIntvl="{sample}.runNanoSeq/tmpNanoSeq/dsa/job_to_intvl.txt",
    params:
        jobs=jobs_partitioned
    benchmark:
        "benchmarks/start_dsa/{sample}.txt"
    log:
        "logs/start_dsa/{sample}.log"
    resources:
        mem_mb=5000,
        runtime=20,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        touch {output.jobToIntvl}
        for i in {params.jobs}; do
            intvl=$(head -$i {input.intvl_list} | tail -1)
            echo -e "$intvl" > {wildcards.sample}.runNanoSeq/tmpNanoSeq/dsa/$i.start
            echo -e "$i\t$intvl" >> {output.jobToIntvl}
        done
        """

rule add_dsa_args:
    # modified to add `nfiles` files
    input:
        intvl_list=rules.list_intervals.output.intvl_list,
    output:
        nfiles="{sample}.runNanoSeq/tmpNanoSeq/dsa/nfiles",
        argsJson="{sample}.runNanoSeq/tmpNanoSeq/dsa/args.json",
    params:
        jobs=jobs_partitioned
    benchmark:
        "benchmarks/add_dsa_args/{sample}.txt"
    log:
        "logs/add_dsa_args/{sample}.log"
    resources:
        mem_mb=5000,
        runtime=10,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """        

# {"out": ".", "index": null, "max_index": null, "threads": 300, "ref": "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh37d5/human_g1k_v37_decoy.fasta", "normal": "../control_bam/S11239_BA9_NeuN_0_05fmol_CKDL250013506-1A_22T3K2LT4_L8.ctrl.bam", "duplex": "../readBundle_duplex/S11239_BA9_NeuN_0_05fmol_CKDL250013506-1A_22T3K2LT4_L8.filtered.bam", "subcommand": "dsa", "snp": "../test/SNP.sorted.bed.gz", "mask": "../test/NOISE.sorted.bed.gz", "d": 2, "q": 30, "no_test": false}
rule dsa_bed_per_partition:
    input:
        indir                                   =   rules.coverage_histogram_controlBam.output.runNanoSeqDir,
        duplex_bam                              =   rules.check_efficiency.input.duplex_bam,
        ctrl_bam                                =   rules.dilute_normal.output.diluted_control,
        coverage                                =   rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
        job                                     =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start", # something wrong with this job not existing...
        allJobs                                 =   rules.start_dsa.output.jobToIntvl,
        nfiles_dsa_arg                          =   rules.add_dsa_args.output.nfiles,
        argsJson                                =   rules.add_dsa_args.output.argsJson,
    output:
       coverage                                 =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.dsa.bed.gz", # job=n_jobs_partitioned,allow_missing=True),
       jobDone                                  =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.done", # job=n_jobs_partitioned,allow_missing=True),
    benchmark:
        "benchmarks/dsa_bed_per_partition/{sample}.{job}.txt"
    log:
        "logs/dsa_bed_per_partition/{sample}.{job}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta                                   =   config["fasta"],
        snp                                     =   SNP,
        noise                                   =   NOISE,
    resources:
        mem_mb                                  =   5000,
        runtime                                 =   60,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        # set path to file
        PATH=$PATH:$PWD/bin/

        # get line number of the interval -- this will be our chunk
        jobIndex={wildcards.job}
        intvl=$(head -1 {input.job})
        intvl=$(basename $intvl)
        intvl=${{intvl%.intervals*}}

        # format the interval
        chrom=$(echo -e $intvl | cut -d'_' -f1)
        startPos=$(echo -e $intvl | cut -d'_' -f2)
        endPos=$(echo -e $intvl | cut -d'_' -f3)

        # navigate to output directory
        cd {input.indir}
        
        echo -e "running dsa..."
        dsa \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -C {params.snp} \
        -D {params.noise} \
        -R {params.fasta} \
        -d 2 \
        -Q 30 \
        -M 0 \
        -r $chrom \
        -b $startPos \
        -e $endPos \
        -O ./tmpNanoSeq/dsa/$jobIndex.dsa.bed &&
        touch ./tmpNanoSeq/dsa/$jobIndex.done
        echo -e "Job {wildcards.job} over $intvl is done"
        # touch ./tmpNanoSeq/dsa/$intvl.done
        """


# variant calling step
rule start_varCall:
    # modified to add `nfiles` files
    input:
        intvl_list=rules.list_intervals.output.intvl_list,
    output:
        nfiles="{sample}.runNanoSeq/tmpNanoSeq/var/nfiles",
        argsJson="{sample}.runNanoSeq/tmpNanoSeq/var/args.json",
    params:
        jobs=jobs_partitioned
    benchmark:
        "benchmarks/start_varCall/{sample}.txt"
    log:
        "logs/start_varCall/{sample}.log"
    resources:
        mem_mb=5000,
        runtime=2,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """        

varPath="{sample}.runNanoSeq/tmpNanoSeq/var/{job}"
rule varCall_per_partition:
    input:
        dsa                                 =   rules.dsa_bed_per_partition.output.coverage,
        nfiles                              =   rules.start_varCall.output.nfiles,
        indir                               =   rules.coverage_histogram_controlBam.output.runNanoSeqDir,
        duplex_bam                          =   rules.check_efficiency.input.duplex_bam,
        ctrl_bam                            =   rules.dilute_normal.output.diluted_control,
        coverage                            =   rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
        job                                 =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",
    output:
       coverage                             =   varPath+".cov.bed.gz",
       var                                  =   varPath+".var",
       discarded_var                        =   varPath+".discarded_var",
       doneFile                             =   varPath+".done",
    benchmark:
        "benchmarks/varCall_per_partition/{sample}.{job}.txt"
    log:
        "logs/varCall_per_partition/{sample}.{job}.log"
    params: # a preset for the ultrashear or covaris libraries
        min_as_xs                           =   50,
        min_bulk_reads_per_strand           =   3, # 5 formerly, but per Gilad: 3 for non-Nanoseq germline library, set to 0 for undiluted nanoseq library
        min_number_dplx_reads_per_strand    =   2, # toggle this option to go from a4s2 (default) to a2s1 or a6s3. 
        min_fraction_reads_consensus        =   0.9,
        max_frac_reads_w_indel              =   1.0,
        min_cycle_num                       =   80,
        max_num_mismatches                  =   3,
        min_fraction_reads_properPairs      =   0,
        min_consensus_base_quality          =   60,
        read_length_post_5prime_trimming    =   144,
        max_bulk_vaf                        =   0.01,
        max_cycle_num                       =   8,
        min_num_bulkReads_total             =   12,
        max_frac_clips                      =   0.02,
    resources:
        mem_mb                              =   500,
        runtime                             =   5,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        # set path to file
        PATH=$PATH:$PWD/bin/

        # navigate to output directory
        cd {input.indir}

        # run variant caller
        variantcaller \
        -B ./tmpNanoSeq/dsa/{wildcards.job}.dsa.bed.gz \
        -U ./tmpNanoSeq/var/{wildcards.job}.cov.bed \
        -O ./tmpNanoSeq/var/{wildcards.job}.var \
        -D ./tmpNanoSeq/var/{wildcards.job}.discarded_var \
        -a {params.min_as_xs} \
        -b {params.min_bulk_reads_per_strand} \
        -c {params.max_frac_clips} \
        -d {params.min_number_dplx_reads_per_strand} \
        -f {params.min_fraction_reads_consensus} \
        -i {params.max_frac_reads_w_indel} \
        -m {params.min_cycle_num} \
        -n {params.max_num_mismatches} \
        -p {params.min_fraction_reads_properPairs} \
        -q {params.min_consensus_base_quality} \
        -r {params.read_length_post_5prime_trimming} \
        -v {params.max_bulk_vaf} \
        -x {params.max_cycle_num} \
        -z {params.min_num_bulkReads_total} &&
        touch ../{output.doneFile}
        echo -e "variant calling for {wildcards.job} is done"
        """

# Indel calling
indelPath="{sample}.runNanoSeq/tmpNanoSeq/indel"
rule start_indelCall:
    # modified to add `nfiles` files
    input:
        intvl_list              =   rules.list_intervals.output.intvl_list,
    output:
        nfiles                  =   indelPath+"/nfiles",
        argsJson                =   indelPath+"/args.json",
    params:
        jobs                    =   jobs_partitioned
    benchmark:
        "benchmarks/start_indelCall/{sample}.txt"
    log:
        "logs/start_indelCall/{sample}.log"
    resources:
        mem_mb                  =   100,
        # runtime                 =   10,
        runtime                 =   2,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """        

indelPath_wJob=indelPath+"/{job}"
rule indelCall_per_partition:
    input:
        dsa                     =   rules.dsa_bed_per_partition.output.coverage,
        indir                   =   rules.coverage_histogram_controlBam.output.runNanoSeqDir,
        duplex_bam              =   rules.check_efficiency.input.duplex_bam,
        ctrl_bam                =   rules.dilute_normal.output.diluted_control,
        coverage                =   rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
        job                     =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",
    output:
        indel_bed               =   indelPath_wJob+".indel.bed.gz",
        indel_vcf               =   indelPath_wJob+".indel.vcf.gz",
        indel_filtered_vcf      =   indelPath_wJob+".indel.filtered.vcf.gz",
        indel_filtered_vcf_tbi  =   indelPath_wJob+".indel.filtered.vcf.gz.tbi",
        doneFile                =   indelPath_wJob+".done",
    benchmark:
        "benchmarks/indelCall_per_partition/{sample}.{job}.txt"
    log:
        "logs/indelCall_per_partition/{sample}.{job}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta                   =   config["fasta"],
        max_reads_bundle        =   2,
        trim_from_3p_at_pos     =   135,
        trim_from_5p_at_pos     =   10,
        min_as_xs               =   50,
        max_frac_clips          =   0.02,
        min_normal_coverage     =    16,
        max_bulk_vaf            =   0.2, # this should have been 0.2!
    resources:
        mem_mb                  =   20000,
        runtime                 =   5,
    threads: 1
    conda:
        "../envs/nanoseq_snakemake.yml"
    shell:
        """
        # set path to file
        PATH=$PATH:$PWD/bin/

        # navigate to output directory
        cd {input.indir}

        echo -e "starting indel call for job #{wildcards.job}"
        # indel call step 1
        # ./tmpNanoSeq/indel/{wildcards.job}.indel.bed.gz
        indelCaller_step1.pl \
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
        echo -e "indel calling job {wildcards.job} is done"

        """

# I should add options to specify the final sample name of the files we are analyzing