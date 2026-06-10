# here, we associate each sample with its control, and then run the entire NanoSeq pipeline for each sample-control pair.

def get_control_bam(wildcards):
    inTable = pd.read_csv("input.tsv",sep="\t",header=0)
    control = inTable.loc[inTable["TestFastqID"]==wildcards.sample,"ControlFastqID"].values[0]
    control_bam = f"control_duplex/{control}.diluted.ctrl.bam"
    return control_bam


# Retry-aware runtime helpers. Each rule's `resources.runtime` should be a
# lambda over (wildcards, attempt) so that Snakemake's `--retries N` actually
# requests a larger wall-clock from SLURM on subsequent attempts.
def _get_dynamic_runtime(attempt, basetime, increment, label="rule"):
    total = int(basetime + increment * (attempt - 1))
    print(f"[runtime] attempt {attempt} for {label}: requesting {total} min (base={basetime}, +{increment}/attempt)")
    return total

def get_dynamic_runtime_efficiency(attempt):
    return _get_dynamic_runtime(attempt, basetime=int(60*0.5), increment=60, label="check_efficiency")

def get_dynamic_runtime_coverage(attempt):
    return _get_dynamic_runtime(attempt, basetime=int(60*1/6), increment=60, label="coverage_histogram_controlBam")

def get_dynamic_runtime_diag_coverage(attempt):
    return _get_dynamic_runtime(attempt, basetime=60*2, increment=60, label="run_diagnostic_coverage")

def get_dynamic_runtime_dsa(attempt):
    return _get_dynamic_runtime(attempt, basetime=60*10, increment=60*4, label="dsa_bed_per_partition")

def get_dynamic_runtime_varcall(attempt):
    return _get_dynamic_runtime(attempt, basetime=5*10, increment=60, label="varCall_per_partition")

def get_dynamic_runtime_indelcall(attempt):
    return _get_dynamic_runtime(attempt, basetime=60*15, increment=60*4, label="indelCall_per_partition")

def get_dynamic_runtime_post(attempt):
    return _get_dynamic_runtime(attempt, basetime=120, increment=120, label="post")

rule check_efficiency:
    input:
        duplex_bam=rules.mark_read_bundles.output.outbam,
        # duplex_bam="readBundle_duplex/{sample}.filtered.bam",
        ctrl_bam=get_control_bam,
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
        runtime=lambda wildcards, attempt: get_dynamic_runtime_efficiency(attempt),
    threads: 20
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "check_efficiency"
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

# specify chromosomes that must finish first...go from 1 to 24
rule coverage_histogram_controlBam:
    input:
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=get_control_bam,
    output:
        # coverage="cov/gIntervals.dat"
        # runNanoSeqDir=directory("{sample}.runNanoSeq"),
        coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",chroms=CHROMS,allow_missing=True),
    benchmark:
        "benchmarks/coverage_histogram_controlBam/{sample}.txt"
    log:
        "logs/coverage_histogram_controlBam/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        runNanoSeqDir = lambda wildcards: f"{wildcards.sample}.runNanoSeq",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: get_dynamic_runtime_coverage(attempt),
    threads: 10
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "coverage_histogram_controlBam"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        
        module load gcc/9.2.0 samtools/1.14

        # newDir="{wildcards.sample}.runNanoSeq/"

        mkdir -p {params.runNanoSeqDir}

        cd {params.runNanoSeqDir}

        runNanoSeq.py \
        -t {threads} \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -R {params.fasta} \
        cov \
        -Q 0 \
        --exclude "MT,GL%,NC_%,hs37d5"

        """

# we can use the coverage diagnostic to figure out what percent of the reference genome was assessed. 
rule run_diagnostic_coverage:
    input:
        coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",chroms=CHROMS,allow_missing=True),
    output:
        coverage="{sample}.runNanoSeq/tmpNanoSeq/diagnostic/{sample}.tsv",
        coverage_table="{sample}.runNanoSeq/tmpNanoSeq/diagnostic/{sample}.bulk_coverage_diagnostic.logistic_pred.txt",
        coverage_model="{sample}.runNanoSeq/tmpNanoSeq/diagnostic/{sample}.bulk_coverage_diagnostic.logistic_fit.rds",
        threshold="{sample}.runNanoSeq/tmpNanoSeq/diagnostic/{sample}.bulk_coverage_diagnostic.threshold.txt",
    benchmark:
        "benchmarks/run_diagnostic_coverage/{sample}.txt"
    log:
        "logs/run_diagnostic_coverage/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        maxCovAssess=20, # next time, set to 30 instead of 40. more time-efficient, and 30 is sufficient for estimating the coverage threshold.
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: get_dynamic_runtime_diag_coverage(attempt),
    threads: 10
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "run_diagnostic_coverage"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        
        nanoseqDir={wildcards.sample}.runNanoSeq/tmpNanoSeq

        for i in $(seq 0 1 {params.maxCovAssess}); do 
            for j in $nanoseqDir/cov/*.gz; do 
                zcat $j | awk -v threshold=$i '$NF/($3 - $2) >= threshold' | grep -c "^" | paste <(echo -e "$j\t$i") -
            done
        done > {output.coverage}

        # Rscript that fits a logistic curve to the coverage diagnostic and estimates the minimum coverage to use.
        outputName=$nanoseqDir/diagnostic/{wildcards.sample}.bulk_coverage_diagnostic
        Rscript $PWD/R/estimate_bulk_min_cov.R \
        {output.coverage} \
        $outputName

        """

rule partition_coverage:
    input:
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=get_control_bam,
        histCoverage=expand("{sample}.runNanoSeq/tmpNanoSeq/cov/{chroms}.done",chroms=CHROMS,allow_missing=True),
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
        indir=lambda wildcards: f"{wildcards.sample}.runNanoSeq",
    resources:
        mem_mb=10000,
        # runtime=240,
        runtime=5, # next time, set to 5 minutes instead of 4 hours. more time-efficient
    threads: 10
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "partition_coverage"
    shell:
        """
        PATH=$PATH:$PWD/bin/

        cd {params.indir}

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
        intvl_list=".intervals.txt"
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
        # runtime=20,
        runtime=20, # next time, set to 20 minutes instead of 20 minutes. more time-efficient
    threads: 1
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "list_intervals"
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
        # jobDone=temp(expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",job=jobs_partitioned,allow_missing=True)), # job=n_jobs_partitioned,allow_missing=True),
        jobDone=expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",job=jobs_partitioned,allow_missing=True), # job=n_jobs_partitioned,allow_missing=True),
        jobToIntvl="{sample}.runNanoSeq/tmpNanoSeq/dsa/job_to_intvl.txt",
        # nfiles="{sample}.runNanoSeq/tmpNanoSeq/dsa/nfiles"
    params:
        jobs=jobs_partitioned
    # benchmark:
    #     "benchmarks/start_dsa/all.txt"
    # log:
    #     "logs/start_dsa/all.log"
    resources:
        mem_mb=5000,
        # runtime=20,
        runtime=5, # next time, set to 5 minutes instead of 20 minutes. more time-efficient
    threads: 1
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "start_dsa"
    shell:
        """
        touch {output.jobToIntvl}
        for i in {params.jobs}; do
            intvl=$(head -$i {input.intvl_list} | tail -1)
            echo -e "$intvl" > {wildcards.sample}.runNanoSeq/tmpNanoSeq/dsa/$i.start
            echo -e "$i\t$intvl" >> {output.jobToIntvl}
        done
        """
        #grep -c "^" {input.intvl_list} > {wildcards.sample}.runNanoSeq/tmpNanoSeq/dsa/nfiles

rule add_dsa_args:
    # modified to add `nfiles` files
    input:
        intvl_list=rules.list_intervals.output.intvl_list,
        jobToIntvl=rules.start_dsa.output.jobToIntvl,
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
        # runtime=10,
        runtime=2, # next time, set to 2 minutes instead of 10 minutes. more time-efficient
    threads: 1
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "add_dsa_args"
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """        

# HOW DO WE GET THE NANO-SEQ RULES TO RECOGNIZE THAT WE CREATED A JOB WITH THE CORRESPONDING NUMBER?
# checkpoints: "A more practical example building on the previous one is a clustering process with an unknown number of clusters for different samples, where each cluster shall be saved into a separate file. In this example the clusters are being processed by an intermediate rule before being aggregated:..."
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution 
# checkpoint check_start_dsa:
# SEPT 21, 2025: NOT NECESSARY

# {"out": ".", "index": null, "max_index": null, "threads": 300, "ref": "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh37d5/human_g1k_v37_decoy.fasta", "normal": "../control_bam/S11239_BA9_NeuN_0_05fmol_CKDL250013506-1A_22T3K2LT4_L8.ctrl.bam", "duplex": "../readBundle_duplex/S11239_BA9_NeuN_0_05fmol_CKDL250013506-1A_22T3K2LT4_L8.filtered.bam", "subcommand": "dsa", "snp": "../test/SNP.sorted.bed.gz", "mask": "../test/NOISE.sorted.bed.gz", "d": 2, "q": 30, "no_test": false}
rule dsa_bed_per_partition:
    input:
        duplex_bam                              =   rules.mark_read_bundles.output.outbam,
        ctrl_bam                                =   get_control_bam,
        # ctrl_bam                                =   rules.dilute_normal.output.diluted_control,
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
        indir                                   =   lambda wildcards: f"{wildcards.sample}.runNanoSeq",
        # jobs=n_jobs,
    resources:
        mem_mb                                  =   5000,
        runtime                                 =   lambda wildcards, attempt: get_dynamic_runtime_dsa(attempt),
    threads: 1
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "dsa_bed_per_partition"
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
        cd {params.indir}
        
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
        -O ./tmpNanoSeq/dsa/$jobIndex.dsa.bed
        if [ -f ./tmpNanoSeq/dsa/$jobIndex.dsa.bed.gz ]; then
            touch ./tmpNanoSeq/dsa/$jobIndex.done
            echo -e "Job {wildcards.job} over $intvl is done"
        else
            echo -e "DSA failed for job {wildcards.job} over $intvl"
            exit 1
        fi
        """
        # jobIndex=$(grep -n -o {input.intvl} {input.intvl_list})
        # jobIndex=$(grep -s {wildcards.interval} {input.allJobs} | cut -f1)
        # inputJob={wildcards.sample}.runNanoSeq/tmpNanoSeq/dsa/$jobIndex.start
        # intvl=$(head -1 $inputJob)


# variant calling step
rule start_varCall:
    # modified to add `nfiles` files
    input:
        intvl_list=rules.list_intervals.output.intvl_list,
        nfiles_dsa_arg                          =   rules.add_dsa_args.output.nfiles,
        argsJson                                =   rules.add_dsa_args.output.argsJson,
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
        # runtime=10,
        runtime=2,
    threads: 1
    group:
        "start_varCall"
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
        duplex_bam                          =   rules.mark_read_bundles.output.outbam,
        ctrl_bam                            =   get_control_bam,
        coverage                            =   rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
        job                                 =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",
        min_num_bulkReads_total_file        =   rules.run_diagnostic_coverage.output.threshold, # formerly 12 -- empirically assessed # should probably be 6, but reducing it all the way to 0 so that we can use external SNP filters
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
        min_as_xs                           =   50, # as a last resort, co-pilot suggested dropping this to 30
        min_bulk_reads_per_strand           =   0, # 5 formerly, but per Gilad: 3 for non-Nanoseq germline library, set to 0 for undiluted nanoseq library
        max_frac_clips                      =   0, # formerly 0.02, but the Sanger default is 0
        min_number_dplx_reads_per_strand    =   2, # toggle this option to go from a4s2 (default) to a2s1 or a6s3. 
        min_fraction_reads_consensus        =   0.9,
        max_frac_reads_w_indel              =   1.0,
        min_cycle_num                       =   8, # formerly 80, which was too high
        max_num_mismatches                  =   3, # co-pilot suggested raising this to 6
        min_fraction_reads_properPairs      =   0,
        min_consensus_base_quality          =   60,
        read_length_post_5prime_trimming    =   144,
        max_bulk_vaf                        =   0.02, # 0.01, # 0.01 may be way too strict
        max_cycle_num                       =   8,
        # min_num_bulkReads_total             =   6, # formerly 12 -- empirically assessed # should probably be 6, but reducing it all the way to 0 so that we can use external SNP filters
        indir                                   =   lambda wildcards: f"{wildcards.sample}.runNanoSeq",
    resources:
        mem_mb                              =   5000,
        runtime                             =   lambda wildcards, attempt: get_dynamic_runtime_varcall(attempt),
    threads: 1
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "varCall_per_partition"
    shell:
        """
        # set path to file
        PATH=$PATH:$PWD/bin/

        # # get the threshold for min_num_bulkReads_total from the output of the coverage diagnostic
        min_num_bulkReads_total=$(cut -f1 {input.min_num_bulkReads_total_file} | tail -1)
        echo -e "Using a threshold of $min_num_bulkReads_total for min_num_bulkReads_total"

        # navigate to output directory
        cd {params.indir}

        # for segmentation faults
        ulimit -c unlimited
        ulimit -s unlimited

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
        -z $min_num_bulkReads_total &&
        touch ../{output.doneFile}
        echo -e "variant calling for {wildcards.job} is done"
        """
        # -z {params.min_num_bulkReads_total} &&

# Indel calling
indelPath="{sample}.runNanoSeq/tmpNanoSeq/indel"
rule start_indelCall:
    # modified to add `nfiles` files
    input:
        intvl_list              =   rules.list_intervals.output.intvl_list,
        nfiles_dsa_arg                          =   rules.add_dsa_args.output.nfiles,
        argsJson                                =   rules.add_dsa_args.output.argsJson,
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
        mem_mb                  =   5000,
        # runtime                 =   10,
        runtime                 =   10,
    threads: 1
    # conda:
    #     workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "start_indelCall"
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """        

indelPath_wJob=indelPath+"/{job}"
rule indelCall_per_partition:
    input:
        dsa                     =   rules.dsa_bed_per_partition.output.coverage,
        duplex_bam              =   rules.mark_read_bundles.output.outbam,
        ctrl_bam                =   get_control_bam,
        coverage                =   rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
        job                     =   "{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.start",
        nfiles                  =   rules.start_indelCall.output.nfiles,
        argsJson                =   rules.start_indelCall.output.argsJson,
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
        # DEFAULTS FROM NANOSEQ SOFTWARE `runNanoSeq.py`
        fasta                   =   config["fasta"],
        max_reads_bundle        =   2,
        min_as_xs               =   50,
        min_normal_coverage     =    15,
        max_bulk_vaf            =   0.02, # 0.01, # 0.01 may be way too strict
        ## FROM GILAD'S 2025 CEREBELLUM PAPER: https://www.biorxiv.org/content/10.1101/2025.09.29.679392v1.full.pdf 
        max_frac_clips          =   0,
        trim_from_3p_at_pos     =   135,
        trim_from_5p_at_pos     =   10,
        indir                   =   lambda wildcards: f"{wildcards.sample}.runNanoSeq",
        ## FROM THE SCALEUP EXPERIMENT
        # max_frac_clips          =   0.02,
        # max_bulk_vaf            =   0.2,
        # trim_from_3p_at_pos     =   135,
        # trim_from_5p_at_pos     =   10,
        # min_normal_coverage     =    16,
        ## FOR FUTURE?
        # max_bulk_vaf            =   0.1, # to be very strict for future runs...
        # trim_from_3p_at_pos     =   136, # DEFAULT FROM `runNanoSeq.py`
        # trim_from_5p_at_pos     =   8, # DEFAULT FROM `runNanoSeq.py`
    resources:
        mem_mb                  =   20000,
        runtime                 =   lambda wildcards, attempt: get_dynamic_runtime_indelcall(attempt),
    threads: 1
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "indelCall_per_partition"
    shell:
        """
        # set path to file
        PATH=$PATH:$PWD/bin/

        # navigate to output directory
        cd {params.indir}

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

# include: "workflows/ssIndel_call.smk"

# I should add options to specify the final sample name of the files we are analyzing
# post
# postPath="{sample}.runNanoSeq/tmpNanoSeq/post/"
rule post:
    input:
        duplex_bam                          =   rules.mark_read_bundles.output.outbam,
        ctrl_bam                            =   get_control_bam,
        snv_doneFile                        =   expand(indelPath_wJob+".done",job=jobs_partitioned,allow_missing=True),
        indel_doneFile                      =   expand(varPath+".done",job=jobs_partitioned,allow_missing=True),
        nfiles                              =   rules.start_indelCall.output.nfiles,
        argsJson                            =   rules.start_indelCall.output.argsJson,
    output:
        variants                            =   postPath+"results.muts.vcf.gz",
        doneFile                            =   postPath+"1.done",
    benchmark:
        "benchmarks/post/{sample}.txt"
    log:
        "logs/post/{sample}.log"
    params:
        fasta=config["fasta"],
        indir=lambda wildcards: f"{wildcards.sample}.runNanoSeq",
    threads: 10
    resources:
        mem_mb                  =   5000,
        runtime                 =   lambda wildcards, attempt: get_dynamic_runtime_post(attempt),
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "post"
    shell:
        """
        cd {params.indir}
        mkdir -p ./tmpNanoSeq/post
        touch ./tmpNanoSeq/post/args.json
        runNanoSeq.py \
        -t {threads} \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -R {params.fasta} \
        post
        """

# SUMMARY_RESULT_FILES = ["burdens", "callvsqpos", "coverage", "pyrvsmask", "readbundles"]
rule annotate_post_csv_with_interval:
    input:
        variants                            =   postPath+"results.muts.vcf.gz",
    output:
        new_burden                          = postPath+"burdens.annotated.tsv",
        new_callvsqpos                      = postPath+"callvsqpos.annotated.tsv",
        new_coverage                        = postPath+"coverage.annotated.tsv",
        new_pyrvsmask                       = postPath+"pyrvsmask.annotated.tsv",
        new_readbundles                     = postPath+"readbundles.annotated.tsv",
        # annotated_files                     =   expand(postPath+"{summary_results}.annotated.tsv",summary_results=SUMMARY_RESULT_FILES),
        doneFile                            =   postPath+"1.annotate_interval.done",
    benchmark:
        "benchmarks/post/{sample}.txt"
    log:
        "logs/post/{sample}.log"
    params:
        fasta=config["fasta"],
        indir=lambda wildcards: f"{wildcards.sample}.runNanoSeq",
    threads: 10
    resources:
        mem_mb                  =   5000,
        runtime                 =   lambda wildcards, attempt: get_dynamic_runtime_post(attempt),
    conda:
        workflow.source_path("envs/nanoseq_snakemake.yml")
    group:
        "post"
    shell:
        """
        cd {params.indir}
        postDir="./tmpNanoSeq/post"
        varDir="./tmpNanoSeq/var"
        dsaDir="./tmpNanoSeq/dsa"

        new_burdens=$postDir"/burdens.annotated.tsv"
        new_callvsqpos=$postDir"/callvsqpos.annotated.tsv"
        new_coverage=$postDir"/coverage.annotated.tsv"
        new_pyrvsmask=$postDir"/pyrvsmask.annotated.tsv"
        new_readbundles=$postDir"/readbundles.annotated.tsv"
        
        paste <(head -1 $postDir/burdens.csv) <(echo -e "intvl") | sed "s/,/\t/g" > $new_burdens
        paste <(head -1 $postDir/callvsqpos.csv) <(echo -e "intvl") | sed "s/,/\t/g" > $new_callvsqpos
        paste <(head -1 $postDir/coverage.csv) <(echo -e "intvl") | sed "s/,/\t/g" > $new_coverage
        paste <(head -1 $postDir/pyrvsmask.csv) <(echo -e "intvl") | sed "s/,/\t/g" > $new_pyrvsmask
        paste <(head -1 $postDir/readbundles.csv) <(echo -e "intvl") | sed "s/,/\t/g" > $new_readbundles

        for i in $varDir/*.var; do
            intvlChunk=$(basename $i)
            intvlChunk=${{intvlChunk%.var}}
            dsaIntvlChunk=$(cat $dsaDir/$intvlChunk.start)
            grep -s "Burdens" $i | cut -f2- | awk -v chunk=$dsaIntvlChunk 'BEGIN {{FS=OFS="\t"}} {{print $0,chunk}}' >> $new_burdens
            grep -s "CallVsQpos" $i | cut -f2- | awk -v chunk=$dsaIntvlChunk 'BEGIN {{FS=OFS="\t"}} {{print $0,chunk}}' >> $new_callvsqpos
            grep -s "Coverage" $i | cut -f2- | awk -v chunk=$dsaIntvlChunk 'BEGIN {{FS=OFS="\t"}} {{print $0,chunk}}' >> $new_coverage
            grep -s "PyrVsMask" $i | cut -f2- | awk -v chunk=$dsaIntvlChunk 'BEGIN {{FS=OFS="\t"}} {{print $0,chunk}}' >> $new_pyrvsmask
            grep -s "ReadBundles" $i | cut -f2- | awk -v chunk=$dsaIntvlChunk 'BEGIN {{FS=OFS="\t"}} {{print $0,chunk}}' >> $new_readbundles
        done

        touch $postDir"/1.annotate_interval.done"
        """

# rule filter_keep_passing_remove_likely_germline:
#     input:
#         variants                            =   postPath+"results.muts.vcf.gz",
#         indels                              =   postPath+"results.indels.vcf.gz",
#     output:
#         variants                            =   postPath+"results.muts.remove_noise_remove_bulkNM0_keepPass.vcf.gz",
#         indels                              =   postPath+"results.indels.remove_noise_remove_bulkNM0_keepPass.vcf.gz",
#     benchmark:
#         "benchmarks/post/{sample}.txt"
#     log:
#         "logs/post/{sample}.log"
#     params:
#         fasta=config["fasta"],
#         noise_wgns=config["noise_wgns"],
#         snp_wgns=config["snp_wgns"],
#         gnomad=config["gnomad"],
#         dbsnp=config["dbsnp"],
#     threads: 10
#     resources:
#         mem_mb                  =   5000,
#         runtime                 =   120,
#     conda:
#         workflow.source_path("envs/nanoseq_snakemake.yml")
#     group:
#         "post"
#     shell:
#         """
#         # here, filter out variants that might fall within noise and snp masks, just in case you accidentally ran rens
#         # also overlap mutations with gnomAD and retain the gnomAD AFs
#         # keep PASS and non-PASS mutations so that we can separately model the rate of overlap
#         """

# # add a workflow here to call germline mutations from control BAM, check if they occur within the covered regions (at a4s2), and determine if they are present within the NanoSeq BAMs