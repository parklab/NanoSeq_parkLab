# function to check if control BAMs exist

CHROMS = [i+1 for i in range(24)]

# REWRITE THIS WROFKLOW TO SPECIFY THOUSANDS OF CHUNKS AND INSTEAD ASSESS 10KB INTERVALS AT A TIME. 

# I need to give this a run manually in a separate directory outside of Snakefile s thatI know what the workflow steps are like

# for this analysis, we will split BAMs by chromosome to conduct variant calling
# INTERVALS, = glob_wildcards("intervals/{interval}.intervals.list")

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
rule check_efficiency:
    input:
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=check_ctrl_bam_exists,
    output:
        efficiency_stats="efficiency/{sample}.efficiency_stats",
        # read_bundle_read_counts="efficiency/{sample}.efficiency_stats",
    benchmark:
        "benchmarks/extract_tags/{sample}.txt"
    log:
        "logs/check_efficiency/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 10
    shell:
        """
        ln -s ./bin/efficiency_nanoseq.pl
        ln -s ./bin/efficiency_nanoseq.R
        perl efficiency_nanoseq.pl \
        -t {threads} \
        -duplex {input.duplex_bam} \
        -dedup {input.ctrl_bam} \
        -o {output.efficiency_stats} \
        -r {params.fasta}
        """
        

# specify chromosomes that must finish first...go from 1 to 24
rule coverage_histogram_controlBam:
    input:
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=check_ctrl_bam_exists,
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
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=check_ctrl_bam_exists,
    output:
        coverage="{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
    benchmark:
        "benchmarks/partition_coverage/{sample}.txt"
    log:
        "logs/partition_coverage/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        n_partitions=60,
        jobs=n_jobs,
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 10
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


rule dsa_beds:
    input:
        indir=rules.coverage_histogram_controlBam.output.runNanoSeqDir,
        duplex_bam=rules.mark_read_bundles.output.outbam,
        ctrl_bam=check_ctrl_bam_exists,
        coverage=rules.partition_coverage.output.coverage, # "{sample}.runNanoSeq/tmpNanoSeq/part/args.json",
    output:
        coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.dsa.bed.gz",job=dsa_threads,allow_missing=True),
        jobDone=expand("{sample}.runNanoSeq/tmpNanoSeq/dsa/{job}.done",job=dsa_threads,allow_missing=True),
    benchmark:
        "benchmarks/dsa_beds/{sample}.txt"
    log:
        "logs/dsa_beds/{sample}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        snp="test/SNP.sorted.bed.gz",
        noise="test/NOISE.sorted.bed.gz",
        jobs=n_jobs,
    resources:
        mem_mb=20000,
        runtime=60*48,
    threads: 1
    shell:
        """
        PATH=$PATH:$PWD/bin/

        cd {input.indir}

        runNanoSeq.py \
        -t {params.jobs} \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -R {params.fasta} \
        dsa \
        -C ../{params.snp} \
        -D ../{params.noise} \
        -d 2 \
        -q 30 \

        # dsa \
        # -A ../{input.ctrl_bam} \
        # -B ../{input.duplex_bam} \
        # -C ../{params.snp} \
        # -D ../{params.noise} \
        # -R {params.fasta} \
        # -d 2 \
        # -Q 30 \
        # -M 0 \
        # -r "1" \
        # -b 0 \
        # -e 249250620 > ./tmpNanoSeq/dsa/1.dsa.bed

        """

# # NOT CLEAR WHAT THE OUTPUTS ARE HERE...
# rule variant_tables:
#     input:
#         indir=rules.coverage_histogram_controlBam.runNanoSeqDir,
#         duplex_bam=rules.mark_read_bundles.output.outbam,
#         ctrl_bam=check_ctrl_bam_exists,
#     output:
#         coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.var.bed",job=dsa_threads,allow_missing=True),
#     benchmark:
#         "benchmarks/extract_tags/{sample}.txt"
#     log:
#         "logs/extract_tags/{sample}.log"
#     params: # a preset for the ultrashear or covaris libraries
#         fasta=config["fasta"],
#         snp="test/SNP.sorted.bed.gz",
#         noise="test/NOISE.sorted.bed.gz",
#         jobs=n_jobs,
#     resources:
#         mem_mb=10000,
#         runtime=240,
#     threads: 1
#     shell:
#         """
        
#         cd $indir

#         runNanoSeq.py \
#         -t {params.jobs} \
#         -A ../{input.ctrl_bam} \
#         -B ../{input.duplex_bam} \
#         -R {params.fasta} \
#         var \
#         -a 50 \
#         -b 5 \
#         -c 0 \
#         -f 0.9 \
#         -i 1 \
#         -m 8 \
#         -n 3 \
#         -p 0 \
#         -q 60 \
#         -r 144 \
#         -v 0.01 \
#         -x 8 \
#         -z 12
#         """


# rule indel_vcf:
#     input:
#         indir=rules.coverage_histogram_controlBam.runNanoSeqDir,
#         duplex_bam=rules.mark_read_bundles.output.outbam,
#         ctrl_bam=check_ctrl_bam_exists,
#     output:
#         coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/indel/{job}.var.bed",job=dsa_threads,allow_missing=True),
#     benchmark:
#         "benchmarks/extract_tags/{sample}.txt"
#     log:
#         "logs/extract_tags/{sample}.log"
#     params: # a preset for the ultrashear or covaris libraries
#         fasta=config["fasta"],
#         snp="test/SNP.sorted.bed.gz",
#         noise="test/NOISE.sorted.bed.gz",
#         jobs=n_jobs,
#     resources:
#         mem_mb=10000,
#         runtime=240,
#     threads: 1
#     shell:
#         """
        
#         cd $indir

#         runNanoSeq.py \
#         -t {params.jobs} \
#         -A ../{input.ctrl_bam} \
#         -B ../{input.duplex_bam} \
#         -R {params.fasta} \
#         indel \
#         -s sample \
#         --rb 2 \
#         --t3 135 \
#         --t5 10 \
#         --mc 16

#         """

# rule post_process_FINAL:
#     input:
#         indir=rules.coverage_histogram_controlBam.runNanoSeqDir,
#         duplex_bam=rules.mark_read_bundles.output.outbam,
#         ctrl_bam=check_ctrl_bam_exists,
#     output:
#         coverage=expand("{sample}.runNanoSeq/tmpNanoSeq/var/{job}.var.bed",job=dsa_threads,allow_missing=True),
#     benchmark:
#         "benchmarks/extract_tags/{sample}.txt"
#     log:
#         "logs/extract_tags/{sample}.log"
#     params: # a preset for the ultrashear or covaris libraries
#         fasta=config["fasta"],
#         snp="test/SNP.sorted.bed.gz",
#         noise="test/NOISE.sorted.bed.gz",
#         jobs=n_jobs,
#     resources:
#         mem_mb=10000,
#         runtime=240,
#     threads: 1
#     shell:
#         """
        
#         cd $indir

#         runNanoSeq.py \
#         -t {params.jobs} \
#         -A ../{input.ctrl_bam} \
#         -B ../{input.duplex_bam} \
#         -R {params.fasta} \
#         var \
#         -a 50 \
#         -b 5 \
#         -c 0 \
#         -f 0.9 \
#         -i 1 \
#         -m 8 \
#         -n 3 \
#         -p 0 \
#         -q 60 \
#         -r 144 \
#         -v 0.01 \
#         -x 8 \
#         -z 12
#         """




# rule extract_tags:
#     input:
#         expand("input_duplex_fq/{sample}.{pe}.fastq.gz",pe=PAIRED_ENDS,allow_missing=True)
#         # "input_duplex_fq/{sample}.R1.fastq.gz",
#         # "input_duplex_fq/{sample}.R2.fastq.gz",
#     output:
#         expand("extracted_tags_duplex_fq/{sample}.{pe}.fastq.gz",pe=PAIRED_ENDS,allow_missing=True)
#         # "extracted_tags_duplex_fq/{sample}.R1.fastq.gz",
#         # "extracted_tags_duplex_fq/{sample}.R2.fastq.gz",
#     benchmark:
#         "benchmarks/extract_tags/{sample}.txt"
#     log:
#         "logs/extract_tags/{sample}.log"
#     params: # a preset for the ultrashear or covaris libraries
#         bases_trim=3,
#         bases_skip=2,
#         read_len=150,
#     resources:
#         mem_mb=10000,
#         runtime=240,
#     threads: 1
#     shell:
#         """
#         python ./python/extract_tags.py \
#         -a {input[0]} \
#         -b {input[1]} \
#         -c {output[0]} \
#         -d {output[1]} \
#         -m {params.bases_trim} \
#         -s {params.bases_skip} \
#         -l {params.read_len}
#         """

# rule align:
#     input:
#         expand(rules.extract_tags.output,allow_missing=True),
#     output:
#         outsam="aligned_duplex/{sample}.sam"
#     benchmark:
#         "benchmarks/align/{sample}.txt"
#     log:
#         "logs/align/{sample}.log"
#     params:
#         fasta=config["fasta"],
#     resources:
#         mem_mb=30000,
#         runtime=120,
#     threads: 36
#     shell:
#         """
#         sentieon \
#         bwa mem \
#         -t {threads} \
#         -K 10000000 \
#         -C \
#         {config[fasta]} \
#         {input} > {output.outsam} || exit 1
#         """

# rule prepare_RcMcOd_tags:
#     input:
#         rules.align.output.outsam
#     output:
#         tmpdir=temp(directory("{sample}_tmp/")),
#         outbam="rcMcOd_duplex/{sample}.od.bam"
#     benchmark:
#         "benchmarks/prepare_RcMcOd_tags/{sample}.txt"
#     log:
#         "logs/prepare_RcMcOd_tags/{sample}.log"
#     params:
#         fasta=config["fasta"],
#     resources:
#         mem_mb=30000,
#         runtime=480,
#     threads: 12
#     shell:
#         """
#         ulimit -c unlimited
#         ulimit -n 8192
#         tmpdir={wildcards.sample}_tmp
#         mkdir -p $tmpdir
#         echo -e "{wildcards.sample}" >> {log}
#         bamsormadup inputformat=sam rcsupport=1 threads=1 tmpfile=$tmpdir/{wildcards.sample} < {input} > {output.outbam}
#         """
#         # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 

# rule mark_optical_duplicates:
#     input:
#         rules.prepare_RcMcOd_tags.output.outbam
#     output:
#         outbam="rcMcOd_duplex/{sample}.marked_od.bam",
#         metrics="rcMcOd_duplex/{sample}.mark_od_metrics.txt",
#     benchmark:
#         "benchmarks/mark_optical_duplicates/{sample}.txt"
#     log:
#         "logs/mark_optical_duplicates/{sample}.log"
#     params:
#         fasta=config["fasta"],
#         outduplicates="rcMcOd_duplex/{sample}.duplicates.bam",
#         tmpdir=temp(directory("{sample}_tmp2/")),
#     resources:
#         mem_mb=30000,
#         runtime=480,
#     threads: 12
#     shell:
#         """
#         ulimit -c unlimited
#         ulimit -n 8192
#         tmpdir={params.tmpdir}
#         mkdir -p $tmpdir
#         echo -e "{wildcards.sample}" >> {log}
#         bammarkduplicatesopt \
#         I={input} \
#         O={output.outbam} \
#         M={output.metrics} \
#         D={params.outduplicates} \
#         tmpfile=$tmpdir/{wildcards.sample} \
#         index=1 \
#         optminpixeldif=2500 \
#         inputformat=bam \
#         outputformat=bam \
#         inputthreads={threads} \
#         outputthreads={threads}
#         """
#         # bamsormadup inputformat=sam rcsupport=1 threads={threads} < {input} > {output.outbam} 


# rule mark_read_bundles:
#     input:
#         rules.mark_optical_duplicates.output.outbam,
#     output:
#         outbam="readBundle_duplex/{sample}.filtered.bam",
#         outbamInd="readBundle_duplex/{sample}.filtered.bam.bai",
#     benchmark:
#         "benchmarks/mark_read_bundles/{sample}.txt"
#     log:
#         "logs/mark_read_bundles/{sample}.log"
#     params:
#         fasta=config["fasta"],
#     resources:
#         mem_mb=30000,
#         runtime=480,
#     threads: 36
#     shell:
#         """
#         bin/bamaddreadbundles -I {input} -O {output.outbam} || exit 1
#         samtools index {output.outbam} || exit 1
#         """

