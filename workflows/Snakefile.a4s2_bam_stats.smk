rule sort_bam_by_rb:
    input:
        rb_bam="readBundle_duplex/{sample}.filtered.bam"
    output:
        byRb="a4s2_bundles/{sample}.sort_by_byRb.bam",
    benchmark:
        "benchmarks/sort_bam_by_rb/{sample}.txt"
    log:
        "logs/sort_bam_by_rb/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        #
        samtools view -bh -d RB -F 780 {input.rb_bam} | \
        samtools sort -t RB - | \
        samtools view -CT {params.fasta} -o {output.byRb} - || exit 1
        """

rule get_covstats_rb_bam:
    input:
        rb_bam="readBundle_duplex/{sample}.filtered.bam",
    output:
        covStatsAll="a4s2_bundles/{sample}.all.covStats.tsv",
    benchmark:
        "benchmarks/get_covstats_rb_bam/{sample}.txt"
    log:
        "logs/get_covstats_rb_bam/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        samtools stats {input.rb_bam} > {output.covStatsAll} || exit 1 # stats for the original BAM; consider duplicates to assess raw coverage
        """

rule isolate_a4s2_bam:
    input:
        rb_bam="readBundle_duplex/{sample}.filtered.bam",
        byRb=rules.sort_bam_by_rb.output.byRb,
    output:
        outbam="a4s2_bundles/{sample}.a4s2.bam",
        outbamInd="a4s2_bundles/{sample}.a4s2.bam.bai",
        # outbam_downSampled="a4s2_bundles/{sample}.a4s2.non_duplicate_read_per_bundle.bam",
        # outbam_downSampled_ind="a4s2_bundles/{sample}.a4s2.non_duplicate_read_per_bundle.bam.bai",        
    benchmark:
        "benchmarks/isolate_a4s2_bam/{sample}.txt"
    log:
        "logs/isolate_a4s2_bam/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        # python script that tracks whether an RB tag is present and if so, adds an "a4s2" tag to the read, then outputs a BAM; the python script indexes it
        python workflows/filter_for_a4s2.py --bam {input.byRb} --out {output.outbam} --status {log} || exit 1
        """
        # samtools view -bhF 1024 -o {output.outbam_downSampled} {output.outbam} || exit 1
        # samtools index {output.outbam_downSampled} || exit 1

rule get_covstats_a4s2_bam:
    input:
        inbam=rules.isolate_a4s2_bam.output.outbam,
    output:
        covStatsA2S4="a4s2_bundles/{sample}.a4s2.covStats.tsv",
    benchmark:
        "benchmarks/get_covstats_a4s2_bam/{sample}.txt"
    log:
        "logs/get_covstats_a4s2_bam/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        samtools stats -d {input.inbam} > {output.covStatsA2S4} || exit 1 # stats for the a4s2 BAM; ignore duplicates to assess cosensus coverage
        """

rule summarize_covStats_a4s2_originalRb_Bams:
    input:
        covStatsA2S4=rules.get_covstats_a4s2_bam.output.covStatsA2S4,
        covStatsAll=rules.get_covstats_rb_bam.output.covStatsAll,
    output:
        estimated_depth="a4s2_bundles/{sample}.a4s2.avg_depth.tsv",
    benchmark:
        "benchmarks/summarize_covStats_a4s2_originalRb_Bams/{sample}.txt"
    log:
        "logs/summarize_covStats_a4s2_originalRb_Bams/{sample}.log"
    params:
        fasta=config["fasta"],
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        #
        echo -e "summarize the stats files to get the average depth across all genomic positions covered by the BAM file, and output a table with the number of bases, number of genomic positions covered, and average depth for both the a4s2 BAM and the original BAM"
        cat <(echo -e "Number of bases sequenced\tNumber of genomic positions covered\tAverage depth") \
        <(grep ^COV {input.covStatsA2S4} | cut -f 3- | awk 'BEGIN {{FS=OFS="\t"}} {{sum_depths += $1 * $2; sum_weights += $2}} END {{print sum_depths,sum_weights,sum_depths/sum_weights}}') \
        <(grep ^COV {input.covStatsAll} | cut -f 3- | awk 'BEGIN {{FS=OFS="\t"}} {{sum_depths += $1 * $2; sum_weights += $2}} END {{print sum_depths,sum_weights,sum_depths/sum_weights}}') | \
        paste <(echo -e "BAM\na2s4 duplex consensus\nall reads") - > {output.estimated_depth}
        """
