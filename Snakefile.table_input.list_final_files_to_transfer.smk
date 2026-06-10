# rule cram_convert_original_alignment_sample:
#     input:
#         rules.align_samples.output.tempSam
#     output:
#         md5=FINAL_SAMPLE_CRAM+"{sample}.originalSam.md5",
#         outcram=FINAL_SAMPLE_CRAM+"{sample}.cram",
#         outcram_index=FINAL_SAMPLE_CRAM+"{sample}.cram.crai",
#         fasta_used=FINAL_SAMPLE_CRAM+"{sample}.cram_fasta.txt",
#         md5cram=FINAL_SAMPLE_CRAM+"{sample}.cram.md5",
#     log:
#         "logs/cram_convert_original_alignment_sample/{sample}.txt",
#     group:
#         "cram_convert_original_alignment_sample"
#     benchmark:
#         "benchmarks/cram_convert_original_alignment_sample/{sample}.txt"
#     log:
#         "logs/cram_convert_original_alignment_sample/{sample}.log"
#     params:
#         fasta=config["fasta"],
#     resources:
#         mem_mb=30000,
#         runtime=480,
#     threads:
#         36
#     shell:
#         """
#         md5sum {input} > {output.md5} || exit 1
#         samtools sort -@ {threads} {input} | \
#         samtools view -@ {threads} -C -T {params.fasta} \
#         -o {output.outcram} || exit 1
#         samtools index -@ {threads} {output.outcram} || exit 1
#         echo -e "{params.fasta}" > {output.fasta_used}
#         md5sum {output.outcram} > {output.md5cram} || exit 1
#         """

# rule cram_convert_original_alignment_control:
#     input:
#         rules.align_controls.output.tempSam
#     output:
#         md5=FINAL_CONTROL_CRAM+"{control}.originalSam.md5",
#         outcram=FINAL_CONTROL_CRAM+"{control}.cram",
#         outcram_index=FINAL_CONTROL_CRAM+"{control}.cram.crai",
#         fasta_used=FINAL_CONTROL_CRAM+"{control}.cram_fasta.txt",
#         md5cram=FINAL_CONTROL_CRAM+"{control}.cram.md5",
#     log:
#         "logs/cram_convert_original_alignment_control/{control}.txt",
#     group:
#         "cram_convert_original_alignment_control"
#     benchmark:
#         "benchmarks/cram_convert_original_alignment_control/{control}.txt"
#     log:
#         "logs/cram_convert_original_alignment_control/{control}.log"
#     params:
#         fasta=config["fasta"],
#     resources:
#         mem_mb=30000,
#         runtime=480,
#     threads:
#         36
#     shell:
#         """
#         md5sum {input} > {output.md5} || exit 1
#         samtools sort -@ {threads} {input} | \
#         samtools view -@ {threads} -C -T {params.fasta} \
#         -o {output.outcram} || exit 1
#         samtools index -@ {threads} {output.outcram} || exit 1
#         echo -e "{params.fasta}" > {output.fasta_used}
#         md5sum {output.outcram} > {output.md5cram} || exit 1
#         """


rule list_final_files_to_transfer:
    input:
        "readBundle_duplex/{sample}.filtered.bam",
        "readBundle_duplex/{sample}.filtered.bam.bai",
        "control_duplex/{control}.filtered.bam",
        "control_duplex/{control}.filtered.bam.bai",
        "control_duplex/{control}.diluted.ctrl.bam",
        "control_duplex/{control}.diluted.ctrl.bam.bai",
        postPath,
    output:
        "{sample}.final_files_to_transfer.txt"
    log:
        "logs/list_final_files_to_transfer.txt"
    group:
        "list_final_files_to_transfer"
    benchmark:
        "benchmarks/list_final_files_to_transfer/{sample}.txt"
    log:
        "logs/list_final_files_to_transfer/{sample}.log"
    resources:
        mem_mb=30000,
        runtime=480,
    run:
        with open(output[0], "w") as f:
            for file in input:
                f.write(f"{file}\n")
