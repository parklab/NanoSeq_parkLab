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