# here, we associate each sample with its control, and then run the entire NanoSeq pipeline for each sample-control pair.

# need to evaluate this...

# def get_donor_expected_sample_bam(wildcards):
#     # import DONOR/{wildcards.donor}.txt and get the expected control BAM file
#     if not os.path.exists(f"DONOR/{wildcards.donor}.txt"):
#         raise ValueError(f"Expected file DONOR/{wildcards.donor}.txt not found")
#     inTable = pd.read_csv(f"DONOR/{wildcards.donor}.txt",sep="\t",header=0)
#     sample_bam = inTable.loc[inTable["Type"] == "ExpectedTestBAM", "Value"].values[0]
#     return sample_bam

# def get_donor_expected_control_bam(wildcards):
#     # import DONOR/{wildcards.donor}.txt and get the expected control BAM file
#     if not os.path.exists(f"DONOR/{wildcards.donor}.txt"):
#         raise ValueError(f"Expected file DONOR/{wildcards.donor}.txt not found")
#     inTable = pd.read_csv(f"DONOR/{wildcards.donor}.txt",sep="\t",header=0)
#     control_bam = inTable.loc[inTable["Type"] == "ExpectedControlBAM", "Value"].values[0]
#     return control_bam

# ---- entry-point compatibility ---------------------------------------------
# These modules are donor-keyed, but they are included from entry points that
# key on either {donor} (Snakefile.table_input.mito_nanoseq.smk, the by_donor
# workflow) or {sample} (Snakefile.testSample_bam.ctrlSample_bam.table_input).
# The sample-keyed parents do not define donor_expected_bams, so provide it
# here rather than requiring every parent to carry a copy.
if "donor_expected_bams" not in globals():
    def donor_expected_bams(donor):
        """Per-donor duplex/control BAM paths, resolved from input.tsv."""
        _t = pd.read_csv("input.tsv", sep="\t", header=0)
        _r = _t.loc[_t["Donor"] == donor]
        if not len(_r):
            raise ValueError(f"donor {donor!r} not found in input.tsv")
        return {
            "test": f"readBundle_duplex/{_r['TestBamID'].values[0]}.filtered.bam",
            "control": f"control_duplex/{_r['ControlBamID'].values[0]}.diluted.ctrl.bam",
        }

if "DONOR" not in globals():
    DONOR = set(pd.read_csv("input.tsv", sep="\t", header=0)["Donor"].values)

if "SUMMARY_RESULT_FILES" not in globals():
    SUMMARY_RESULT_FILES = ["burdens", "callvsqpos", "coverage",
                            "pyrvsmask", "readbundles"]


def get_donor_expected_sample_bam(wildcards):
    return donor_expected_bams(wildcards.donor)["test"]
    # inTable = pd.read_csv("input.tsv", sep="\t", header=0)
    # testbam = inTable.loc[inTable["Donor"] == wildcards.donor, "TestBamID"].values[0]
    # return f"readBundle_duplex/{testbam}.filtered.bam"

def get_donor_expected_control_bam(wildcards):
    return donor_expected_bams(wildcards.donor)["control"]
    # inTable = pd.read_csv("input.tsv", sep="\t", header=0)
    # control = inTable.loc[inTable["Donor"] == wildcards.donor, "ControlBamID"].values[0]
    # return f"control_duplex/{control}.filtered.bam"


def get_control_bam(wildcards):
    inTable = pd.read_csv("input.tsv",sep="\t",header=0)
    control = inTable.loc[inTable["Donor"]==wildcards.donor,"ControlBamID"].values[0]
    control_bam = f"control_duplex/{control}.diluted.ctrl.bam"
    return control_bam

def get_sample_bam(wildcards):
    inTable = pd.read_csv("input.tsv",sep="\t",header=0)
    sample = inTable.loc[inTable["Donor"]==wildcards.donor,"TestBamID"].values[0]
    sample_bam = f"readBundle_duplex/{sample}.filtered.bam"
    return sample_bam


# Retry-aware runtime helpers. Each rule's `resources.runtime` should be a
# lambda over (wildcards, attempt) so that Snakemake's `--retries N` actually
# requests a larger wall-clock from SLURM on subsequent attempts.
def _get_dynamic_runtime(attempt, basetime, increment, label="rule"):
    total = int(basetime + increment * (attempt - 1))
    print(f"[runtime] attempt {attempt} for {label}: requesting {total} min (base={basetime}, +{increment}/attempt)")
    return total

def get_dynamic_runtime_efficiency(attempt):
    return _get_dynamic_runtime(attempt, basetime=int(60*0.5), increment=60, label="check_mito_efficiency")

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


# Name of the mitochondrial contig in the reference. UCSC-style references
# (hg38_no_alt.fa) call it "chrM"; Ensembl/GRCh37d5 calls it "MT".
MITO_CONTIG = config.get("mito_contig", "chrM")

# Optional BED of NUMT segments in mito coordinates. When unset, the coverage
# plots simply omit that ring and the low-MAPQ track stands in for the same
# concern. Derive one (e.g. by aligning the mito contig against the nuclear
# contigs) and point this key at it to switch the ring on.
NUMT_BED = config.get("numt_bed", None)

# Coverage tracks are keyed on TestBamID rather than Donor: donors differing
# only in their matched control share a test BAM, so keying on Donor would
# compute and plot the same mitochondrial tracks two or more times.
def _mito_sample_maps():
    tbl = pd.read_csv("input.tsv", sep="\t", header=0)
    donor_of, frag_of = {}, {}
    for _, row in tbl.iterrows():
        sample = row["TestBamID"]
        donor_of.setdefault(sample, row["Donor"])
        if "Fragmentation" in tbl.columns:
            frag_of.setdefault(sample, str(row["Fragmentation"]))
        else:
            frag_of.setdefault(sample, "NA")
    return donor_of, frag_of

MITO_DONOR_OF_SAMPLE, MITO_FRAG_OF_SAMPLE = _mito_sample_maps()
MITO_SAMPLES = sorted(MITO_DONOR_OF_SAMPLE)


def get_mito_a4s2_bam(wildcards):
    return f"mito_a4s2_bundles/{MITO_DONOR_OF_SAMPLE[wildcards.sample]}.a4s2.bam"


# verify consistency of input files
rule verify_donor_bam_consistency:
    input:
        sheet = "DONOR/{donor}.txt",
    output:
        "DONOR/{donor}.verified",
    resources:
        mem_mb=2000, runtime=5,
    run:
        sheet = pd.read_csv(input.sheet, sep="\t", header=0)
        sheet_test = sheet.loc[sheet["Type"] == "ExpectedTestBAM", "Value"].values[0]
        sheet_ctrl = sheet.loc[sheet["Type"] == "ExpectedControlBAM", "Value"].values[0]
        fn_test = get_donor_expected_sample_bam(wildcards)
        fn_ctrl = get_donor_expected_control_bam(wildcards)
        mismatches = []
        if sheet_test != fn_test:
            mismatches.append(f"test: sheet={sheet_test!r} fn={fn_test!r}")
        if sheet_ctrl != fn_ctrl:
            mismatches.append(f"control: sheet={sheet_ctrl!r} fn={fn_ctrl!r}")
        if mismatches:
            raise ValueError(f"{wildcards.donor}: " + "; ".join(mismatches))
        with open(output[0], "w") as f:
            f.write(f"{wildcards.donor}\t{fn_test}\t{fn_ctrl}\tOK\n")


rule check_mito_efficiency:
    # FOR THIS JOB, THE SANGER'S EFFICIENCY SCRIPT WILL ASSESS JUST THE FIRST CONTIG FROM THE FASTA
    input:
        duplex_bam=get_donor_expected_sample_bam,
        ctrl_bam=get_donor_expected_control_bam,
    output:
        tempFasta=temp("mito_efficiency/{donor}.temp.fa"),
        tempFastaIndex=temp("mito_efficiency/{donor}.temp.fa.fai"),
        read_bundles="mito_efficiency/{donor}.RBs",
        read_bundles_gc_inserts="mito_efficiency/{donor}.RBs.GC_inserts.tsv",
        read_bundles_pdf="mito_efficiency/{donor}.RBs.pdf",
        efficiency_stats="mito_efficiency/{donor}.tsv",
    benchmark:
        "benchmarks/check_mito_efficiency/{donor}.txt"
    log:
        "logs/check_mito_efficiency/{donor}.log"
    params: # a preset for the ultrashear or covaris libraries
        fasta=config["fasta"],
        mito_contig_name=MITO_CONTIG,
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: get_dynamic_runtime_efficiency(attempt),
    threads: 20
    conda:
        workflow.source_path("../envs/nanoseq_snakemake.yml")
    group:
        "check_mito_efficiency"
    shell:
        """
        PATH=$PATH:$PWD/bin/:$PWD/perl/
        # efficiency_nanoseq.pl restricts its read-bundle metrics to the FIRST
        # contig of the supplied reference index, so hand it a reference that
        # contains only the mitochondrial contig.
        samtools faidx {params.fasta} {params.mito_contig_name} > {output.tempFasta} || exit 1
        samtools faidx {output.tempFasta} || exit 1
        perl $PWD/perl/efficiency_nanoseq.pl \
        -t {threads} \
        -duplex {input.duplex_bam} \
        -dedup {input.ctrl_bam} \
        -out mito_efficiency/{wildcards.donor} \
        -ref {output.tempFasta}
        """
        # 2> {logs} 1> {logs}

# IN THE FOLLOWING RULES, WE WILL ISOLATE JUST THE A4S2 DUPLEXES FROM MITOCHONDRIAL GENOME
rule sort_mito_bam_by_rb:
    input:
        # rb_bam="readBundle_duplex/{sample}.filtered.bam"
        rb_bam=get_donor_expected_sample_bam,
    output:
        byRb="mito_a4s2_bundles/{donor}.sort_by_byRb.bam",
    benchmark:
        "benchmarks/sort_mito_bam_by_rb/{donor}.txt"
    log:
        "logs/sort_mito_bam_by_rb/{donor}.log"
    params:
        fasta=config["fasta"],
        mito_contig_name=MITO_CONTIG,
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        # Sorted by RB tag, not by coordinate, so keep it as BAM -- a CRAM here
        # would need a reference lookup on every read for no benefit.
        samtools view -bh -d RB -F 780 {input.rb_bam} {params.mito_contig_name} | \
        samtools sort -t RB -@ {threads} -o {output.byRb} - || exit 1
        """

rule get_mito_covstats_rb_bam:
    input:
        # rb_bam="readBundle_duplex/{sample}.filtered.bam",
        rb_bam=get_donor_expected_sample_bam,
    output:
        covStatsAll="mito_a4s2_bundles/{donor}.all.covStats.tsv",
    benchmark:
        "benchmarks/get_mito_covstats_rb_bam/{donor}.txt"
    log:
        "logs/get_mito_covstats_rb_bam/{donor}.log"
    params:
        fasta=config["fasta"],
        mito_contig_name=MITO_CONTIG,
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        # -c raises the coverage histogram ceiling from the 1000x default: mitochondrial
        # depth runs orders of magnitude higher, and every position landing in the
        # overflow bin pins the COV-derived mean depth at exactly 1000.
        samtools stats -c 1,1000000,1 {input.rb_bam} {params.mito_contig_name} > {output.covStatsAll} || exit 1 # stats for the original BAM; consider duplicates to assess raw coverage
        """

rule isolate_mito_a4s2_bam:
    input:
        # rb_bam="readBundle_duplex/{sample}.filtered.bam",
        rb_bam=get_donor_expected_sample_bam,
        byRb=rules.sort_mito_bam_by_rb.output.byRb,
    output:
        outbam="mito_a4s2_bundles/{donor}.a4s2.bam",
        outbamInd="mito_a4s2_bundles/{donor}.a4s2.bam.bai",
        status="mito_a4s2_bundles/{donor}.a4s2.progress.txt",
        # outbam_downSampled="mito_a4s2_bundles/{sample}.a4s2.non_duplicate_read_per_bundle.bam",
        # outbam_downSampled_ind="mito_a4s2_bundles/{sample}.a4s2.non_duplicate_read_per_bundle.bam.bai",
    benchmark:
        "benchmarks/isolate_mito_a4s2_bam/{donor}.txt"
    log:
        "logs/isolate_mito_a4s2_bam/{donor}.log"
    params:
        fasta=config["fasta"],
        mito_contig_name=MITO_CONTIG,
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        # python script that tracks whether an RB tag is present and if so, adds an "a4s2" tag to the read, then outputs a BAM; the python script indexes it
        # NB: --status is the progress counter, {log} captures stdout/stderr including tracebacks
        python workflows/filter_for_a4s2.py --bam {input.byRb} --out {output.outbam} --status {output.status} > {log} 2>&1 || exit 1
        """
        # samtools view -bhF 1024 -o {output.outbam_downSampled} {output.outbam} || exit 1
        # samtools index {output.outbam_downSampled} || exit 1

rule get_mito_covstats_a4s2_bam:
    input:
        inbam=rules.isolate_mito_a4s2_bam.output.outbam,
    output:
        covStatsA2S4="mito_a4s2_bundles/{donor}.a4s2.covStats.tsv",
    benchmark:
        "benchmarks/get_mito_covstats_a4s2_bam/{donor}.txt"
    log:
        "logs/get_mito_covstats_a4s2_bam/{donor}.log"
    params:
        fasta=config["fasta"],
        mito_contig_name=MITO_CONTIG,
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        samtools stats -d -c 1,1000000,1 {input.inbam} {params.mito_contig_name} > {output.covStatsA2S4} || exit 1 # stats for the a4s2 BAM; ignore duplicates to assess cosensus coverage
        """

rule summarize_mito_covStats_a4s2_originalRb_Bams:
    input:
        covStatsA2S4=rules.get_mito_covstats_a4s2_bam.output.covStatsA2S4,
        covStatsAll=rules.get_mito_covstats_rb_bam.output.covStatsAll,
        a4s2_bam=rules.isolate_mito_a4s2_bam.output.outbam,
        rb_bam=get_donor_expected_sample_bam,
    output:
        estimated_depth="mito_a4s2_bundles/{donor}.a4s2.avg_depth.tsv",
    benchmark:
        "benchmarks/summarize_mito_covStats_a4s2_originalRb_Bams/{donor}.txt"
    log:
        "logs/summarize_mito_covStats_a4s2_originalRb_Bams/{donor}.log"
    params:
        fasta=config["fasta"],
        mito_contig_name=MITO_CONTIG,
    resources:
        mem_mb=30000,
        runtime=480,
    threads:
        36
    shell:
        """
        # Depth is taken from `samtools coverage`, which reports mean depth over the
        # whole contig and has no binning ceiling. The COV-derived columns from
        # `samtools stats` are kept alongside it, but note they average only over
        # positions with non-zero coverage.
        #
        # all reads:  duplicates retained    -> raw sequencing depth over the contig
        # a4s2:       default --ff drops DUP -> one read pair per read bundle,
        #                                       i.e. duplex consensus depth
        summarize() {{
            local label="$1" covstats="$2"
            paste \
              <(echo -e "$label") \
              <(grep ^COV "$covstats" | cut -f 3- | awk 'BEGIN {{FS=OFS="\t"}} {{sum_depths += $1 * $2; sum_weights += $2}} END {{if (sum_weights > 0) print sum_depths,sum_weights,sum_depths/sum_weights; else print 0,0,0}}')
        }}

        {{
          echo -e "BAM\tbases_sequenced\tpositions_covered\tmean_depth_over_covered\tcontig_len\tcovbases\tpct_contig_covered\tmean_depth_over_contig\tmean_baseq\tmean_mapq"
          paste \
            <(summarize "a4s2 duplex consensus" {input.covStatsA2S4}) \
            <(samtools coverage -r {params.mito_contig_name} {input.a4s2_bam} | awk 'BEGIN {{FS=OFS="\t"}} NR>1 {{print $3,$5,$6,$7,$8,$9}}')
          paste \
            <(summarize "all reads" {input.covStatsAll}) \
            <(samtools coverage --ff UNMAP,SECONDARY,QCFAIL -r {params.mito_contig_name} {input.rb_bam} | awk 'BEGIN {{FS=OFS="\t"}} NR>1 {{print $3,$5,$6,$7,$8,$9}}')
        }} > {output.estimated_depth}
        """



# for i in *filtered.bam; do samtools view -q 30 $i chrM | awk '{
#        as=""; xs="";
#        for(i=1; i<=NF; i++) {
#            if($i ~ /^AS:/) as=$i;
#            if($i ~ /^XS:/) xs=$i;
#        }
#        print $1, as, xs;
#    }' | sed "s/ /\t/g" | sed "s/AS:i://g" | sed "s/XS:i://g" | awk '{print $0,$2-$3}' | awk '$4 > 100' | grep -c "^" | paste <(echo -e "$i") -; done

# ---------------------------------------------------------------------------
# Per-base mitochondrial coverage tracks and Circos-style plots.
#
# Two depth notions are carried side by side, because for duplex data they
# differ by more than an order of magnitude: read depth counts every read,
# while read-bundle depth counts distinct source molecules. The latter is what
# duplex sensitivity actually scales with.
# ---------------------------------------------------------------------------

rule mito_per_base_tracks:
    input:
        raw_bam=READBUNDLE_DIR + "{sample}.filtered.bam",
        raw_bai=READBUNDLE_DIR + "{sample}.filtered.bam.bai",
        a4s2_bam=get_mito_a4s2_bam,
    output:
        tracks="mito_tracks/{sample}.per_base.tsv.gz",
    benchmark:
        "benchmarks/mito_per_base_tracks/{sample}.txt"
    log:
        "logs/mito_per_base_tracks/{sample}.log"
    params:
        mito_contig_name=MITO_CONTIG,
        fragmentation=lambda wildcards: MITO_FRAG_OF_SAMPLE[wildcards.sample],
        min_mapq=30,
    resources:
        mem_mb=8000,
        runtime=60,
    threads: 2
    shell:
        """
        python workflows/mito_per_base_tracks.py \
        --raw-bam {input.raw_bam} \
        --a4s2-bam {input.a4s2_bam} \
        --contig {params.mito_contig_name} \
        --sample {wildcards.sample} \
        --fragmentation {params.fragmentation} \
        --min-mapq {params.min_mapq} \
        --out {output.tracks} > {log} 2>&1 || exit 1
        """


rule plot_mito_circos_sample:
    input:
        tracks="mito_tracks/{sample}.per_base.tsv.gz",
        numt=[NUMT_BED] if NUMT_BED else [],
    output:
        pdf="mito_plots/{sample}.circos.pdf",
    benchmark:
        "benchmarks/plot_mito_circos_sample/{sample}.txt"
    log:
        "logs/plot_mito_circos_sample/{sample}.log"
    params:
        bin_size=10,
        numt_arg=("--numt-bed " + NUMT_BED) if NUMT_BED else "",
    resources:
        mem_mb=8000,
        runtime=30,
    threads: 1
    shell:
        """
        Rscript workflows/plot_mito_circos.R \
        --mode sample \
        --tracks {input.tracks} \
        --bin {params.bin_size} {params.numt_arg} \
        --out {output.pdf} > {log} 2>&1 || exit 1
        """


rule plot_mito_circos_combined:
    input:
        tracks=expand("mito_tracks/{sample}.per_base.tsv.gz", sample=MITO_SAMPLES),
        numt=[NUMT_BED] if NUMT_BED else [],
    output:
        pdf="mito_plots/all_samples.circos.pdf",
        summary="mito_tracks/summary.tsv",
    benchmark:
        "benchmarks/plot_mito_circos_combined/all_samples.txt"
    log:
        "logs/plot_mito_circos_combined/all_samples.log"
    params:
        bin_size=10,
        numt_arg=("--numt-bed " + NUMT_BED) if NUMT_BED else "",
    resources:
        mem_mb=16000,
        runtime=60,
    threads: 1
    shell:
        """
        Rscript workflows/plot_mito_circos.R \
        --mode combined \
        --tracks "$(echo {input.tracks} | tr ' ' ',')" \
        --bin {params.bin_size} {params.numt_arg} \
        --summary {output.summary} \
        --out {output.pdf} > {log} 2>&1 || exit 1
        """


# Mitochondrial variant calling (first pass, unmasked). Included last so the
# helpers and MITO_CONTIG defined above are in scope.
include: "Snakefile.table_input.call_mito.from_bam.smk"


# ---------------------------------------------------------------------------
# Empirical noise mask from per-base low-MAPQ fraction.
#
# Built only from the deep, even libraries (the same measured-signature filter
# the calling workflow uses), because a sparse restriction library has no
# coverage over ~70% of chrM and cannot contribute a mappability estimate
# there. Output is bgzip+tabix BED, which is the format dsa -D requires.
# ---------------------------------------------------------------------------
MITO_MASK_THRESHOLD = float(config.get("mito_lowmapq_threshold", 0.20))
MITO_MASK_MIN_DEPTH = int(config.get("mito_lowmapq_min_depth", 20))


def _mito_mask_track_inputs(wildcards=None):
    # Hand over every track and let the script decide eligibility from measured
    # breadth. Filtering here via the donor table looked equivalent but was not:
    # MITO_DONOR_OF_SAMPLE resolves a shared sample to whichever donor appears
    # first in input.tsv, which for four of the seven deep libraries is a
    # cross-donor "-vs-" row that the calling workflow excludes -- silently
    # building the mask from three samples instead of seven.
    return expand("mito_tracks/{sample}.per_base.tsv.gz", sample=MITO_SAMPLES)


rule mito_lowmapq_noise_mask:
    input:
        tracks=_mito_mask_track_inputs,
    output:
        bed="mito_mask/lowmapq_noise.chrM.bed",
        bedgz="mito_mask/lowmapq_noise.chrM.bed.gz",
        tbi="mito_mask/lowmapq_noise.chrM.bed.gz.tbi",
    log:
        "logs/mito_lowmapq_noise_mask/all.log",
    params:
        contig=MITO_CONTIG,
        threshold=MITO_MASK_THRESHOLD,
        min_depth=MITO_MASK_MIN_DEPTH,
        min_breadth=float(config.get("mito_lowmapq_min_breadth", 0.99)),
    resources:
        mem_mb=4000, runtime=15,
    threads: 1
    shell:
        """
        mkdir -p mito_mask
        python workflows/mito_lowmapq_mask.py \
        --tracks {input.tracks} \
        --contig {params.contig} \
        --threshold {params.threshold} \
        --min-depth {params.min_depth} \
        --min-breadth {params.min_breadth} \
        --out {output.bed} > {log} 2>&1 || exit 1
        # dsa -D needs a tabix-indexed BED, not a plain one.
        bgzip -c {output.bed} > {output.bedgz}
        tabix -p bed -f {output.bedgz}
        """


# ---------------------------------------------------------------------------
# NUMT annotation track (reference-derived, independent of any sample).
#
# Deliberately NOT the mask: ~80% of chrM carries some NUMT >=200 bp, so NUMT
# presence is nearly uninformative for masking. Identity is what predicts
# mismapping, and the annotation carries it so the two can be reasoned about
# separately. See docs/mito_mask_methods.md.
# ---------------------------------------------------------------------------
MITO_TILE_WIN  = int(config.get("mito_numt_tile_window", 300))
MITO_TILE_STEP = int(config.get("mito_numt_tile_step", 100))


rule mito_numt_tiles:
    output:
        fa=temp("mito_mask/numt_tiles.fa"),
    params:
        fasta=config["fasta"],
        contig=MITO_CONTIG,
        win=MITO_TILE_WIN,
        step=MITO_TILE_STEP,
    resources:
        mem_mb=2000, runtime=10,
    threads: 1
    run:
        import subprocess
        length = _mito_contig_length()
        os.makedirs("mito_mask", exist_ok=True)
        with open(output.fa, "w") as fh:
            # Overlapping windows so a NUMT boundary cannot fall in a gap.
            for start in range(1, length + 1, params.step):
                end = min(start + params.win - 1, length)
                region = f"{params.contig}:{start}-{end}"
                fh.write(subprocess.run(
                    ["samtools", "faidx", params.fasta, region],
                    capture_output=True, text=True, check=True).stdout)
                if end == length:
                    break


rule mito_numt_align:
    input:
        fa=rules.mito_numt_tiles.output.fa,
    output:
        sam=temp("mito_mask/numt_tiles.sam"),
    log:
        "logs/mito_numt_align/all.log",
    params:
        fasta=config["fasta"],
    resources:
        # bwa holds the whole-genome index in memory (~5.5 GB for GRCh38);
        # this will be OOM-killed on a login node, so it must go to the cluster.
        mem_mb=24000, runtime=60,
    threads: 4
    shell:
        """
        # -a reports every alignment, not just the best: a NUMT is by definition
        # a secondary hit and would otherwise be discarded.
        bwa mem -a -t {threads} {params.fasta} {input.fa} > {output.sam} 2> {log}
        """


rule mito_numt_annotation:
    input:
        sam=rules.mito_numt_align.output.sam,
    output:
        bed="mito_mask/numt_blocks.chrM.bed",
        bedgz="mito_mask/numt_blocks.chrM.bed.gz",
        tbi="mito_mask/numt_blocks.chrM.bed.gz.tbi",
    log:
        "logs/mito_numt_annotation/all.log",
    params:
        contig=MITO_CONTIG,
        min_block=int(config.get("mito_numt_min_block", 200)),
    resources:
        mem_mb=4000, runtime=15,
    threads: 1
    shell:
        """
        python workflows/mito_numt_annotation.py \
        --sam {input.sam} \
        --contig {params.contig} \
        --min-block {params.min_block} \
        --out {output.bed} > {log} 2>&1 || exit 1
        bgzip -c {output.bed} > {output.bedgz}
        tabix -p bed -f {output.bedgz}
        """


# ---------------------------------------------------------------------------
# Homopolymer track: a stratifier for indel calls, not a mask.
# See docs/mito_mask_methods.md.
# ---------------------------------------------------------------------------
rule mito_homopolymer_track:
    output:
        bed="mito_mask/homopolymers.chrM.bed",
        bedgz="mito_mask/homopolymers.chrM.bed.gz",
        tbi="mito_mask/homopolymers.chrM.bed.gz.tbi",
    log:
        "logs/mito_homopolymer_track/all.log",
    params:
        fasta=config["fasta"],
        contig=MITO_CONTIG,
        min_run=int(config.get("mito_homopolymer_min_run", 5)),
        pad=int(config.get("mito_homopolymer_pad", 1)),
    resources:
        mem_mb=2000, runtime=10,
    threads: 1
    shell:
        """
        mkdir -p mito_mask
        python workflows/mito_homopolymer_track.py \
        --fasta {params.fasta} --contig {params.contig} \
        --min-run {params.min_run} --pad {params.pad} \
        --out {output.bed} > {log} 2>&1 || exit 1
        bgzip -c {output.bed} > {output.bedgz}
        tabix -p bed -f {output.bedgz}
        """
