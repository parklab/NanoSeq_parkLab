# Mitochondrial variant calling over the mito contig.
#
# Mirrors the rule chain of Snakefile.table_input.run_nanoseq.from_bam.smk but
# keyed on {donor} (as analyze_mito.from_bam.smk already is) and restricted to
# the mitochondrial contig. Included from analyze_mito.from_bam.smk, so the
# helpers defined there -- MITO_CONTIG, get_donor_expected_sample_bam,
# get_donor_expected_control_bam, list_donor_info, get_dynamic_runtime_* --
# are in scope here and are not redefined.
#
# Four deliberate departures from the genome-wide workflow:
#
#  1. cov / run_diagnostic_coverage / partition_coverage are BYPASSED. The
#     genome-wide chain runs `runNanoSeq.py cov` with --exclude "MT,..." over
#     CHROMS 1..24, so the mito contig never gets a coverage file at all --
#     one of three independent reasons mito has never been called here (the
#     others: chrM is absent from the 322 genome-wide interval files, and the
#     Abascal NOISE/SNP masks contain zero chrM intervals). The diagnostic
#     rule fits a logistic to nuclear coverage to estimate a bulk-read floor;
#     extrapolating that to 16,000-23,000x chrM bulk depth is meaningless, so
#     -z is set explicitly from config instead. `runNanoSeq.py post` only
#     validates dsa/var/indel nfiles and .done markers, never cov/ or part/,
#     so nothing downstream misses them.
#
#  2. MASKS are config-driven. -D takes config mito_noise_bed (the empirical
#     low-MAPQ mask); leaving it unset gives a genuinely unmasked run, because
#     src/dsa.cc defaults beds[0..1] to "\0" and Bed::Load returns immediately
#     on an empty filename -- no empty-BED-plus-tabix placeholder needed.
#     -C (germline/SNP) is still unset: no mito germline mask exists yet, so
#     germline heteroplasmy is handled only by -v against the matched bulk.
#     A population filter (e.g. gnomAD v3.1 mtDNA, which is called against
#     rCRS and so needs no liftover here) can be applied post hoc or wired to
#     -C later.
#
#  3. Output root is {donor}.runNanoSeqMito/, NOT {donor}.runNanoSeq/, which
#     already holds genome-wide results and would otherwise be clobbered.
#     Every rule is mito_-prefixed for the same reason.
#
#  4. Depth. Bulk chrM runs ~16,000-23,000x and duplex in the hundreds-to-
#     thousands, versus ~30x nuclear, so dsa/varCall memory is raised well
#     above the genome-wide settings. Run one interval first (see
#     mito_call_single_interval_test in the README notes) before fanning out.

import os

# ---- tunables --------------------------------------------------------------
MITO_N_INTERVALS   = int(config.get("mito_n_intervals", 10))
# -z: floor on total bulk reads. At mito bulk depth this is only a sanity
# guard; the germline/heteroplasmy discrimination is done by -v below.
MITO_MIN_BULK_READS = int(config.get("mito_min_bulk_reads", 100))
# -v: THE parameter to sweep. Rejects sites whose bulk VAF exceeds this, which
# is what removes germline heteroplasmies -- but a clonally expanded true
# somatic mito variant can also exceed it and be lost the same way.
MITO_MAX_BULK_VAF  = float(config.get("mito_max_bulk_vaf", 0.02))
# Sole NUMT defence while unmasked: how far the best alignment must beat the
# second-best for a read to be trusted.
MITO_MIN_AS_XS     = int(config.get("mito_min_as_xs", 50))
MITO_MIN_DPLX_PER_STRAND = int(config.get("mito_min_dplx_per_strand", 2))  # a4s2

# Noise mask for dsa -D. Unset => genuinely unmasked (Bed::Load short-circuits
# on the empty default). Point it at the empirical low-MAPQ mask to switch
# masking on; it becomes a rule input so Snakemake builds it first.
MITO_NOISE_BED = config.get("mito_noise_bed", None)

# Mito-specific wall-clock. The genome-wide helpers request 10 h for dsa and
# 15 h for indel calling, sized for 100 Mb intervals; chrM intervals are 1.7 kb
# and measured at 7 min (dsa), 12 s (varCall), 12 s (indel chain). Asking for
# 10 h x 70 jobs would sit in the queue for no reason.
def get_mito_runtime_dsa(attempt):
    return _get_dynamic_runtime(attempt, basetime=60, increment=60, label="mito_dsa")

def get_mito_runtime_varcall(attempt):
    return _get_dynamic_runtime(attempt, basetime=30, increment=30, label="mito_varCall")

def get_mito_runtime_indelcall(attempt):
    return _get_dynamic_runtime(attempt, basetime=30, increment=30, label="mito_indelCall")

def get_mito_runtime_post(attempt):
    return _get_dynamic_runtime(attempt, basetime=60, increment=60, label="mito_post")


MITO_RUNDIR   = "{donor}.runNanoSeqMito"
MITO_TMPDIR   = MITO_RUNDIR + "/tmpNanoSeq"
MITO_INTVLDIR = "mitoIntervals"   # no underscore: interval names are parsed on '_'


# ---- interval layout (resolved at parse time, so no checkpoint is needed) ---
def _mito_contig_length():
    fai = config["fasta"] + ".fai"
    with open(fai) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == MITO_CONTIG:
                return int(f[1])
    raise ValueError(f"{MITO_CONTIG} not found in {fai}")


def _mito_interval_bounds(length, n):
    step = -(-length // n)          # ceiling division
    bounds, start = [], 1
    while start <= length:
        end = min(start + step - 1, length)
        bounds.append((start, end))
        start = end + 1
    return bounds


MITO_LEN     = _mito_contig_length()
MITO_BOUNDS  = _mito_interval_bounds(MITO_LEN, MITO_N_INTERVALS)
# dsa parses chrom/start/end back out of the basename with cut -d'_', so the
# contig must not contain an underscore (chrM/MT do not).
MITO_INTERVAL_NAMES = [f"{MITO_CONTIG}_{s}_{e}" for s, e in MITO_BOUNDS]
MITO_JOBS = list(range(1, len(MITO_BOUNDS) + 1))


# ---- which donors get called -----------------------------------------------
# Selection is on the MEASURED coverage signature, never on the Fragmentation
# label or the donor suffix: two donors (ST001-1D, ST003-1Q) have their
# RENS/WGNS assignment transposed, so a label-based filter would pick their
# sparse restriction library and discard the deepest sample in the cohort.
MITO_SUMMARY_TSV   = "mito_tracks/summary.tsv"
MITO_MIN_BREADTH   = float(config.get("mito_min_breadth_pct", 99.0))
MITO_MAX_CV        = float(config.get("mito_max_cv", 0.5))
MITO_CALL_CROSS    = bool(config.get("mito_call_cross_donor", False))


def _mito_callable_donors():
    explicit = config.get("mito_call_donors", None)
    if explicit:
        return sorted(explicit)
    if not os.path.exists(MITO_SUMMARY_TSV):
        print(f"NOTE: {MITO_SUMMARY_TSV} not found -- no mito calling targets. "
              f"Run the per-base track rules first, or set config mito_call_donors.")
        return []
    summ = pd.read_csv(MITO_SUMMARY_TSV, sep="\t", header=0)
    ok = set(summ.loc[(summ["breadth_a4s2_pct"] >= MITO_MIN_BREADTH) &
                      (summ["cv_read_depth"] < MITO_MAX_CV), "sample"].values)
    tbl = pd.read_csv("input.tsv", sep="\t", header=0)
    donors = []
    for _, row in tbl.iterrows():
        if row["TestBamID"] not in ok:
            continue
        if not MITO_CALL_CROSS and "-vs-" in str(row["Donor"]):
            continue          # cross-donor normals are a control, not a callset
        donors.append(row["Donor"])
    return sorted(set(donors))


MITO_CALL_DONORS = _mito_callable_donors()
print(f"MITO CALLING: {len(MITO_CALL_DONORS)} donor(s) over {len(MITO_BOUNDS)} "
      f"interval(s) of {MITO_CONTIG} ({MITO_LEN} bp), "
      f"noise mask={MITO_NOISE_BED or 'NONE (unmasked)'}, SNP mask=NONE, "
      f"-v {MITO_MAX_BULK_VAF} -z {MITO_MIN_BULK_READS} -a {MITO_MIN_AS_XS}")
for d in MITO_CALL_DONORS:
    print(f"  - {d}")


wildcard_constraints:
    job=r"\d+",


# ---- interval + run-directory scaffolding ----------------------------------
rule mito_make_intervals:
    output:
        intervals=expand(MITO_INTVLDIR + "/{interval}.intervals.list",
                         interval=MITO_INTERVAL_NAMES),
    log:
        "logs/mito_make_intervals/all.log",
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    run:
        os.makedirs(MITO_INTVLDIR, exist_ok=True)
        for (start, end), name in zip(MITO_BOUNDS, MITO_INTERVAL_NAMES):
            with open(f"{MITO_INTVLDIR}/{name}.intervals.list", "w") as fh:
                # samtools-style region, matching the genome-wide interval files
                fh.write(f"{MITO_CONTIG}:{start}-{end}\n")
        with open(log[0], "w") as fh:
            fh.write(f"{MITO_CONTIG} {MITO_LEN} bp -> {len(MITO_BOUNDS)} intervals\n")
            for (s, e), n in zip(MITO_BOUNDS, MITO_INTERVAL_NAMES):
                fh.write(f"{n}\t{s}\t{e}\t{e - s + 1}\n")


rule mito_list_intervals:
    input:
        intvl=expand(MITO_INTVLDIR + "/{interval}.intervals.list",
                     interval=MITO_INTERVAL_NAMES),
    output:
        intvl_list=".mito_intervals.txt",
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    shell:
        # Input order is the parse-time (ascending) order, so unlike the
        # genome-wide rule this needs no sed/sort round-trip -- which would
        # have mangled the path anyway, since it splits on '_'.
        """
        printf "%s\\n" {input.intvl} > {output.intvl_list}
        """


rule mito_setup_rundir:
    # `runNanoSeq.py cov` normally creates this tree; we bypass cov, so make it here.
    input:
        duplex_bam=get_donor_expected_sample_bam,
        ctrl_bam=get_donor_expected_control_bam,
    output:
        rundir=directory(MITO_RUNDIR),
        marker=MITO_TMPDIR + "/.setup.done",
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    shell:
        """
        for d in cov part dsa var indel post; do
            mkdir -p {output.rundir}/tmpNanoSeq/$d
        done
        touch {output.marker}
        """


rule mito_start_dsa:
    input:
        intvl_list=rules.mito_list_intervals.output.intvl_list,
        marker=rules.mito_setup_rundir.output.marker,
    output:
        jobStart=expand(MITO_TMPDIR + "/dsa/{job}.start", job=MITO_JOBS,
                        allow_missing=True),
        jobToIntvl=MITO_TMPDIR + "/dsa/job_to_intvl.txt",
    params:
        jobs=MITO_JOBS,
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    shell:
        """
        : > {output.jobToIntvl}
        for i in {params.jobs}; do
            intvl=$(head -$i {input.intvl_list} | tail -1)
            echo -e "$intvl" > {wildcards.donor}.runNanoSeqMito/tmpNanoSeq/dsa/$i.start
            echo -e "$i\t$intvl" >> {output.jobToIntvl}
        done
        """


rule mito_add_dsa_args:
    input:
        intvl_list=rules.mito_list_intervals.output.intvl_list,
        marker=rules.mito_setup_rundir.output.marker,
    output:
        nfiles=MITO_TMPDIR + "/dsa/nfiles",
        argsJson=MITO_TMPDIR + "/dsa/args.json",
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """


# ---- dsa (unmasked) --------------------------------------------------------
rule mito_dsa_bed_per_partition:
    input:
        # NB: no donor sheet input. It was vestigial once masking moved to
        # config (mito_noise_bed) -- the shell never read it -- and depending on
        # rules.list_donor_info tied this module to the by-donor entry point,
        # breaking `include:` from the sample-keyed Snakefile.
        duplex_bam=get_donor_expected_sample_bam,
        ctrl_bam=get_donor_expected_control_bam,
        job=MITO_TMPDIR + "/dsa/{job}.start",
        allJobs=rules.mito_start_dsa.output.jobToIntvl,
        # dsa reads the mask through tabix, so the index is a real dependency.
        noise=[MITO_NOISE_BED, MITO_NOISE_BED + ".tbi"] if MITO_NOISE_BED else [],
        nfiles=rules.mito_add_dsa_args.output.nfiles,
        argsJson=rules.mito_add_dsa_args.output.argsJson,
    output:
        dsa_bed=MITO_TMPDIR + "/dsa/{job}.dsa.bed.gz",
        jobDone=MITO_TMPDIR + "/dsa/{job}.done",
    benchmark:
        "benchmarks/mito_dsa_bed_per_partition/{donor}.{job}.txt"
    log:
        "logs/mito_dsa_bed_per_partition/{donor}.{job}.log"
    params:
        fasta=config["fasta"],
        rundir=MITO_RUNDIR,
        # dsa runs from inside the rundir, so relative paths need ../
        noise_arg=("-D ../" + MITO_NOISE_BED) if MITO_NOISE_BED else "",
    resources:
        # Measured on chrM:1-1657 at ~1,500x duplex / ~16,000x bulk: dsa peaked
        # at 85 MB and 6m54s. It streams per position, so depth drives runtime,
        # not memory -- the initial 24 GB guess was ~280x over.
        mem_mb=4000,
        runtime=lambda wildcards, attempt: get_mito_runtime_dsa(attempt),
    threads: 1
    group:
        "mito_dsa"
    shell:
        """
        PATH=$PATH:$PWD/bin/

        intvl=$(head -1 {input.job})
        intvl=$(basename $intvl)
        intvl=${{intvl%.intervals*}}
        chrom=$(echo -e $intvl | cut -d'_' -f1)
        startPos=$(echo -e $intvl | cut -d'_' -f2)
        endPos=$(echo -e $intvl | cut -d'_' -f3)
        echo -e "dsa job {wildcards.job}: $chrom:$startPos-$endPos (UNMASKED)"

        cd {params.rundir}

        # -C (SNP mask) is still omitted: no mito germline mask exists yet, and
        # Bed::Load short-circuits on the empty default. -D carries the
        # empirical low-MAPQ noise mask when config mito_noise_bed is set.
        dsa \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        {params.noise_arg} \
        -R {params.fasta} \
        -d 2 \
        -Q 30 \
        -M 0 \
        -r $chrom \
        -b $startPos \
        -e $endPos \
        -O ./tmpNanoSeq/dsa/{wildcards.job}.dsa.bed
        if [ -f ./tmpNanoSeq/dsa/{wildcards.job}.dsa.bed.gz ]; then
            touch ./tmpNanoSeq/dsa/{wildcards.job}.done
        else
            echo -e "dsa failed for job {wildcards.job} over $intvl"
            exit 1
        fi
        """


# ---- SNV calling -----------------------------------------------------------
rule mito_start_varCall:
    input:
        intvl_list=rules.mito_list_intervals.output.intvl_list,
        nfiles=rules.mito_add_dsa_args.output.nfiles,
    output:
        nfiles=MITO_TMPDIR + "/var/nfiles",
        argsJson=MITO_TMPDIR + "/var/args.json",
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """


mitoVarPath = MITO_TMPDIR + "/var/{job}"

rule mito_varCall_per_partition:
    input:
        dsa=rules.mito_dsa_bed_per_partition.output.dsa_bed,
        nfiles=rules.mito_start_varCall.output.nfiles,
        duplex_bam=get_donor_expected_sample_bam,
        ctrl_bam=get_donor_expected_control_bam,
        job=MITO_TMPDIR + "/dsa/{job}.start",
    output:
        coverage=mitoVarPath + ".cov.bed.gz",
        var=mitoVarPath + ".var",
        discarded_var=mitoVarPath + ".discarded_var",
        doneFile=mitoVarPath + ".done",
    benchmark:
        "benchmarks/mito_varCall_per_partition/{donor}.{job}.txt"
    log:
        "logs/mito_varCall_per_partition/{donor}.{job}.log"
    params:
        rundir=MITO_RUNDIR,
        min_as_xs=MITO_MIN_AS_XS,
        min_bulk_reads_per_strand=0,
        max_frac_clips=0,
        min_number_dplx_reads_per_strand=MITO_MIN_DPLX_PER_STRAND,
        min_fraction_reads_consensus=0.9,
        max_frac_reads_w_indel=1.0,
        min_cycle_num=8,
        max_num_mismatches=3,
        min_fraction_reads_properPairs=0,
        min_consensus_base_quality=60,
        read_length_post_5prime_trimming=144,
        max_bulk_vaf=MITO_MAX_BULK_VAF,
        max_cycle_num=8,
        # Explicit, not read from the bypassed coverage diagnostic.
        min_num_bulkReads_total=MITO_MIN_BULK_READS,
    resources:
        # Measured: 11 MB peak, 12 s for the same interval.
        mem_mb=2000,
        runtime=lambda wildcards, attempt: get_mito_runtime_varcall(attempt),
    threads: 1
    group:
        "mito_varCall"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        cd {params.rundir}
        ulimit -c unlimited
        ulimit -s unlimited

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
        touch ./tmpNanoSeq/var/{wildcards.job}.done
        """


# ---- indel calling ---------------------------------------------------------
mitoIndelPath = MITO_TMPDIR + "/indel"
mitoIndelJob  = mitoIndelPath + "/{job}"

rule mito_start_indelCall:
    input:
        intvl_list=rules.mito_list_intervals.output.intvl_list,
        nfiles=rules.mito_add_dsa_args.output.nfiles,
    output:
        nfiles=mitoIndelPath + "/nfiles",
        argsJson=mitoIndelPath + "/args.json",
    resources:
        mem_mb=1000, runtime=5,
    threads: 1
    shell:
        """
        grep -c "^" {input.intvl_list} > {output.nfiles}
        touch {output.argsJson}
        """


rule mito_indel_propose_per_partition:
    input:
        dsa=rules.mito_dsa_bed_per_partition.output.dsa_bed,
        nfiles=rules.mito_start_indelCall.output.nfiles,
        argsJson=rules.mito_start_indelCall.output.argsJson,
    output:
        indel_bed=mitoIndelJob + ".indel.bed.gz",
    benchmark:
        "benchmarks/mito_indel_propose_per_partition/{donor}.{job}.txt"
    log:
        "logs/mito_indel_propose_per_partition/{donor}.{job}.log"
    params:
        rundir=MITO_RUNDIR,
        max_reads_bundle=2,
        min_as_xs=MITO_MIN_AS_XS,
        min_normal_coverage=15,
        max_bulk_vaf=MITO_MAX_BULK_VAF,
        max_frac_clips=0,
        trim_from_3p_at_pos=135,
        trim_from_5p_at_pos=10,
    resources:
        mem_mb=8000, runtime=60,
    threads: 1
    group:
        "mito_indel_propose"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        cd {params.rundir}
        indelCaller.py propose \
        -o ./tmpNanoSeq/indel/{wildcards.job}.indel.bed.gz \
        -rb {params.max_reads_bundle} \
        -t3 {params.trim_from_3p_at_pos} \
        -t5 {params.trim_from_5p_at_pos} \
        -mc {params.min_normal_coverage} \
        -vaf {params.max_bulk_vaf} \
        -a {params.min_as_xs} \
        -c {params.max_frac_clips} \
        ./tmpNanoSeq/dsa/{wildcards.job}.dsa.bed.gz
        """


rule mito_indel_call_per_partition:
    input:
        duplex_bam=get_donor_expected_sample_bam,
        indel_bed=rules.mito_indel_propose_per_partition.output.indel_bed,
    output:
        indel_vcf=mitoIndelJob + ".indel.vcf.gz",
    benchmark:
        "benchmarks/mito_indel_call_per_partition/{donor}.{job}.txt"
    log:
        "logs/mito_indel_call_per_partition/{donor}.{job}.log"
    params:
        fasta=config["fasta"],
        rundir=MITO_RUNDIR,
    resources:
        mem_mb=12000,
        runtime=lambda wildcards, attempt: get_mito_runtime_indelcall(attempt),
    threads: 1
    group:
        "mito_indel_call"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        cd {params.rundir}
        indelCaller.py call \
        -t \
        -o ./tmpNanoSeq/indel/{wildcards.job}.indel \
        -r {params.fasta} \
        -b ../{input.duplex_bam} \
        ./tmpNanoSeq/indel/{wildcards.job}.indel.bed.gz
        """


rule mito_indel_verify_per_partition:
    input:
        ctrl_bam=get_donor_expected_control_bam,
        indel_vcf=rules.mito_indel_call_per_partition.output.indel_vcf,
    output:
        indel_filtered_vcf=mitoIndelJob + ".indel.filtered.vcf.gz",
        indel_filtered_vcf_tbi=mitoIndelJob + ".indel.filtered.vcf.gz.tbi",
        doneFile=mitoIndelJob + ".done",
    benchmark:
        "benchmarks/mito_indel_verify_per_partition/{donor}.{job}.txt"
    log:
        "logs/mito_indel_verify_per_partition/{donor}.{job}.log"
    params:
        fasta=config["fasta"],
        rundir=MITO_RUNDIR,
        max_bulk_vaf=MITO_MAX_BULK_VAF,
    resources:
        mem_mb=12000, runtime=60 * 2,
    threads: 1
    group:
        "mito_indel_verify"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        cd {params.rundir}
        indelCaller.py verify \
        {params.fasta} \
        ./tmpNanoSeq/indel/{wildcards.job}.indel.vcf.gz \
        ../{input.ctrl_bam} \
        {params.max_bulk_vaf} &&
        touch ./tmpNanoSeq/indel/{wildcards.job}.done
        """


# ---- post ------------------------------------------------------------------
mitoPostPath = MITO_TMPDIR + "/post/"

rule mito_post:
    input:
        duplex_bam=get_donor_expected_sample_bam,
        ctrl_bam=get_donor_expected_control_bam,
        var_done=expand(mitoVarPath + ".done", job=MITO_JOBS, allow_missing=True),
        indel_done=expand(mitoIndelJob + ".done", job=MITO_JOBS, allow_missing=True),
        nfiles=rules.mito_start_indelCall.output.nfiles,
    output:
        variants=mitoPostPath + "results.muts.vcf.gz",
        doneFile=mitoPostPath + "1.done",
    benchmark:
        "benchmarks/mito_post/{donor}.txt"
    log:
        "logs/mito_post/{donor}.log"
    params:
        fasta=config["fasta"],
        rundir=MITO_RUNDIR,
    resources:
        mem_mb=8000,
        runtime=lambda wildcards, attempt: get_mito_runtime_post(attempt),
    threads: 4
    group:
        "mito_post"
    shell:
        """
        PATH=$PATH:$PWD/bin/
        cd {params.rundir}
        mkdir -p ./tmpNanoSeq/post
        touch ./tmpNanoSeq/post/args.json
        runNanoSeq.py \
        -t {threads} \
        -A ../{input.ctrl_bam} \
        -B ../{input.duplex_bam} \
        -R {params.fasta} \
        post
        """


rule mito_annotate_post_csv_with_interval:
    input:
        variants=rules.mito_post.output.variants,
    output:
        annotated=expand(mitoPostPath + "{summary_results}.annotated.tsv",
                         summary_results=SUMMARY_RESULT_FILES, allow_missing=True),
        doneFile=mitoPostPath + "1.annotate_interval.done",
    benchmark:
        "benchmarks/mito_annotate_post/{donor}.txt"
    log:
        "logs/mito_annotate_post/{donor}.log"
    params:
        rundir=MITO_RUNDIR,
    resources:
        mem_mb=4000, runtime=30,
    threads: 1
    group:
        "mito_post"
    shell:
        """
        cd {params.rundir}
        postDir="./tmpNanoSeq/post"
        varDir="./tmpNanoSeq/var"
        dsaDir="./tmpNanoSeq/dsa"

        for f in burdens callvsqpos coverage pyrvsmask readbundles; do
            paste <(head -1 $postDir/$f.csv) <(echo -e "intvl") | sed "s/,/\t/g" \
                > $postDir/$f.annotated.tsv
        done

        for i in $varDir/*.var; do
            intvlChunk=$(basename $i); intvlChunk=${{intvlChunk%.var}}
            dsaIntvlChunk=$(cat $dsaDir/$intvlChunk.start)
            for f in Burdens:burdens CallVsQpos:callvsqpos Coverage:coverage \
                     PyrVsMask:pyrvsmask ReadBundles:readbundles; do
                tag=${{f%%:*}}; out=${{f##*:}}
                grep -s "$tag" $i | cut -f2- \
                  | awk -v chunk=$dsaIntvlChunk 'BEGIN {{FS=OFS="\t"}} {{print $0,chunk}}' \
                  >> $postDir/$out.annotated.tsv
            done
        done

        touch $postDir/1.annotate_interval.done
        """


# Final targets for the mito calling pass. Exposed as a function because
# `rule all` in the parent Snakefile is defined BEFORE this file is included,
# so MITO_CALL_DONORS does not exist yet at that point; referencing this
# lazily (via a lambda in rule all) defers evaluation to DAG-build time.
# ---- consolidation: the terminal step of a calling run ---------------------
# Rolls every donor's post/ output into one call table, one per-sample burden
# table, and NB rate estimates. Built as a rule rather than a script run by
# hand afterwards so the consolidated tables are guaranteed to correspond to
# the run that produced them, and so a scaled-up cohort needs no manual step.
MITO_CONSOLIDATE_PREFIX = config.get("mito_consolidate_prefix", "mito_consolidated/cohort")

rule mito_consolidate:
    input:
        muts=lambda wildcards: expand(MITO_TMPDIR + "/post/results.muts.vcf.gz",
                                      donor=MITO_CALL_DONORS),
        burden=lambda wildcards: expand(MITO_TMPDIR + "/post/results.mut_burden.tsv",
                                        donor=MITO_CALL_DONORS),
    output:
        calls=MITO_CONSOLIDATE_PREFIX + ".calls.tsv.gz",
        burden=MITO_CONSOLIDATE_PREFIX + ".burden.tsv",
        rates=MITO_CONSOLIDATE_PREFIX + ".rates.tsv",
    benchmark:
        "benchmarks/mito_consolidate/all.txt"
    log:
        "logs/mito_consolidate/all.log",
    params:
        prefix=MITO_CONSOLIDATE_PREFIX,
        # Contaminated libraries are NOT dropped from the call/burden tables --
        # they are reported with their diagnostics so exclusion stays a
        # documented downstream decision. They are excluded only from the NB
        # rate fits, where a single contaminated library can dominate a group.
        exclude=",".join(config.get("mito_exclude_samples", [])),
    resources:
        mem_mb=8000,
        runtime=30,
    threads: 1
    conda:
        workflow.source_path("../envs/nanoseq_snakemake.yml")
    shell:
        """
        mkdir -p $(dirname {params.prefix}) logs/mito_consolidate
        python3 bin/compile_mito_master_table.py -d . -o {params.prefix} > {log} 2>&1
        Rscript bin/mito_burden_nb.R {output.burden} {output.rates} "{params.exclude}" >> {log} 2>&1
        """


def mito_call_targets(wildcards=None):
    if not MITO_CALL_DONORS:
        return []
    targets = expand(MITO_TMPDIR + "/post/results.muts.vcf.gz",
                     donor=MITO_CALL_DONORS)
    targets += expand(MITO_TMPDIR + "/post/{summary_results}.annotated.tsv",
                      donor=MITO_CALL_DONORS,
                      summary_results=SUMMARY_RESULT_FILES)
    targets += [MITO_CONSOLIDATE_PREFIX + ".calls.tsv.gz",
                MITO_CONSOLIDATE_PREFIX + ".burden.tsv",
                MITO_CONSOLIDATE_PREFIX + ".rates.tsv"]
    return targets
