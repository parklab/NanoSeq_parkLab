#!/usr/bin/env bash
# Snakemake --cluster-status hook: report a Slurm job as success/failed/running.
#
# Why this exists: with plain --cluster, Snakemake infers completion from a
# marker file the job wrapper writes on exit. A job killed by the scheduler
# (SIGTERM, node failure, preemption) writes NEITHER jobfinished NOR jobfailed,
# so Snakemake waits on it forever -- no error, no retry, no exit. That cost a
# 4h40m silent hang on 2026-08-15 over a single cancelled sort_mito_bam_by_rb.
# Polling sacct instead makes such a job visibly failed, so --retries applies.
#
# Contract: receives the external job id as $1, prints exactly one of
# "success", "failed", "running".

set -uo pipefail

# Snakemake passes back whatever `sbatch` printed. Without --parsable that is
# "Submitted batch job 12345", so take the trailing integer either way.
jobid="$(printf '%s' "${1:-}" | grep -oE '[0-9]+$')"
if [[ -z "$jobid" ]]; then
    echo failed
    exit 0
fi

state="$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null \
         | head -n1 | tr -d ' ')"

# Empty means sacct has not registered the job yet (or slurmdbd hiccuped).
# Report running: a false "failed" here would kill a healthy job, whereas a
# false "running" only delays detection until the next poll.
if [[ -z "$state" ]]; then
    echo running
    exit 0
fi

case "$state" in
    COMPLETED)
        echo success ;;
    PENDING|RUNNING|COMPLETING|CONFIGURING|SUSPENDED|REQUEUED|RESIZING)
        echo running ;;
    *)
        # FAILED, CANCELLED*, TIMEOUT, NODE_FAIL, PREEMPTED, OUT_OF_MEMORY,
        # BOOT_FAIL, DEADLINE, REVOKED -- all terminal and unsuccessful.
        echo failed ;;
esac
