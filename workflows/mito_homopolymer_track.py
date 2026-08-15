#!/usr/bin/env python3
"""Homopolymer runs on the mitochondrial contig, as a stratifier for indel calls.

Polymerase slippage in homopolymer tracts is the dominant source of spurious
mtDNA indels, and it concentrates in a very small footprint: runs of >=5 bp
cover only ~3.7% of chrM. Any claim about mitochondrial indel burden should
therefore be reported stratified by this track, because a global indel:SNV
ratio is dominated by whether homopolymer calls passed filtering.

Duplex sequencing is unusually well suited to settling this: a slippage event
arising during library PCR would have to occur identically on both strands to
survive duplex consensus, so indels called here from a4s2 bundles are far more
credible than single-strand calls. The track exists to make that comparison
explicit rather than assumed.

BED score carries the run length; name carries base and length (e.g. C7).
"""

import argparse
import sys

import pysam


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--contig", default="chrM")
    p.add_argument("--min-run", type=int, default=5,
                   help="minimum homopolymer length to report")
    p.add_argument("--pad", type=int, default=1,
                   help="bases of flank to include; an indel at a run edge is "
                        "ambiguous in placement, so the run alone under-captures")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = _parse_args()
    seq = pysam.FastaFile(args.fasta).fetch(args.contig).upper()
    n = len(seq)

    runs, i = [], 0
    while i < n:
        j = i
        while j + 1 < n and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= args.min_run and seq[i] in "ACGT":
            runs.append((i + 1, j + 1, seq[i], j - i + 1))
        i = j + 1

    total = 0
    with open(args.out, "w") as out:
        for start, end, base, length in runs:
            s = max(start - args.pad, 1)
            e = min(end + args.pad, n)
            total += e - s + 1
            out.write(f"{args.contig}\t{s - 1}\t{e}\t{base}{length}\t{length}\n")

    print(f"{len(runs)} homopolymer run(s) >= {args.min_run} bp", file=sys.stderr)
    print(f"  padded footprint: {total} bp ({100.0 * total / n:.2f}% of "
          f"{args.contig}, {n} bp)", file=sys.stderr)
    print(f"  => under a null of uniformly distributed indels, this is the "
          f"fraction expected to fall in homopolymers", file=sys.stderr)


if __name__ == "__main__":
    main()
