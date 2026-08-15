#!/usr/bin/env python3
"""Derive an empirical noise mask for the mitochondrial contig from per-base
low-MAPQ fraction.

This is a MAPPABILITY mask, not a germline mask and not literally a NUMT
catalogue. It marks where reads are placed ambiguously, which is the main
symptom of NUMT homology but also picks up ordinary repeats and low-complexity
sequence. Treat it as a first-pass stand-in until a sequence-derived NUMT set
or an assembly-based germline mask is available.

Consensus is taken as the MEDIAN across samples so that one aberrant library
cannot create or erase an interval: a region is masked only if it is
ambiguous in most samples, which is what makes it a property of the reference
rather than of a library.

Caveat worth stating plainly: the mask is derived from the same data it will
be applied to. Because the signal is near-identical across independent
samples, that is defensible, but it is not independent validation.
"""

import argparse
import sys

import numpy as np
import pandas as pd


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tracks", nargs="+", required=True,
                   help="per-base track TSVs (one per sample)")
    p.add_argument("--contig", default="chrM")
    p.add_argument("--threshold", type=float, default=0.20,
                   help="median low-MAPQ fraction above which a base is masked")
    p.add_argument("--min-depth", type=int, default=20,
                   help="ignore bases below this raw read depth: a fraction over "
                        "a tiny denominator is noise, not a measurement")
    p.add_argument("--merge-gap", type=int, default=15,
                   help="join flagged runs separated by at most this many bases")
    p.add_argument("--min-length", type=int, default=15,
                   help="drop intervals shorter than this")
    p.add_argument("--min-breadth", type=float, default=0.99,
                   help="only use samples covering at least this fraction of the "
                        "contig at --min-depth; a restriction-digested library "
                        "has no coverage over most of chrM and so cannot report "
                        "on mappability there")
    p.add_argument("--out", required=True, help="output BED (uncompressed)")
    return p.parse_args()


def main():
    args = _parse_args()

    # Sample eligibility is decided here, from the tracks themselves, rather
    # than upstream: routing it through the donor table silently dropped four
    # of seven deep libraries, because samples shared by a cross-donor row map
    # to a donor that the calling workflow excludes.
    cols, dropped = {}, []
    for path in args.tracks:
        d = pd.read_csv(path, sep="\t")
        name = d["sample"].iloc[0]
        breadth = float((d["read_depth_raw"] >= args.min_depth).mean())
        if breadth < args.min_breadth:
            dropped.append((name, breadth))
            continue
        d.loc[d["read_depth_raw"] < args.min_depth, "lowmq_frac"] = np.nan
        cols[name] = d.set_index("pos")["lowmq_frac"]
    for name, breadth in dropped:
        print(f"  skipping {name}: breadth {100 * breadth:.1f}% < "
              f"{100 * args.min_breadth:.0f}%", file=sys.stderr)
    for name in cols:
        print(f"  using    {name}", file=sys.stderr)
    if not cols:
        sys.exit("no eligible tracks: all samples below --min-breadth")
    med = pd.DataFrame(cols).median(axis=1)
    print(f"{len(cols)} sample(s); median low-MAPQ fraction baseline "
          f"{med.median():.4f}, max {med.max():.3f}", file=sys.stderr)

    flagged = np.where((med > args.threshold).values)[0] + 1
    if not len(flagged):
        print("WARNING: no bases exceeded the threshold; writing empty mask",
              file=sys.stderr)

    runs, start, prev = [], None, None
    for pos in flagged:
        if start is None:
            start = prev = pos
            continue
        if pos - prev > args.merge_gap:
            runs.append((start, prev))
            start = pos
        prev = pos
    if start is not None:
        runs.append((start, prev))
    runs = [(a, b) for a, b in runs if b - a + 1 >= args.min_length]

    total = 0
    with open(args.out, "w") as out:
        for a, b in runs:
            total += b - a + 1
            peak = float(med[a:b + 1].max())
            # BED is 0-based half-open; positions here are 1-based inclusive.
            out.write(f"{args.contig}\t{a - 1}\t{b}\tlowmapq_peak_{peak:.2f}\t"
                      f"{int(round(1000 * peak))}\n")

    span = len(med)
    print(f"{len(runs)} interval(s), {total} bp masked "
          f"({100.0 * total / span:.2f}% of {args.contig}) at threshold "
          f"{args.threshold}", file=sys.stderr)


if __name__ == "__main__":
    main()
