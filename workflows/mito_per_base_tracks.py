#!/usr/bin/env python3
"""Per-base coverage tracks over the mitochondrial contig.

Emits one row per base of the mito contig with both read-level and
read-bundle-level depth, so duplex coverage can be judged on the quantity that
actually matters: how many distinct molecules (read bundles) cover each base,
not how many reads.

Why bundle depth is derived from the RB tag rather than from a deduplicated
BAM: the RB tag carries the *fragment* interval, so a bundle contributes
coverage across its whole insert. Counting one representative read pair per
bundle instead would leave the unsequenced middle of long inserts uncovered.
The RB tag is already in BED coordinates -- RB:Z:chrM,0,141 corresponds to a
read at POS=1 with TLEN=141 -- so intervals are used verbatim, no shifting.

All accumulation is via difference arrays, so each BAM is a single O(n_reads)
pass with no per-base inner loop.
"""

import argparse
import gzip
import sys

import numpy as np
import pysam

# Reads that should never contribute to any track.
_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x800  # unmapped, secondary, qcfail, supplementary


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-bam", required=True,
                   help="read-bundled BAM for the sample (all reads, duplicates retained)")
    p.add_argument("--a4s2-bam", required=True,
                   help="BAM restricted to reads of a4s2-qualifying duplex bundles")
    p.add_argument("--contig", default="chrM", help="mitochondrial contig name")
    p.add_argument("--sample", required=True, help="sample label emitted in column 1")
    p.add_argument("--fragmentation", default="NA", help="fragmentation protocol label")
    p.add_argument("--min-mapq", type=int, default=30,
                   help="reads below this MAPQ count toward the ambiguity track")
    p.add_argument("--out", required=True, help="output .tsv.gz")
    return p.parse_args()


def _contig_length(bam, contig):
    try:
        return bam.lengths[bam.references.index(contig)]
    except ValueError:
        sys.exit(f"ERROR: contig {contig!r} not present in {bam.filename.decode()}")


def _add(diff, start, end, length):
    """Accumulate a half-open [start, end) interval into a difference array."""
    start = max(int(start), 0)
    end = min(int(end), length)
    if end > start:
        diff[start] += 1
        diff[end] -= 1


def _scan_reads(path, contig, length, min_mapq):
    """One pass over the contig: read depth, low-MAPQ depth, and the RB tag set."""
    depth = np.zeros(length + 1, dtype=np.int64)
    lowmq = np.zeros(length + 1, dtype=np.int64)
    bundles = set()
    with pysam.AlignmentFile(path, "r") as bam:
        for read in bam.fetch(contig):
            if read.flag & _SKIP_FLAGS:
                continue
            start, end = read.reference_start, read.reference_end
            if end is None:
                continue
            _add(depth, start, end, length)
            if read.mapping_quality < min_mapq:
                _add(lowmq, start, end, length)
            if read.has_tag("RB"):
                bundles.add(read.get_tag("RB"))
    return depth, lowmq, bundles


def _bundle_depth(bundles, contig, length):
    """Per-base count of distinct read bundles, from the fragment interval in each RB tag."""
    diff = np.zeros(length + 1, dtype=np.int64)
    skipped = 0
    for rb in bundles:
        fields = rb.split(",")
        if len(fields) < 3 or fields[0] != contig:
            skipped += 1
            continue
        try:
            _add(diff, int(fields[1]), int(fields[2]), length)
        except ValueError:
            skipped += 1
    if skipped:
        print(f"  note: {skipped} RB tag(s) unparseable or off-contig, skipped",
              file=sys.stderr)
    return diff


def main():
    args = _parse_args()

    with pysam.AlignmentFile(args.raw_bam, "r") as bam:
        length = _contig_length(bam, args.contig)
    print(f"{args.sample}: {args.contig} is {length} bp", file=sys.stderr)

    raw_depth_d, lowmq_d, raw_bundles = _scan_reads(
        args.raw_bam, args.contig, length, args.min_mapq)
    _, _, a4s2_bundles = _scan_reads(
        args.a4s2_bam, args.contig, length, args.min_mapq)
    print(f"  bundles: {len(raw_bundles)} total, {len(a4s2_bundles)} a4s2", file=sys.stderr)

    rb_all_d = _bundle_depth(raw_bundles, args.contig, length)
    rb_a4s2_d = _bundle_depth(a4s2_bundles, args.contig, length)

    read_depth = np.cumsum(raw_depth_d)[:length]
    lowmq_depth = np.cumsum(lowmq_d)[:length]
    rb_all = np.cumsum(rb_all_d)[:length]
    rb_a4s2 = np.cumsum(rb_a4s2_d)[:length]

    # Guard the denominators: uncovered bases give an undefined ratio, not zero.
    with np.errstate(divide="ignore", invalid="ignore"):
        duplex_frac = np.where(rb_all > 0, rb_a4s2 / rb_all, np.nan)
        lowmq_frac = np.where(read_depth > 0, lowmq_depth / read_depth, np.nan)

    with gzip.open(args.out, "wt") as out:
        out.write("sample\tfragmentation\tpos\tread_depth_raw\trb_depth_all"
                  "\trb_depth_a4s2\tduplex_frac\tlowmq_frac\n")
        for i in range(length):
            df = duplex_frac[i]
            lf = lowmq_frac[i]
            out.write(
                f"{args.sample}\t{args.fragmentation}\t{i + 1}\t{read_depth[i]}"
                f"\t{rb_all[i]}\t{rb_a4s2[i]}"
                f"\t{'NA' if np.isnan(df) else f'{df:.6g}'}"
                f"\t{'NA' if np.isnan(lf) else f'{lf:.6g}'}\n")

    covered = int((rb_a4s2 > 0).sum())
    print(f"  a4s2 bundle depth: mean {rb_a4s2.mean():.2f}, "
          f"breadth {100.0 * covered / length:.2f}%", file=sys.stderr)


if __name__ == "__main__":
    main()
