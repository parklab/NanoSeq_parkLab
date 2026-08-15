#!/usr/bin/env python3
"""Derive a NUMT annotation track for the mitochondrial contig.

Tiles the mito contig into overlapping windows, aligns them to the nuclear
contigs, and merges collinear hits into NUMT blocks reported in MITO
coordinates.

This is an ANNOTATION, deliberately not a mask. Roughly 80% of chrM has some
NUMT of >=200 bp, so NUMT presence carries almost no information for masking;
what predicts mismapping is sequence IDENTITY, and sharply so. Measured
against the observed low-MAPQ fraction:

    best nuclear identity   median low-MAPQ fraction
    < 97%                   0.0002 - 0.0005   (baseline)
    97 - 99%                0.0045
    99 - 99.5%              0.0833
    >= 99.5%                0.3054

Below ~97% there are enough mismatches per read to place it unambiguously, so
a 5 kb NUMT at 95% identity is harmless. Use the identity column here to
reason about which blocks matter; use the low-MAPQ mask for actual masking.

Sensitivity caveat: alignment is by bwa mem, tuned for short reads, so highly
diverged NUMTs (below roughly 85-90%) will be under-reported. That is
acceptable for this purpose precisely because those are the ones that do not
cause mismapping.
"""

import argparse
import re
import statistics
import sys


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sam", required=True,
                   help="SAM of tiled mito windows aligned to the whole genome")
    p.add_argument("--contig", default="chrM")
    p.add_argument("--min-aln", type=int, default=80,
                   help="ignore alignments shorter than this")
    p.add_argument("--min-block", type=int, default=200,
                   help="report blocks spanning at least this much of the contig")
    p.add_argument("--offset-tol", type=int, default=50,
                   help="collinearity tolerance when merging windows into a block")
    p.add_argument("--join-gap", type=int, default=600,
                   help="join collinear windows separated by at most this much")
    p.add_argument("--out", required=True, help="output BED (uncompressed)")
    return p.parse_args()


def main():
    args = _parse_args()
    hits = []
    contig_len = None
    for line in open(args.sam):
        if line.startswith("@"):
            if line.startswith("@SQ") and f"SN:{args.contig}\t" in line + "\t":
                m = re.search(r"LN:(\d+)", line)
                if m:
                    contig_len = int(m.group(1))
            continue
        f = line.rstrip("\n").split("\t")
        rname = f[2]
        if rname in ("*", args.contig):
            continue
        flag, pos = int(f[1]), int(f[3])
        nm = next((int(t.split(":")[2]) for t in f[11:] if t.startswith("NM:i:")), 0)
        alen = sum(int(n) for n, op in re.findall(r"(\d+)([MID])", f[5]) if op == "M")
        if alen < args.min_aln:
            continue
        qstart = int(f[0].split(":")[1].split("-")[0])
        rev = bool(flag & 16)
        # Collinearity key: for a forward hit the nuclear-minus-mito offset is
        # constant along a NUMT; for a reverse hit the sum is.
        offset = pos - qstart if not rev else pos + qstart
        hits.append((rname, rev, offset, qstart, qstart + alen - 1, pos, alen,
                     100.0 * (1 - nm / max(alen, 1))))

    if not hits:
        sys.exit("no nuclear alignments found")

    hits.sort(key=lambda h: (h[0], h[1], h[2], h[3]))
    blocks = []
    for rname, rev, offset, qs, qe, rpos, alen, ident in hits:
        for b in blocks:
            if (b["rname"] == rname and b["rev"] == rev
                    and abs(b["offset"] - offset) <= args.offset_tol
                    and qs <= b["qe"] + args.join_gap):
                b["qs"], b["qe"] = min(b["qs"], qs), max(b["qe"], qe)
                b["rs"] = min(b["rs"], rpos)
                b["re"] = max(b["re"], rpos + alen - 1)
                b["ident"].append(ident)
                break
        else:
            blocks.append(dict(rname=rname, rev=rev, offset=offset, qs=qs, qe=qe,
                               rs=rpos, re=rpos + alen - 1, ident=[ident]))

    blocks = [b for b in blocks if b["qe"] - b["qs"] + 1 >= args.min_block]
    blocks.sort(key=lambda b: b["qs"])

    covered = set()
    with open(args.out, "w") as out:
        for b in blocks:
            covered.update(range(b["qs"], b["qe"] + 1))
            top = max(b["ident"])
            med = statistics.median(b["ident"])
            name = (f"NUMT_{b['rname']}:{b['rs']}-{b['re']}"
                    f"_maxid{top:.1f}_medid{med:.1f}")
            out.write(f"{args.contig}\t{b['qs'] - 1}\t{b['qe']}\t{name}\t"
                      f"{int(round(10 * top))}\t{'-' if b['rev'] else '+'}\n")

    # Percentages are of the whole contig. Using the last NUMT end as the
    # denominator (as an earlier version did) inflates them: it silently
    # excludes the NUMT-free tail of the contig from the total.
    span = contig_len or max(b["qe"] for b in blocks)
    hi = {p for b in blocks if max(b["ident"]) >= 99
          for p in range(b["qs"], b["qe"] + 1)}
    print(f"{len(blocks)} NUMT block(s) >= {args.min_block} bp", file=sys.stderr)
    print(f"  contig bases with any NUMT : {len(covered)} "
          f"({100.0 * len(covered) / max(span, 1):.1f}% of {span} bp contig)",
          file=sys.stderr)
    print(f"  ... with >=99% identity    : {len(hi)} "
          f"({100.0 * len(hi) / max(span, 1):.1f}%)  <- the subset that mismaps",
          file=sys.stderr)


if __name__ == "__main__":
    main()
