#!/usr/bin/env python3
"""indelCaller.py -- Python port of the Sanger NanoSeq three-stage indel caller.

The classical NanoSeq indel caller is a chain of three programs:

    indelCaller_step1.pl  ->  indelCaller_step2.pl  ->  indelCaller_step3.R

This script reimplements that chain as one readable program with three
subcommands, named for their role (the classical step they replace is noted
in each subcommand's help):

    propose   (= classical indelCaller_step1.pl)
              Filter the per-site `dsa.bed.gz` table to candidate indel sites.
              Default reproduces the canonical double-stranded caller. Flags add
              the single-stranded-indel behaviour and a rich all-candidate emitter.

    call      (= classical indelCaller_step2.pl)
              Group candidates by read bundle, extract the reads, and call indels
              per bundle with samtools/bcftools.

    verify    (= classical indelCaller_step3.R, ported to pysam)
              Check each call against the matched normal and flag indel-rich or
              missing sites.

Each subcommand keeps the flag names of its perl/R original, so existing
workflows only need to swap the executable and add the subcommand token.
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# dsa.bed.gz column layout (0-based indices), matching indelCaller_step1.pl:70-78
# ---------------------------------------------------------------------------
COLS = [
    "chrom", "chromBeg", "chromEnd", "context", "commonSNP", "shearwater",
    "bulkASXS", "bulkNM",
    "bulkForwardA", "bulkForwardC", "bulkForwardG", "bulkForwardT", "bulkForwardIndel",
    "bulkReverseA", "bulkReverseC", "bulkReverseG", "bulkReverseT", "bulkReverseIndel",
    "dplxBreakpointBeg", "dplxBreakpointEnd", "dplxBarcode", "dplxOri",
    "dplxASXS", "dplxCLIP", "dplxNM",
    "dplxForwardA", "dplxForwardC", "dplxForwardG", "dplxForwardT", "dplxForwardIndel",
    "dplxReverseA", "dplxReverseC", "dplxReverseG", "dplxReverseT", "dplxReverseIndel",
    "dplxCQForwardA", "dplxCQForwardC", "dplxCQForwardG", "dplxCQForwardT",
    "dplxCQReverseA", "dplxCQReverseC", "dplxCQReverseG", "dplxCQReverseT",
    "bulkProperPair", "dplxProperPair",
]


def _open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def reverse_signature(sig):
    """Reverse-complement the trinucleotide of a `NNN>indel` signature.

    Port of the perl `reverse_signature` sub (indelCaller_step1.pl:199-204):
    uppercase, take the first 3 chars, complement ACGT->TGCA, reverse, append
    '>indel'.
    """
    sig = sig.upper()
    tri = sig[:3]
    comp = tri.translate(str.maketrans("ACGT", "TGCA"))
    return comp[::-1] + ">indel"


# ===========================================================================
# propose  (= indelCaller_step1.pl)
# ===========================================================================

def _parse_dsa_line(line):
    """Split a dsa.bed.gz data line into the named columns (raw strings)."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < len(COLS):
        return None
    return dict(zip(COLS, fields))


def _totals(rec):
    """Return the per-strand duplex and bulk totals as ints."""
    di = int
    dplx_fwd = (di(rec["dplxForwardA"]) + di(rec["dplxForwardC"]) + di(rec["dplxForwardG"])
                + di(rec["dplxForwardT"]) + di(rec["dplxForwardIndel"]))
    dplx_rev = (di(rec["dplxReverseA"]) + di(rec["dplxReverseC"]) + di(rec["dplxReverseG"])
                + di(rec["dplxReverseT"]) + di(rec["dplxReverseIndel"]))
    bulk_fwd = (di(rec["bulkForwardA"]) + di(rec["bulkForwardC"]) + di(rec["bulkForwardG"])
                + di(rec["bulkForwardT"]) + di(rec["bulkForwardIndel"]))
    bulk_rev = (di(rec["bulkReverseA"]) + di(rec["bulkReverseC"]) + di(rec["bulkReverseG"])
                + di(rec["bulkReverseT"]) + di(rec["bulkReverseIndel"]))
    return dplx_fwd, dplx_rev, bulk_fwd, bulk_rev


def _passes_quality(rec, bulk_fwd, bulk_rev, args):
    """Hard quality gates shared by all propose modes (step1:110-113)."""
    if bulk_fwd + bulk_rev < args.min_coverage:
        return False
    if float(rec["dplxCLIP"]) > args.max_clip:
        return False
    if float(rec["dplxNM"]) > 20:              # hardcoded in the perl
        return False
    if float(rec["dplxASXS"]) < args.min_as_xs or float(rec["bulkASXS"]) < args.min_as_xs:
        return False
    return True


def _qpos(rec):
    """Trimmed read position and orientation (step1:143-151)."""
    chrom_beg = int(rec["chromBeg"])
    chrom_end = int(rec["chromEnd"])
    bp_beg = int(rec["dplxBreakpointBeg"])
    bp_end = int(rec["dplxBreakpointEnd"])
    qpos_f = chrom_beg - bp_beg + 1
    qpos_r = bp_end - chrom_end
    if qpos_r < qpos_f:
        return qpos_r, "Rev"
    return qpos_f, "Fwd"


def _signature_trinuc(rec, rbase):
    sig = rec["context"].replace(".", rbase, 1) + ">indel"
    if re.search(r"[AGag]", rbase):
        sig = reverse_signature(sig)
    return sig


def _classify(f_indel, r_indel, f_tot, r_tot):
    """Classify a candidate strand family for the emit-all table."""
    if f_tot == 0 and r_tot == 0:
        return "no_reads"
    if f_tot == 0 or r_tot == 0:
        return "single_strand_family"
    if f_indel == f_tot and r_indel == r_tot:
        return "ds_consensus"
    if f_indel == f_tot and r_indel == 0:
        return "ss_consensus_fwd"
    if r_indel == r_tot and f_indel == 0:
        return "ss_consensus_rev"
    if (f_indel > 0) != (r_indel > 0):
        return "nonconsensus_one_strand"
    return "nonconsensus_two_strand"


def _repeat_context(fa, chrom, pos1, window=40):
    """Heuristic reference repeat context at a 1-based position.

    Returns (homopolymer_base, homopolymer_len, str_unit, str_len). This is a
    first-cut annotation from the reference alone (the emitter cannot see the
    resolved indel allele); the modelling layer may refine it.
    """
    try:
        start0 = max(0, pos1 - 1 - window)
        seq = fa.fetch(chrom, start0, pos1 + window).upper()
    except (KeyError, ValueError):
        return ".", 0, ".", 0
    if not seq:
        return ".", 0, ".", 0
    idx = (pos1 - 1) - start0
    idx = min(max(idx, 0), len(seq) - 1)

    # Homopolymer: longest single-base run touching the site or its 3' neighbour.
    best_base, best_len = ".", 0
    for anchor in (idx, idx + 1):
        if anchor < 0 or anchor >= len(seq):
            continue
        base = seq[anchor]
        if base not in "ACGT":
            continue
        lo = anchor
        while lo - 1 >= 0 and seq[lo - 1] == base:
            lo -= 1
        hi = anchor
        while hi + 1 < len(seq) and seq[hi + 1] == base:
            hi += 1
        run = hi - lo + 1
        if run > best_len:
            best_base, best_len = base, run

    # Short tandem repeat (unit length 2..6) anchored around the site.
    best_unit, best_str = ".", 0
    for ulen in range(2, 7):
        for anchor in (idx, idx + 1):
            u0 = anchor
            if u0 < 0 or u0 + ulen > len(seq):
                continue
            unit = seq[u0:u0 + ulen]
            if any(b not in "ACGT" for b in unit):
                continue
            # extend left
            lo = u0
            while lo - ulen >= 0 and seq[lo - ulen:lo] == unit:
                lo -= ulen
            # extend right
            hi = u0 + ulen
            while hi + ulen <= len(seq) and seq[hi:hi + ulen] == unit:
                hi += ulen
            copies = (hi - lo) // ulen
            if copies >= 2 and copies * ulen > best_str:
                best_unit, best_str = unit, copies * ulen
    return best_base, best_len, best_unit, best_str


EMIT_ALL_HEADER = [
    "chrom", "pos0", "pos1", "rb_id", "dplxBarcode", "dplxOri",
    "dplxBreakpointBeg", "dplxBreakpointEnd",
    "dplxFwd_A", "dplxFwd_C", "dplxFwd_G", "dplxFwd_T", "dplxFwd_Indel", "r1",
    "dplxRev_A", "dplxRev_C", "dplxRev_G", "dplxRev_T", "dplxRev_Indel", "r2",
    "bulkFwd_A", "bulkFwd_C", "bulkFwd_G", "bulkFwd_T", "bulkFwd_Indel",
    "bulkRev_A", "bulkRev_C", "bulkRev_G", "bulkRev_T", "bulkRev_Indel",
    "indel_frac_fwd", "indel_frac_rev", "indel_frac_pooled",
    "candidate_class", "both_strands_present",
    "ref_base", "trinuc", "homopolymer_base", "homopolymer_len", "str_unit", "str_len",
    "qpos", "orientation",
    "dplxASXS", "dplxCLIP", "dplxNM", "bulkASXS", "bulkNM",
    "shearwater", "commonSNP", "bulk_seen",
]


def _emit_all_row(rec, fa):
    f_tot, r_tot, bulk_fwd, bulk_rev = _totals(rec)
    f_indel = int(rec["dplxForwardIndel"])
    r_indel = int(rec["dplxReverseIndel"])
    if f_indel + r_indel < 1:                      # candidate = >=1 indel-supporting read
        return None
    rbase = rec["context"][1]
    qpos, orient = _qpos(rec)
    cls = _classify(f_indel, r_indel, f_tot, r_tot)
    both = int(f_tot > 0 and r_tot > 0)
    frac_f = f_indel / f_tot if f_tot else 0.0
    frac_r = r_indel / r_tot if r_tot else 0.0
    frac_p = (f_indel + r_indel) / (f_tot + r_tot) if (f_tot + r_tot) else 0.0
    hp_base, hp_len, str_unit, str_len = _repeat_context(fa, rec["chrom"], int(rec["chromEnd"]))
    bulktotal = bulk_fwd + bulk_rev
    bulk_indel = int(rec["bulkForwardIndel"]) + int(rec["bulkReverseIndel"])
    bulk_seen = int(bulktotal > 0 and bulk_indel > 0)
    rb_id = "%s:%s-%s:%s" % (rec["chrom"], rec["dplxBreakpointBeg"],
                             rec["dplxBreakpointEnd"], rec["dplxBarcode"])
    row = [
        rec["chrom"], str(int(rec["chromEnd"]) - 1), rec["chromEnd"], rb_id,
        rec["dplxBarcode"], rec["dplxOri"],
        rec["dplxBreakpointBeg"], rec["dplxBreakpointEnd"],
        rec["dplxForwardA"], rec["dplxForwardC"], rec["dplxForwardG"], rec["dplxForwardT"],
        rec["dplxForwardIndel"], str(f_tot),
        rec["dplxReverseA"], rec["dplxReverseC"], rec["dplxReverseG"], rec["dplxReverseT"],
        rec["dplxReverseIndel"], str(r_tot),
        rec["bulkForwardA"], rec["bulkForwardC"], rec["bulkForwardG"], rec["bulkForwardT"],
        rec["bulkForwardIndel"],
        rec["bulkReverseA"], rec["bulkReverseC"], rec["bulkReverseG"], rec["bulkReverseT"],
        rec["bulkReverseIndel"],
        "%.6g" % frac_f, "%.6g" % frac_r, "%.6g" % frac_p,
        cls, str(both),
        rbase, rec["context"], hp_base, str(hp_len), str_unit, str(str_len),
        str(qpos), orient,
        rec["dplxASXS"], rec["dplxCLIP"], rec["dplxNM"], rec["bulkASXS"], rec["bulkNM"],
        rec["shearwater"], rec["commonSNP"], str(bulk_seen),
    ]
    return "\t".join(row)


def cmd_propose(args):
    if args.emit_all_candidates and not args.ref:
        sys.exit("--emit_all_candidates requires --ref (for repeat-context annotation)")

    fa = None
    if args.emit_all_candidates:
        import pysam
        fa = pysam.FastaFile(args.ref)

    out = _open_maybe_gzip(args.out, "wt")
    if args.emit_all_candidates:
        out.write("\t".join(EMIT_ALL_HEADER) + "\n")

    with _open_maybe_gzip(args.input, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            rec = _parse_dsa_line(line)
            if rec is None:
                continue

            f_tot, r_tot, bulk_fwd, bulk_rev = _totals(rec)
            if not _passes_quality(rec, bulk_fwd, bulk_rev, args):
                continue

            # ---- emit-all mode: rich table, no consensus/fraction/trim gates ----
            if args.emit_all_candidates:
                row = _emit_all_row(rec, fa)
                if row is not None:
                    out.write(row + "\n")
                continue

            # ---- standard modes require both strands (step1:116) ----
            if not (f_tot >= args.reads_bundle and r_tot >= args.reads_bundle):
                continue

            f_indel = int(rec["dplxForwardIndel"])
            r_indel = int(rec["dplxReverseIndel"])
            total = f_tot + r_tot

            if args.hardCutoff_single_strand_indels:
                # ssIndel fork: >=0.5 support + case1/case2 strand-consensus logic
                if f_indel + r_indel < 0.5 * total:
                    continue
                unanimous_both = (f_indel == f_tot and r_indel == r_tot)
                unanimous_one = ((f_indel == f_tot or r_indel == r_tot)
                                 and (f_indel * r_indel == 0)
                                 and (f_tot * r_tot > 0))
                if not (unanimous_both or unanimous_one):
                    continue
            else:
                # canonical double-stranded caller: >=0.9 pooled support (step1:118)
                if f_indel + r_indel < 0.9 * total:
                    continue

            qpos, _orient = _qpos(rec)
            if args.trim5 > 0 and qpos < args.trim5:
                continue
            if args.trim3 > 0 and qpos > args.trim3:
                continue

            rbase = rec["context"][1]
            dp = f_tot + r_tot
            coord1 = rec["dplxBreakpointBeg"]
            coord2 = rec["dplxBreakpointEnd"]
            bulktotal = bulk_fwd + bulk_rev
            bulk_indel = int(rec["bulkForwardIndel"]) + int(rec["bulkReverseIndel"])

            tags = "%s:%s-%s:%s;DP=%d;QPOS=%d;" % (rec["chrom"], coord1, coord2,
                                                   rec["dplxBarcode"], dp, qpos)
            tags += "%s;SW=%s;cSNP=%s" % (_signature_trinuc(rec, rbase),
                                          rec["shearwater"], rec["commonSNP"])
            tags += ";BBEG=%s" % rec["dplxBreakpointBeg"]
            tags += ";BEND=%s" % rec["dplxBreakpointEnd"]
            tags += ";DEPTH_FWD=%d" % f_tot
            tags += ";DEPTH_REV=%d" % r_tot
            tags += ";DEPTH_NORM_FWD=%d" % bulk_fwd
            tags += ";DEPTH_NORM_REV=%d" % bulk_rev
            tags += ";DPLX_ASXS=%s" % rec["dplxASXS"]
            tags += ";DPLX_CLIP=%s" % rec["dplxCLIP"]
            tags += ";DPLX_NM=%s" % rec["dplxNM"]
            tags += ";BULK_ASXS=%s" % rec["bulkASXS"]
            tags += ";BULK_NM=%s" % rec["bulkNM"]
            if bulk_indel > args.max_vaf * bulktotal:
                tags += ";BULK_SEEN(%d+%d/%d)" % (f_indel, r_indel, bulktotal)
            else:
                tags += ";"

            site = int(rec["chromEnd"])
            out.write("%s\t%d\t%d\t%s\n" % (rec["chrom"], site - 1, site, tags))

    out.close()


# ===========================================================================
# call  (= indelCaller_step2.pl)
# ===========================================================================

def _run(cmd):
    """Run a shell command, echoing it to stderr; die on nonzero exit (like runCmd)."""
    sys.stderr.write(cmd + "\n")
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        sys.exit("Error calling %s\n%s" % (cmd, proc.stderr.decode(errors="replace")))


def _cmp_chr_key(chrom):
    """Biological chromosome ordering key (port of cmp_chr, step2:92-113)."""
    c = re.sub(r"^chr", "", chrom)
    if not re.search(r"\D", c):                    # pure number
        return (0, int(c), "")
    if c == "X":
        return (1, 0, "")
    if c == "Y":
        return (1, 1, "")
    if c in ("M", "MT"):
        return (1, 2, "")
    return (2, 0, c)


# INFO/FILTER header lines injected before the first ##FORMAT line, in the exact
# order the perl produces them (step2:245-259; the perl prepends, so the code's
# last prepend appears first).
_STEP2_INJECT = [
    '##INFO=<ID=BULK_NM,Number=1,Type=Float,Description="Normal mean NM">',
    '##INFO=<ID=BULK_ASXS,Number=1,Type=Float,Description="Normal mean AS-XS">',
    '##INFO=<ID=DPLX_NM,Number=1,Type=Float,Description="RB mean NM">',
    '##INFO=<ID=DPLX_CLIP,Number=1,Type=Float,Description="RB mean CLIP">',
    '##INFO=<ID=DPLX_ASXS,Number=1,Type=Float,Description="RB mean AS-XS">',
    '##INFO=<ID=DEPTH_NORM_REV,Number=1,Type=Integer,Description="Depth normal reverse">',
    '##INFO=<ID=DEPTH_NORM_FWD,Number=1,Type=Integer,Description="Depth normal forward">',
    '##INFO=<ID=DEPTH_REV,Number=1,Type=String,Description="Depth RB reverse">',
    '##INFO=<ID=DEPTH_FWD,Number=1,Type=Integer,Description="Depth RB forward">',
    '##INFO=<ID=BEND,Number=1,Type=Integer,Description="End breakpoing">',
    '##INFO=<ID=BBEG,Number=1,Type=Integer,Description="Start breakpoint">',
    '##INFO=<ID=QPOS,Number=1,Type=Integer,Description="Read position">',
    '##INFO=<ID=RB,Number=1,Type=String,Description="Readbundle ID">',
    '##FILTER=<ID=MASKED,Description="Site overlaps with SW or SNP site">',
]


def _load_step1_bed(path):
    """Load the propose/step1 BED, returning (indels, ignore_indels_at).

    indels[rb_id][pos1] = dict of the 15 annotation fields (step2:128-169).
    """
    indels = {}
    ignore = {}
    counts = 0
    with _open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, _pos0, pos1, info = parts[0], parts[1], parts[2], parts[3]
            it = info.split(";")
            if "BULK_SEEN" in info:
                ignore.setdefault(chrom, {})[pos1] = 1
                continue
            (rb_id, dp, qpos, context, sw, snp, bbeg, bend, depth_fwd, depth_rev,
             depth_norm_fwd, depth_norm_rev, dplxASXS, dplxCLIP, dplxNM,
             bulkASXS, bulkNM) = it[0:17]
            # reformat rb_id chr:beg-end:bc1|bc2 -> chr,beg,end,bc1,bc2 (step2:143-149)
            tmp1 = re.split(r"[:\-]", rb_id)
            tmp2 = tmp1[-1].split("|")
            rb_id = ",".join([tmp1[0], tmp1[1], tmp1[2]] + tmp2)
            sw = sw.replace("SW=", "")
            snp = snp.replace("cSNP=", "")
            indels.setdefault(rb_id, {})[pos1] = {
                "dp": dp, "qpos": qpos, "context": context, "sw": sw, "snp": snp,
                "bbeg": bbeg, "bend": bend, "depth_fwd": depth_fwd, "depth_rev": depth_rev,
                "depth_norm_fwd": depth_norm_fwd, "depth_norm_rev": depth_norm_rev,
                "dplxASXS": dplxASXS, "dplxCLIP": dplxCLIP, "dplxNM": dplxNM,
                "bulkASXS": bulkASXS, "bulkNM": bulkNM,
            }
            counts += 1
    return indels, ignore, counts


def _to_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0.0


def cmd_call(args):
    if not shutil.which("samtools"):
        sys.exit("samtools not found in path")
    if not shutil.which("bcftools"):
        sys.exit("bcftools not found in path")
    if not args.ref or not os.path.exists(args.ref):
        sys.exit("Reference %s not found" % args.ref)
    if not args.bam or not os.path.exists(args.bam):
        sys.exit("BAM / CRAM %s not found" % args.bam)

    out_name = os.path.basename(args.out)
    out_dir = os.path.dirname(args.out) or "."
    sample = args.sample or "sample_1"

    indels, ignore, counts = _load_step1_bed(args.input)
    sys.stdout.write("%d indel sites seen\n" % counts)
    sys.stdout.write("%d readbundles with indels\n" % len(indels))

    tempdir = tempfile.mkdtemp(prefix="tmp.", dir=out_dir)
    tmp_bam = os.path.join(tempdir, out_name + ".tmp.bam")
    tmp_bcf = os.path.join(tempdir, out_name + ".bcf")
    tmp_vcf = os.path.join(tempdir, out_name + ".tmp.vcf")
    tmp2_vcf = os.path.join(tempdir, out_name + ".tmp2.vcf")

    header_lines = []          # emitted once
    header_done = [False]
    body_records = []          # selected record strings (one per bundle)

    try:
        for rb_id in list(indels.keys()):
            bundle = indels[rb_id]
            chrom, start_s, end_s = rb_id.split(",")[0:3]
            start = int(start_s)
            if start <= 0:
                start = 1
            end = int(end_s)

            good = [p for p in bundle if bundle[p]["sw"] == "0" and bundle[p]["snp"] == "0"]
            if not good:
                continue                            # all overlap noise / common SNP

            # ---- mini-BAM: reads tagged RB:Z:<rb_id>, duplicate flag cleared ----
            region = "%s:%d-%d" % (chrom, start, end)
            view_in = subprocess.Popen(
                ["samtools", "view", "-h", args.bam, region],
                stdout=subprocess.PIPE, text=True)
            view_out = subprocess.Popen(
                ["samtools", "view", "-bo", tmp_bam, "-"],
                stdin=subprocess.PIPE, text=True)
            rb_tag = "RB:Z:" + rb_id
            for bl in view_in.stdout:
                if bl.startswith("@"):
                    view_out.stdin.write(bl)
                elif rb_tag in bl:
                    fld = bl.split("\t")
                    flag = int(fld[1])
                    if flag > 1024:
                        fld[1] = str(flag - 1024)   # remove duplicate flag
                    view_out.stdin.write("\t".join(fld))
            view_in.stdout.close()
            view_out.stdin.close()
            view_in.wait()
            view_out.wait()

            _run("samtools index %s" % tmp_bam)
            _run("bcftools mpileup --no-BAQ --ignore-RG -L 250 -m 2 -F 0.5 "
                 "-r \"%s\" -O b -a DP,DV,DP4,SP -f %s -o %s %s"
                 % (region, args.ref, tmp_bcf, tmp_bam))
            _run("bcftools index -f %s" % tmp_bcf)
            _run("bcftools call --ploidy 1 --skip-variants snps --multiallelic-caller "
                 "--variants-only -O v %s -o %s" % (tmp_bcf, tmp_vcf))
            _run("bcftools norm -f %s %s > %s" % (args.ref, tmp_vcf, tmp2_vcf))

            get_header = not header_done[0]
            info_records = []
            best_qual = -10.0
            best_count = 0
            count = 0
            sorted_pos = sorted(bundle, key=lambda p: int(p))
            anyisok = sorted_pos[0]
            ann = bundle[anyisok]

            with open(tmp2_vcf) as vh:
                for vline in vh:
                    if vline.startswith("#"):
                        if not get_header:
                            continue
                        if not header_done[0] and vline.startswith("##FORMAT"):
                            for inj in _STEP2_INJECT:
                                header_lines.append(inj + "\n")
                            header_done[0] = True
                        if vline.startswith("#CHROM"):
                            hf = vline.rstrip("\n").split("\t")
                            hf[-1] = sample
                            vline = "\t".join(hf) + "\n"
                        elif vline.startswith("##bcftoolsCommand=mpileup"):
                            vline = ("##bcftoolsCommand=mpileup --no-BAQ --ignore-RG "
                                     "-L 250 -m 2 -F 0.5 -O b -a DP,DV,DP4,SP\n")
                        elif vline.startswith("##bcftools_callCommand=call"):
                            vline = ("##bcftools_callCommand=call --ploidy 1 "
                                     "--skip-variants snps --multiallelic-caller "
                                     "--variants-only -O v\n")
                        elif vline.startswith("##bcftools_normCommand=norm"):
                            vline = "##bcftools_normCommand=norm\n"
                        header_lines.append(vline)
                        continue

                    fields = vline.split("\t")       # last field keeps its "\n"
                    # rewrite INFO (values keep their KEY= prefixes from step1; step2:303)
                    fields[7] = ("%s;%s;RB=%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s"
                                 % (fields[7], ann["qpos"], rb_id, ann["bbeg"], ann["bend"],
                                    ann["depth_fwd"], ann["depth_rev"], ann["depth_norm_fwd"],
                                    ann["depth_norm_rev"], ann["dplxASXS"], ann["dplxCLIP"],
                                    ann["dplxNM"], ann["bulkASXS"], ann["bulkNM"]))
                    fields[6] = "PASS"
                    pos = int(fields[1])
                    ref_len = len(fields[3])
                    overlap = 0
                    for i in range(pos, pos + ref_len + 1):
                        site = bundle.get(str(i))
                        if site is not None and (site["sw"] == "1" or site["snp"] == "1"):
                            overlap += 1
                    if overlap > 0.25 * ref_len:
                        fields[6] = "MASKED"
                    # store record with FILTER PASS/MASKED (BULK_SEEN never stored; step2:317)
                    info_records.append("\t".join(fields))
                    # BULK_SEEN override affects best-selection only, not the stored string
                    filt = fields[6]
                    for i in range(pos, pos + ref_len + 1):
                        if fields[0] in ignore and str(i) in ignore[fields[0]]:
                            filt = "BULK_SEEN"
                    if filt == "PASS" and _to_float(fields[5]) > best_qual:
                        best_qual = _to_float(fields[5])
                        best_count = count
                    count += 1

            if info_records:
                body_records.append(info_records[best_count])
    finally:
        pass

    # ---- assemble, sort, finalize ----
    if not header_lines:
        header_lines = ["##fileformat=VCFv4.2\n",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample]

    unsorted_vcf = os.path.join(tempdir, out_name + ".unsorted.vcf")
    with open(unsorted_vcf, "wt") as fh:
        fh.writelines(header_lines)
        for rline in body_records:
            fh.write(rline if rline.endswith("\n") else rline + "\n")

    sorted_vcf = os.path.join(tempdir, out_name + ".vcf")
    _run("bcftools sort -Ov -T %s -o %s %s" % (tempdir, sorted_vcf, unsorted_vcf))

    final_gz = os.path.join(out_dir, out_name + ".vcf.gz")
    if args.sort:
        chrdic = {}
        with open(args.ref + ".fai") as fai:
            for fl in fai:
                p = fl.split()
                chrdic[p[0]] = p[1]
        new_head = "".join("##contig=<ID=%s,length=%s>\n" % (c, chrdic[c])
                           for c in sorted(chrdic, key=_cmp_chr_key))
        rehead_vcf = os.path.join(tempdir, out_name + ".sorted.vcf")
        wrote_new = [False]
        with open(sorted_vcf) as vin, open(rehead_vcf, "wt") as vout:
            for vl in vin:
                if vl.startswith("##contig="):
                    if not wrote_new[0]:
                        wrote_new[0] = True
                        vout.write(new_head)
                    continue
                vout.write(vl)
        _run("bcftools sort -Oz -T %s -o %s %s" % (tempdir, final_gz, rehead_vcf))
    else:
        _run("bcftools view -Oz -o %s %s" % (final_gz, sorted_vcf))

    if args.index:
        _run("bcftools index -f -t %s" % final_gz)

    if not args.keep:
        shutil.rmtree(tempdir, ignore_errors=True)
    else:
        sys.stdout.write("kept temp dir: %s\n" % tempdir)


# ===========================================================================
# verify  (= indelCaller_step3.R, ported to pysam)
# ===========================================================================

FLANK = 5
BASE_COLS_FWD = ("A", "C", "G", "T")
BASE_COLS_REV = ("a", "c", "g", "t")


def _pileup_counts(bam, ref, chrom, start, end):
    """Per-position base and indel counts over the 1-based inclusive window
    [start, end], replicating deepSNV::bam2R(q=-100, mask=3844, mq=10).

    Returns a list over positions of (n_base, n_indel).
    """
    import pysam
    af = pysam.AlignmentFile(bam, reference_filename=ref)
    npos = end - start + 1
    base = [0] * npos
    ind = [0] * npos
    # mask 3844 = UNMAP|SECONDARY|QCFAIL|DUP|SUPPLEMENTARY ; mq>=10 ; no BQ filter
    MASK = 3844
    for col in af.pileup(chrom, start - 1, end, truncate=True, stepper="all",
                         min_base_quality=0, max_depth=1_000_000):
        p = col.reference_pos + 1        # 1-based
        if p < start or p > end:
            continue
        k = p - start
        for pr in col.pileups:
            aln = pr.alignment
            if aln.flag & MASK:
                continue
            if aln.mapping_quality < 10:
                continue
            if pr.is_refskip:
                continue
            if pr.is_del:
                ind[k] += 1              # deletion covering this position ("-"/"_")
                continue
            base[k] += 1                 # a base is present here
            if pr.indel > 0:
                ind[k] += 1              # insertion follows this position ("INS"/"ins")
    af.close()
    return list(zip(base, ind))


def cmd_verify(args):
    import pysam
    genome_file = args.genome
    vcf_file = args.vcf
    bam_file = args.bam
    max_vaf = float(args.max_vaf)

    for f in (genome_file, vcf_file, bam_file):
        if not os.path.exists(f):
            sys.exit("File not found: %s" % f)

    fa = pysam.FastaFile(genome_file)

    meta, chrom_line, records = [], None, []
    with _open_maybe_gzip(vcf_file, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#CHROM"):
                chrom_line = line
            else:
                records.append(line)

    out_records = []
    for line in records:
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        # Deletion length = nchar(REF)-1 (the correct R/indelCaller_step3.R semantics).
        # NOTE: the production pipeline's bin/indelCaller_step3.R has a bug here
        # (R length() of a scalar string is always 1), so it only ever inspects a 1 bp
        # window. We keep the correct behaviour; this differs from legacy output ONLY for
        # multi-base deletions (all 1 bp indels are identical either way).
        length = max(1, len(ref) - 1)
        start = pos - FLANK
        end = pos + length + FLANK

        counts = _pileup_counts(bam_file, genome_file, chrom, start, end)
        n_bases = sum(b for b, _ in counts)
        n_indels = sum(i for _, i in counts)
        # max indel VAF over covered positions (0/0 positions skipped; see docs)
        per_site = [i / (i + b) for b, i in counts if (i + b) > 0]
        max_per_site = max(per_site) if per_site else 0.0

        if n_bases == 0:
            fields[6] = "MISSINGBULK"
        elif max_per_site > max_vaf:
            fields[6] = "NEI_IND"
        # else: keep existing FILTER (PASS)

        fields[7] = "%s;NN=[%d:%d:%s]" % (fields[7], n_indels, n_bases, repr(max_per_site))
        seq = fa.fetch(chrom, start - 3 - 1, end + 3).upper()   # scanFa(start-3 .. end+3)
        fields[7] = "%s;SEQ=%s" % (fields[7], seq)
        out_records.append("\t".join(fields) + "\n")

    # inject new FILTER/INFO meta before the first ##FILTER=<ID=MASKED line
    new_meta = []
    injected = False
    inject_block = [
        '##INFO=<ID=SEQ,Number=1,Type=String,Description="Sequence of indel plus flanking sequences">\n',
        '##INFO=<ID=NN,Number=1,Type=String,Description="n indels / n bases">\n',
        '##FILTER=<ID=NEI_IND,Description="Site was found in an indel rich region of the matched normal">\n',
        '##FILTER=<ID=MISSINGBULK,Description="Site was not found in the matched normal">\n',
    ]
    for m in meta:
        if not injected and m.startswith("##FILTER=<ID=MASKED"):
            new_meta.extend(inject_block)
            injected = True
        new_meta.append(m)
    if not injected:                              # no MASKED line present
        # Insert AFTER ##fileformat, never at position 0. The VCF spec requires
        # ##fileformat to be the first line, and prepending here produced
        # headerless-looking output that bcftools rejects with "Input is not
        # detected as bcf or vcf format". This only fires when the input has no
        # ##FILTER=<ID=MASKED line, i.e. when `call` emitted zero records -- so
        # it silently corrupts exactly the empty-result partitions, and only
        # surfaces later when `post` tries to bcftools concat them.
        insert_at = 0
        for i, m in enumerate(new_meta):
            if m.startswith("##fileformat"):
                insert_at = i + 1
                break
        new_meta = new_meta[:insert_at] + inject_block + new_meta[insert_at:]

    out_vcf = re.sub(r"\.vcf", ".filtered.vcf", vcf_file)   # e.g. .indel.filtered.vcf.gz
    tmp_plain = out_vcf[:-3] if out_vcf.endswith(".gz") else out_vcf
    with open(tmp_plain, "wt") as fh:
        fh.writelines(new_meta)
        if chrom_line:
            fh.write(chrom_line)
        fh.writelines(out_records)

    _run("bgzip -f %s" % tmp_plain)
    _run("tabix -p vcf -f %s" % out_vcf)


# ===========================================================================
# CLI
# ===========================================================================

def build_parser():
    p = argparse.ArgumentParser(
        prog="indelCaller.py",
        description="NanoSeq indel caller (Python port). Subcommands propose/call/verify "
                    "correspond to the classical Sanger steps "
                    "indelCaller_step1.pl / step2.pl / step3.R.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="subcommand", required=True)

    # ---- propose (step1) ----
    pr = sub.add_parser(
        "propose",
        help="(= classical indelCaller_step1.pl) filter dsa.bed.gz to candidate indels",
        description="propose = classical NanoSeq indelCaller_step1.pl. Filters the "
                    "per-site dsa.bed.gz table to candidate indel sites.")
    pr.add_argument("input", help="input dsa.bed.gz")
    pr.add_argument("-o", "--out", required=True, help="output file (.bed.gz, or .tsv[.gz] for emit-all)")
    pr.add_argument("-rb", "--reads-bundle", dest="reads_bundle", type=int, default=2,
                    help="minimum reads per strand in a bundle (2)")
    pr.add_argument("-t3", "--trim3", type=int, default=135,
                    help="trim reads whose qpos exceeds this (135; 0=off)")
    pr.add_argument("-t5", "--trim5", type=int, default=10,
                    help="trim reads whose qpos is below this (10; 0=off)")
    pr.add_argument("-mc", "--min-coverage", dest="min_coverage", type=int, default=20,
                    help="minimum bulk coverage (20)")
    pr.add_argument("-vaf", "--max-vaf", dest="max_vaf", type=float, default=0.2,
                    help="tag sites with at least this bulk indel VAF (0.2)")
    pr.add_argument("-a", "--min-as-xs", dest="min_as_xs", type=int, default=50,
                    help="minimum AS-XS for duplex and bulk (50)")
    pr.add_argument("-c", "--max-clip", dest="max_clip", type=float, default=0.02,
                    help="maximum clip fraction (0.02)")
    pr.add_argument("--hardCutoff_single_strand_indels", action="store_true",
                    help="emit single-stranded indels too (>=0.5 support + strand-consensus "
                         "logic; reproduces indelCaller_step1.withSingleStrandIndel.pl)")
    pr.add_argument("--emit_all_candidates", action="store_true",
                    help="emit a rich per-strand-family TSV for ALL indel candidates and stop "
                         "(no consensus/fraction/trim gates). Requires --ref.")
    pr.add_argument("--ref", help="reference FASTA (required for --emit_all_candidates)")
    pr.set_defaults(func=cmd_propose)

    # ---- call (step2) ----
    ca = sub.add_parser(
        "call",
        help="(= classical indelCaller_step2.pl) call indels per read bundle",
        description="call = classical NanoSeq indelCaller_step2.pl. Groups candidates by "
                    "read bundle and calls indels with samtools/bcftools.")
    ca.add_argument("input", help="input BED.gz from propose")
    ca.add_argument("-o", "--out", required=True, help="output prefix")
    ca.add_argument("-r", "--ref", required=True, help="reference FASTA")
    ca.add_argument("-b", "--bam", required=True, help="duplex NanoSeq BAM/CRAM")
    ca.add_argument("-s", "--sample", default="sample_1", help="sample name (sample_1)")
    ca.add_argument("-k", "--keep", action="store_true", help="keep intermediate files")
    ca.add_argument("-t", "--sort", action="store_true", help="reorder contigs biologically")
    ca.add_argument("-i", "--index", action="store_true", help="tabix-index the output")
    ca.set_defaults(func=cmd_call)

    # ---- verify (step3) ----
    ve = sub.add_parser(
        "verify",
        help="(= classical indelCaller_step3.R) check calls against the matched normal",
        description="verify = classical NanoSeq indelCaller_step3.R (ported to pysam). "
                    "Flags indel-rich (NEI_IND) or missing (MISSINGBULK) sites in the "
                    "matched normal.")
    ve.add_argument("genome", help="reference FASTA")
    ve.add_argument("vcf", help="VCF from call (step2)")
    ve.add_argument("bam", help="matched-normal BAM/CRAM")
    ve.add_argument("max_vaf", help="maximum neighbouring indel VAF")
    ve.set_defaults(func=cmd_verify)

    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
