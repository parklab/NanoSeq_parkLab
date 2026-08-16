#!/usr/bin/env python3
"""Compile per-donor mito NanoSeq results into a master call table + burden table.

    compile_mito_master_table.py -d <deploy_dir> [-d <deploy_dir> ...] -o <outprefix>

Writes <outprefix>.calls.tsv.gz and <outprefix>.burden.tsv.

Design notes that matter for interpretation:

  RECURRENCE IS ANNOTATED, NEVER FILTERED. mtDNA is multi-copy, so a somatic
  variant that clonally expanded legitimately appears in many independent
  duplex molecules -- recurrence is real heteroplasmy, not automatically an
  artifact. What flags contamination is extreme concentration (few sites, many
  molecules each) combined with a burden outlier. Two columns support that call
  downstream: n_mol_at_pos (within sample) and n_samples_at_pos (across the
  cohort; a site called in many unrelated donors is germline or systematic).

  TWO BURDEN DEFINITIONS, because they answer different questions:
    burden_molecule = mutant molecules / duplex bases -- what NanoSeq's own
        results.mut_burden.tsv reports. Sensitive to clonal expansion.
    burden_site     = distinct mutated positions / duplex bases -- closer to a
        per-base mutation rate, insensitive to expansion.
  For nuclear NanoSeq these coincide; for mtDNA they do not, and the gap is
  itself informative.

  Nothing is excluded here. Contaminated libraries stay in the table with their
  diagnostics attached, so exclusion is a documented downstream decision.
"""
import argparse
import csv
import glob
import gzip
import os
import sys
from collections import Counter, defaultdict

# rCRS control region: 16024-16569 and 1-576 (1-based). muts.tsv chromStart is
# 0-based, so compare against 16023 / 576 directly.
DLOOP_HI, DLOOP_LO = 16023, 576

SHEET_COLS = ["Batch", "Cohort", "CellType", "BiologicalDonor", "Concentration",
              "ControlSource", "Phenotype"]

# Some sample-naming conventions embed the sorted population in the donor
# string, e.g. "<donor>NeuNpos" vs "<donor>NeuNneg", or a "<donor>NeuNp"
# shorthand. Left as-is, donor is perfectly nested within celltype, so a
# donor-adjusted model is aliased rather than paired -- it silently returns a
# cell-type estimate that is not one. donor_id strips the population token so
# one individual's NeuN+ and NeuN− libraries share an id and within-donor
# contrasts become estimable.
_CELLTYPE_TOKENS = ["NeuNpos", "NeuNneg", "NeuNp", "NeuNn", "NeuN", "Neu", "OL"]


def donor_id(biological_donor):
    s = biological_donor or ""
    for tok in _CELLTYPE_TOKENS:            # longest-first, see list order
        s = s.replace(tok, "")
    s = s.rstrip("_")
    for ploidy in ("_2N", "_4N"):
        if s.endswith(ploidy):
            s = s[: -len(ploidy)]
    return s.strip("_") or (biological_donor or "")


def load_sheet(deploy):
    hits = glob.glob(os.path.join(deploy, "input.mito*.tsv"))
    if not hits:
        return {}
    with open(hits[0]) as fh:
        return {r["TestBamID"]: r for r in csv.DictReader(fh, delimiter="\t")}


def load_summary(deploy):
    p = os.path.join(deploy, "mito_tracks", "summary.tsv")
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        return {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def bulk_vaf(r):
    tot = {b: int(r["bulkForward" + b]) + int(r["bulkReverse" + b]) for b in "ACGT"}
    d = sum(tot.values())
    if not d:
        return "", 0
    major = max(tot, key=tot.get)
    return f"{sum(v for b, v in tot.items() if b != major) / d:.6f}", d


def read_indels(path):
    if not os.path.exists(path):
        return []
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1]) - 1, f[3], f[4]))   # to 0-based
    return out


ap = argparse.ArgumentParser()
ap.add_argument("-d", "--deploy", action="append", required=True)
ap.add_argument("-o", "--out", required=True)
args = ap.parse_args()

calls, per_sample = [], {}
for deploy in args.deploy:
    tag = os.path.basename(deploy.rstrip("/"))
    sheet, summ = load_sheet(deploy), load_summary(deploy)
    for d in sorted(glob.glob(os.path.join(deploy, "*.runNanoSeqMito"))):
        s = os.path.basename(d)[:-len(".runNanoSeqMito")]
        post = os.path.join(d, "tmpNanoSeq", "post")
        mb = os.path.join(post, "results.mut_burden.tsv")
        if not os.path.exists(mb):
            print(f"  WARN no mut_burden for {s}", file=sys.stderr)
            continue
        obs = [l.split("\t") for l in open(mb).read().splitlines() if l.startswith("observed")]
        dup_bases = int(obs[0][2]) if obs else 0

        meta = sheet.get(s, {})
        base = dict(deployment=tag, sample=s,
                    **{c.lower(): meta.get(c, "") for c in SHEET_COLS})
        base["donor_id"] = donor_id(meta.get("BiologicalDonor", ""))
        sm = summ.get(s, {})
        base["a4s2_depth"] = sm.get("mean_rb_depth_a4s2", "")
        base["breadth_a4s2_pct"] = sm.get("breadth_a4s2_pct", "")

        rows = []
        mt = os.path.join(post, "results.muts.tsv")
        if os.path.exists(mt):
            with open(mt) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    vaf, bdep = bulk_vaf(r)
                    pos = int(r["chromStart"])
                    rows.append(dict(base, klass="snv", chrom=r["chrom"], pos=pos,
                                     ref_alt=r.get("stdsub", ""), context=r.get("stdcontext", ""),
                                     pyrsub=r.get("pyrsub", ""), call=r.get("call", ""),
                                     bulk_vaf=vaf, bulk_depth=bdep,
                                     dplx_fwd=r.get("dplxfwdTotal", ""), dplx_rev=r.get("dplxrevTotal", ""),
                                     dplx_barcode=r.get("dplxBarcode", ""),
                                     common_snp=r.get("commonSNP", ""), shearwater=r.get("shearwater", ""),
                                     ismasked=r.get("ismasked", "")))
        for chrom, pos, ref, alt in read_indels(os.path.join(post, "results.indel.vcf.gz")):
            rows.append(dict(base, klass="indel", chrom=chrom, pos=pos,
                             ref_alt=f"{ref}>{alt}", context="", pyrsub="", call="",
                             bulk_vaf="", bulk_depth="", dplx_fwd="", dplx_rev="",
                             dplx_barcode="", common_snp="", shearwater="", ismasked=""))

        # within-sample molecular recurrence, SNVs and indels counted separately
        c = Counter((r["klass"], r["pos"]) for r in rows)
        for r in rows:
            n = c[(r["klass"], r["pos"])]
            r["n_mol_at_pos"] = n
            r["is_recurrent"] = int(n > 1)
            r["in_dloop"] = int(r["pos"] < DLOOP_LO or r["pos"] > DLOOP_HI)
        calls.extend(rows)

        snv = [r for r in rows if r["klass"] == "snv"]
        ind = [r for r in rows if r["klass"] == "indel"]
        per_sample[(tag, s)] = dict(base, dup_bases=dup_bases,
                                    n_snv_mol=len(snv), n_snv_site=len({r["pos"] for r in snv}),
                                    n_indel_mol=len(ind), n_indel_site=len({r["pos"] for r in ind}),
                                    n_snv_recurrent=sum(r["is_recurrent"] for r in snv),
                                    n_snv_dloop=sum(r["in_dloop"] for r in snv))

# cross-sample recurrence: a site called in many unrelated donors is germline
# or systematic, and is the single most useful downstream flag.
seen = defaultdict(set)
for r in calls:
    seen[(r["klass"], r["chrom"], r["pos"])].add(r["sample"])
for r in calls:
    r["n_samples_at_pos"] = len(seen[(r["klass"], r["chrom"], r["pos"])])

CALL_COLS = ["deployment", "sample", "batch", "cohort", "celltype", "biologicaldonor",
             "donor_id", "concentration", "controlsource", "phenotype", "a4s2_depth", "breadth_a4s2_pct",
             "klass", "chrom", "pos", "ref_alt", "context", "pyrsub", "call",
             "bulk_vaf", "bulk_depth", "dplx_fwd", "dplx_rev", "dplx_barcode",
             "common_snp", "shearwater", "ismasked",
             "n_mol_at_pos", "is_recurrent", "n_samples_at_pos", "in_dloop"]
with gzip.open(args.out + ".calls.tsv.gz", "wt", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=CALL_COLS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(calls)

BURDEN_COLS = ["deployment", "sample", "batch", "cohort", "celltype", "biologicaldonor",
               "donor_id", "concentration", "controlsource", "phenotype",
               "a4s2_depth", "breadth_a4s2_pct",
               "dup_bases", "n_snv_mol", "n_snv_site", "n_indel_mol", "n_indel_site",
               "n_snv_recurrent", "n_snv_dloop",
               "burden_molecule", "burden_site", "frac_recurrent", "frac_dloop"]
with open(args.out + ".burden.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=BURDEN_COLS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for v in per_sample.values():
        db, nm = v["dup_bases"], v["n_snv_mol"]
        v["burden_molecule"] = f"{nm / db:.6e}" if db else ""
        v["burden_site"] = f"{v['n_snv_site'] / db:.6e}" if db else ""
        v["frac_recurrent"] = f"{v['n_snv_recurrent'] / nm:.4f}" if nm else ""
        v["frac_dloop"] = f"{v['n_snv_dloop'] / nm:.4f}" if nm else ""
        w.writerow(v)

print(f"{args.out}.calls.tsv.gz : {len(calls):,} calls "
      f"({sum(1 for r in calls if r['klass']=='snv'):,} SNV, "
      f"{sum(1 for r in calls if r['klass']=='indel'):,} indel)")
print(f"{args.out}.burden.tsv   : {len(per_sample)} samples")
