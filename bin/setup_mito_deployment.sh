#!/usr/bin/env bash
# Populate a deployment folder for a mitochondrial NanoSeq run.
#
#   setup_mito_deployment.sh -d <deploy_dir> -i <input.tsv> [-r <repo>] [-g hg38|hg19] [--copy-code]
#
# Two things this handles that are easy to get wrong by hand:
#
#  1. Pipeline code (bin/ perl/ R/ workflows/ config/ test/) is SYMLINKED to
#     the repo by default. Existing deployments hold private copies, and that
#     has silently served stale code twice -- once for efficiency_nanoseq.R,
#     once for indelCaller.py. Pass --copy-code to freeze a snapshot instead
#     (better for reproducibility, worse for iteration).
#
#  2. If the input sheet has a DuplexBamFile column, those BAMs are linked
#     straight into readBundle_duplex/ and ALL preprocessing is skipped. Only
#     valid if the BAMs already carry RB tags -- the script verifies this on
#     the first one and refuses if absent, because a BAM without RB tags would
#     silently produce zero read bundles rather than erroring.
#
# Required input.tsv columns: Donor, TestBamID, ControlBamID
# Optional: Fragmentation, ControlType, DuplexBamFile, ControlDuplexBamFile

set -euo pipefail

REPO="/n/data1/hms/dbmi/park/vinay/pipelines/external/NanoSeq_parkLab"
GENOME="hg38"
COPY_CODE=0
DEPLOY=""
SHEET=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -d|--deploy) DEPLOY="$2"; shift 2 ;;
    -i|--input)  SHEET="$2";  shift 2 ;;
    -r|--repo)   REPO="$2";   shift 2 ;;
    -g|--genome) GENOME="$2"; shift 2 ;;
    --copy-code) COPY_CODE=1; shift ;;
    -h|--help)   sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$DEPLOY" && -n "$SHEET" ]] || { echo "need -d <deploy_dir> and -i <input.tsv>" >&2; exit 2; }
[[ -f "$SHEET" ]] || { echo "input sheet not found: $SHEET" >&2; exit 2; }

case "$GENOME" in
  hg38) CFG="config/grch38.yaml"; INTERVALS="/n/data1/hms/dbmi/park/vinay/referenceGenomes/intervals_for_varcall/hg38_smaht/" ;;
  hg19) CFG="config/grch37.yaml"; INTERVALS="/n/data1/hms/dbmi/park/vinay/referenceGenomes/intervals_for_varcall/hg19/" ;;
  *) echo "genome must be hg38 or hg19" >&2; exit 2 ;;
esac

echo "deploy : $DEPLOY"
echo "sheet  : $SHEET"
echo "repo   : $REPO"
echo "genome : $GENOME ($CFG)"
mkdir -p "$DEPLOY"
cd "$DEPLOY"

# ---- pipeline code ---------------------------------------------------------
for d in bin perl R workflows config test; do
  rm -rf "$d"
  if [[ $COPY_CODE -eq 1 ]]; then cp -r "$REPO/$d" "$d"; else ln -s "$REPO/$d" "$d"; fi
done
for f in Snakefile.table_input.mito_nanoseq.smk run_nanoseq_mito_analysis_from_bams.sh; do
  cp -f "$REPO/$f" .
done
chmod +x run_nanoseq_mito_analysis_from_bams.sh 2>/dev/null || true
echo "code   : $([[ $COPY_CODE -eq 1 ]] && echo copied || echo symlinked)"

# The mito workflow does not use the genome-wide interval files, but the parent
# Snakefile globs intervals/ at parse time and derives job counts from it.
[[ -e intervals ]] || ln -s "$INTERVALS" intervals

mkdir -p input_bam input_control_bam readBundle_duplex control_duplex runlogs logs mito_mask

# ---- input sheet -----------------------------------------------------------
[[ -e input.tsv ]] && rm -f input.tsv
ln -s "$(readlink -f "$SHEET")" input.tsv

python3 - "$SHEET" <<'PY'
import csv, os, subprocess, sys
sheet = sys.argv[1]
with open(sheet) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
if not rows:
    sys.exit("input sheet is empty")
cols = rows[0].keys()
for req in ("Donor", "TestBamID", "ControlBamID"):
    if req not in cols:
        sys.exit(f"input sheet missing required column: {req}")

def link(src, dst):
    if not src:
        return "blank"
    src = os.path.expanduser(src)
    if not os.path.exists(src):
        return "MISSING"
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.realpath(src), dst)
    for ext in (".bai", ".crai"):
        if os.path.exists(src + ext):
            if os.path.lexists(dst + ext):
                os.remove(dst + ext)
            os.symlink(os.path.realpath(src + ext), dst + ext)
    return "ok"

def has_rb(bam, scan=2000):
    # Scan the first `scan` reads, not just one: bamaddreadbundles assigns the
    # uppercase RB:Z: bundle id to proper pairs only, so supplementary and
    # secondary alignments carry just the lowercase raw barcodes (rb:Z:/mb:Z:).
    # Those cluster at the start of chr1, where ~30% of leading reads have no
    # RB -- checking read #1 alone rejected a perfectly good BAM.
    #
    # Read incrementally and stop early. capture_output=True here would buffer
    # the entire BAM -- on a multi-GB duplex BAM that is an instant OOM kill.
    try:
        p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True)
        found = False
        for _ in range(scan):
            line = p.stdout.readline()
            if not line:
                break
            if "RB:Z:" in line:
                found = True
                break
        p.stdout.close()
        p.terminate()
        p.wait(timeout=10)
        return found
    except Exception:
        return None

dup_col = "DuplexBamFile" if "DuplexBamFile" in cols else None
ctl_col = "ControlDuplexBamFile" if "ControlDuplexBamFile" in cols else None
stats = {"ok": 0, "MISSING": 0, "blank": 0}
checked = False

for r in rows:
    if dup_col:
        src = (r.get(dup_col) or "").strip()
        if src and not checked:
            rb = has_rb(src)
            if rb is False:
                sys.exit(f"\nERROR: {src}\n  has no RB:Z tag on its first read. Linking it into\n"
                         f"  readBundle_duplex/ would yield zero read bundles silently.\n"
                         f"  Drop the DuplexBamFile column and let the workflow do read bundling,\n"
                         f"  or point it at bamaddreadbundles output.")
            if rb is None:
                print("  note: could not verify RB tags (samtools unavailable?)")
            checked = True
        stats[link(src, f"readBundle_duplex/{r['TestBamID']}.filtered.bam")] += 1
    if ctl_col:
        link((r.get(ctl_col) or "").strip(),
             f"control_duplex/{r['ControlBamID']}.diluted.ctrl.bam")
    for c, d in (("TestBamFile", "input_bam"), ("ControlBamFile", "input_control_bam")):
        if c in cols:
            k = "TestBamID" if c == "TestBamFile" else "ControlBamID"
            link((r.get(c) or "").strip(), f"{d}/{r[k]}.bam")

print(f"\nrows                : {len(rows)}")
print(f"unique TestBamID    : {len({r['TestBamID'] for r in rows})}")
print(f"unique ControlBamID : {len({r['ControlBamID'] for r in rows})}")
print(f"unique Donor        : {len({r['Donor'] for r in rows})}")
if dup_col:
    print(f"duplex BAMs linked  : {stats['ok']} ok, {stats['MISSING']} MISSING, {stats['blank']} blank")
    if stats["MISSING"]:
        print("  ^ resolve MISSING before running: those samples cannot be processed")
PY

echo ""
echo "next:"
echo "  cd $DEPLOY"
echo "  snakemake -n --rerun-triggers mtime -s Snakefile.table_input.mito_nanoseq.smk --configfile $CFG \\"
echo "    2>&1 | sed -n '/^Job stats:/,/^total/p'"
