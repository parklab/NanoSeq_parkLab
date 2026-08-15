#!/usr/bin/env Rscript

# Circos-style coverage plots over the mitochondrial contig.
#
# Usage:
#   plot_mito_circos.R --mode sample   --tracks f1.tsv.gz[,f2...] --out out.pdf [--bin 10] [--numt-bed BED]
#   plot_mito_circos.R --mode combined --tracks f1.tsv.gz,f2...   --out out.pdf [--bin 10] [--numt-bed BED]
#                      [--summary summary.tsv]
#
# Per-sample depth rings are linear and autoscaled per ring: within one sample
# there is no cross-sample range to reconcile, and a log axis would compress a
# 2,000-6,000x band into a near-flat ring. The combined figure uses log10(x + 1)
# instead, because there samples genuinely span orders of magnitude; its axis
# labels are back-transformed to linear depth so they stay readable.
#
# --numt-bed is an optional hook: when a BED of NUMT segments in mito
# coordinates is supplied, it is drawn as an extra ring marking regions that
# should be masked or treated with caution. Absent it, the low-MAPQ ring is the
# data-driven stand-in for the same concern.

suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
})

# ---- args ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
mode      <- getarg("--mode", "sample")
track_arg <- getarg("--tracks")
out_pdf   <- getarg("--out")
bin_size  <- as.integer(getarg("--bin", "10"))
numt_bed  <- getarg("--numt-bed", NA)
summ_out  <- getarg("--summary", NA)

if (is.na(track_arg) || is.na(out_pdf)) {
  stop("must supply --tracks and --out", call. = FALSE)
}
track_files <- strsplit(track_arg, ",", fixed = TRUE)[[1]]
track_files <- track_files[nzchar(track_files)]

# ---- load ------------------------------------------------------------------
dt <- rbindlist(lapply(track_files, fread))
if (!nrow(dt)) stop("no rows loaded from --tracks", call. = FALSE)
setorder(dt, sample, pos)
CONTIG_LEN <- max(dt$pos)
samples <- unique(dt$sample)

# Bin for plotting only; the per-base TSV remains the data product. Binning
# keeps the PDF light without washing out fragmentation periodicity at 10 bp.
dt[, bin := ((pos - 1) %/% bin_size) * bin_size]
binned <- dt[, .(read_depth_raw = mean(read_depth_raw),
                 rb_depth_all   = mean(rb_depth_all),
                 rb_depth_a4s2  = mean(rb_depth_a4s2),
                 duplex_frac    = mean(duplex_frac, na.rm = TRUE),
                 lowmq_frac     = mean(lowmq_frac,  na.rm = TRUE)),
             by = .(sample, fragmentation, start = bin)]
binned[, end := pmin(start + bin_size, CONTIG_LEN)]
binned[is.nan(duplex_frac), duplex_frac := NA_real_]
binned[is.nan(lowmq_frac),  lowmq_frac  := NA_real_]

numt <- NULL
if (!is.na(numt_bed) && nzchar(numt_bed) && file.exists(numt_bed)) {
  numt <- fread(numt_bed, header = FALSE, select = 1:3,
                col.names = c("chr", "start", "end"))
  numt <- numt[, .(chr = "chrM", start, end)]
}

# ---- mitochondrial ideogram ------------------------------------------------
# The mito genome has no cytobands, so its gene map is the meaningful
# equivalent: 13 protein-coding genes, 2 rRNAs, 22 tRNAs and the D-loop.
# Only valid against rCRS (16,569 bp); silently skipped for any other contig
# length so a non-rCRS reference cannot be mislabelled.
RCRS_LEN <- 16569
ideo <- NULL
ideo_file <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                        value = TRUE)[1])),
                       "rCRS_mito_ideogram.tsv")
if (!file.exists(ideo_file)) ideo_file <- "workflows/rCRS_mito_ideogram.tsv"
if (CONTIG_LEN == RCRS_LEN && file.exists(ideo_file)) {
  # skip="start" makes fread seek the header line, stepping over the # preamble
  ideo <- fread(ideo_file, sep = "\t", header = TRUE, skip = "start")
  ideo[, `:=`(start = start - 1L)]           # to BED-style half-open
} else if (CONTIG_LEN != RCRS_LEN) {
  cat(sprintf("note: contig is %d bp, not rCRS (%d) -- omitting ideogram ring\n",
              CONTIG_LEN, RCRS_LEN))
}

IDEO_COL <- c("D-loop" = "#F2C14E", "rRNA" = "#D96C6C",
              "protein" = "#5B8FF9", "tRNA" = "#9AA5B1")

lg <- function(x) log10(pmax(x, 0) + 1)

# Axis labels back-transformed from log10 space to linear depth.
log_axis <- function(ymax) {
  ticks <- c(0, 1, 10, 100, 1000, 10000, 1e5)
  ticks <- ticks[lg(ticks) <= ymax]
  circos.yaxis(side = "left", at = lg(ticks), labels = format(ticks, big.mark = ",",
               trim = TRUE, scientific = FALSE), labels.cex = 0.3, tick.length = 0.1)
}

init_circos <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.after = 10, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.004, 0.004), points.overflow.warning = FALSE)
  circos.genomicInitialize(data.frame(chr = "chrM", start = 0, end = CONTIG_LEN),
                           plotType = c("axis", "labels"), axis.labels.cex = 0.4,
                           labels.cex = 0.7, major.by = 2000)
}

numt_ring <- function() {
  if (is.null(numt) || !nrow(numt)) return(invisible(NULL))
  circos.genomicTrack(as.data.frame(numt[, .(chr, start, end, v = 1)]),
    ylim = c(0, 1), track.height = 0.03, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                         col = "#B23A48", border = NA)
    })
}

# Gene-map ideogram. tRNAs are drawn as unlabelled ticks and the larger genes
# are labelled inside the ring; at 16.5 kb a label per feature would collide.
ideogram_ring <- function(label = TRUE, height = 0.055, min_label_bp = 600) {
  if (is.null(ideo) || !nrow(ideo)) return(invisible(NULL))
  circos.genomicTrack(
    as.data.frame(ideo[, .(chr = "chrM", start, end, type)]),
    ylim = c(0, 1), track.height = height, bg.border = NA,
    panel.fun = function(region, value, ...) {
      tp <- as.character(value[[1]])
      is_trna <- tp == "tRNA"
      # tRNAs sit as short centred ticks so they read as punctuation.
      circos.genomicRect(region, value,
                         ytop = ifelse(is_trna, 0.72, 1),
                         ybottom = ifelse(is_trna, 0.28, 0),
                         col = IDEO_COL[tp], border = "white", lwd = 0.25)
      if (label) {
        w <- region[[2]] - region[[1]]
        sel <- which(!is_trna & w >= min_label_bp)
        if (length(sel)) {
          nm <- as.character(ideo$name[match(region[[1]][sel], ideo$start)])
          circos.text((region[[1]][sel] + region[[2]][sel]) / 2, 0.5, nm,
                      cex = 0.42, col = "white", font = 2,
                      facing = "clockwise", niceFacing = TRUE)
        }
      }
    })
}

ideogram_legend <- function(cex = 0.75) {
  if (is.null(ideo) || !nrow(ideo)) return(invisible(NULL))
  legend("bottomright", bty = "n", cex = cex, title = "rCRS features",
         title.adj = 0, legend = names(IDEO_COL), fill = IDEO_COL,
         border = "white", text.col = "grey20")
}

# ---- per-sample figure -----------------------------------------------------
plot_sample <- function(s) {
  d <- binned[sample == s]
  frag <- d$fragmentation[1]
  gp <- function(col) as.data.frame(d[, .(chr = "chrM", start, end, v = get(col))])

  init_circos()
  ideogram_ring()
  numt_ring()

  # Linear, per-ring autoscale. Each sample is its own panel, so there is no
  # cross-sample range to reconcile here, and a log axis would compress a
  # 2,000-6,000x band into a near-flat ring and hide exactly the structure
  # these plots exist to show. The combined figure keeps log10 because there
  # samples genuinely span orders of magnitude.
  depth_ring <- function(col, fill, border) {
    vals <- d[[col]]
    ymax <- max(vals, na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    circos.genomicTrack(gp(col), ylim = c(0, ymax), track.height = 0.135,
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value[[1]], area = TRUE,
                            col = fill, border = border, lwd = 0.3)
        circos.yaxis(side = "left", at = c(0, ymax),
                     labels = format(round(c(0, ymax)), big.mark = ",",
                                     trim = TRUE, scientific = FALSE),
                     labels.cex = 0.3, tick.length = 0.08)
      })
  }
  depth_ring("read_depth_raw", "#C3CBD4", "#8794A1")
  depth_ring("rb_depth_all",   "#9DB8F2", "#3B5BC0")
  depth_ring("rb_depth_a4s2",  "#8FD3A4", "#2F9E44")

  # Fractions are drawn as points, not filled areas, and only where the
  # denominator supports them. Two reasons: an area fill has to invent a value
  # across uncovered bases (a restriction-digested library leaves ~60% of chrM
  # at zero coverage, and the resulting polygon floods the whole ring), and a
  # ratio over a denominator of 1-4 is 0/1 noise rather than a measurement.
  # Gaps here mean "not enough coverage to say", which is the honest reading.
  frac_ring <- function(col, denom_col, min_denom, pt_col, height) {
    keep <- d[!is.na(get(col)) & get(denom_col) >= min_denom]
    circos.genomicTrack(gp(col), ylim = c(0, 1), track.height = height,
      panel.fun = function(region, value, ...) {
        circos.lines(c(0, CONTIG_LEN), c(0.5, 0.5), col = "grey88", lwd = 0.4)
        if (nrow(keep)) {
          circos.genomicPoints(
            as.data.frame(keep[, .(chr = "chrM", start, end)]),
            as.data.frame(keep[, .(v = get(col))]),
            pch = 16, cex = 0.22, col = pt_col)
        }
        circos.yaxis(side = "left", at = c(0, 0.5, 1), labels = c("0", ".5", "1"),
                     labels.cex = 0.45, tick.length = 0.08)
      })
  }
  frac_ring("duplex_frac", "rb_depth_all",   5,  "#5F3DC4", 0.11)
  frac_ring("lowmq_frac",  "read_depth_raw", 20, "#B23A48", 0.10)

  ds <- dt[sample == s]
  mean_a4s2 <- mean(ds$rb_depth_a4s2)
  breadth <- 100 * mean(ds$rb_depth_a4s2 > 0)

  # Title at the top of the page rather than in the middle of the circle: the
  # centre belongs to the innermost track, and long TestBamIDs collided with it.
  title(main = s, line = 1.4, cex.main = 1.05, font.main = 2)
  mtext(sprintf("%s   |   chrM %s bp   |   a4s2 bundle depth %.0fx, breadth %.1f%%",
                frag, format(CONTIG_LEN, big.mark = ","), mean_a4s2, breadth),
        side = 3, line = 0.2, cex = 0.72, col = "grey30")

  legend("bottomleft", bty = "n", cex = 0.75, y.intersp = 1.15,
         legend = c("raw read depth", "all read bundles", "a4s2 duplex bundles",
                    "duplex fraction", "low MAPQ fraction"),
         fill = c("#C3CBD4", "#9DB8F2", "#8FD3A4", "#5F3DC4", "#B23A48"),
         border = c("#8794A1", "#3B5BC0", "#2F9E44", "#5F3DC4", "#B23A48"),
         text.col = "grey20")
  ideogram_legend(cex = 0.75)
  circos.clear()
}

# ---- combined figure -------------------------------------------------------
plot_combined <- function() {
  wide <- function(col) {
    w <- dcast(binned, start + end ~ sample, value.var = col)
    cbind(data.frame(chr = "chrM"), as.data.frame(w))
  }
  frag_of <- binned[, .(frag = fragmentation[1]), by = sample]
  setkey(frag_of, sample)
  pal <- c(RENS = "#E8590C", WGNS = "#1971C2")
  col_for <- function(nms) {
    f <- frag_of[nms, frag]
    ifelse(is.na(pal[f]), "grey50", pal[f])
  }

  init_circos()
  ideogram_ring()
  numt_ring()

  multi_ring <- function(col, logscale, height) {
    w <- wide(col)
    vcols <- setdiff(names(w), c("chr", "start", "end"))
    vals <- as.matrix(w[, vcols, drop = FALSE])
    ymax <- if (logscale) lg(max(vals, na.rm = TRUE)) else 1
    cols <- col_for(vcols)
    circos.genomicTrack(w, ylim = c(0, ymax), track.height = height,
      panel.fun = function(region, value, ...) {
        for (i in seq_along(vcols)) {
          v <- value[[i]]
          circos.genomicLines(region, if (logscale) lg(v) else v,
                              col = cols[i], lwd = 0.6)
        }
        if (logscale) log_axis(ymax) else
          circos.yaxis(side = "left", at = c(0, 0.5, 1), labels = c("0", ".5", "1"),
                       labels.cex = 0.3, tick.length = 0.1)
      })
  }
  multi_ring("read_depth_raw", TRUE,  0.17)
  multi_ring("rb_depth_all",   TRUE,  0.17)
  multi_ring("rb_depth_a4s2",  TRUE,  0.17)
  multi_ring("duplex_frac",    FALSE, 0.12)

  title(main = sprintf("chrM coverage, %d samples", length(samples)),
        line = 1.4, cex.main = 1.05, font.main = 2)
  mtext("rings, outer to inner: raw reads / all bundles / a4s2 bundles / duplex fraction",
        side = 3, line = 0.2, cex = 0.72, col = "grey30")
  legend("bottomleft", bty = "n", cex = 0.8, lwd = 2.5, y.intersp = 1.15,
         legend = names(pal), col = pal, text.col = "grey20")
  ideogram_legend(cex = 0.72)
  circos.clear()
}

# ---- summary ---------------------------------------------------------------
if (!is.na(summ_out)) {
  summ <- dt[, {
    med <- median(read_depth_raw)
    .(fragmentation      = fragmentation[1],
      mean_read_depth    = mean(read_depth_raw),
      mean_rb_depth_all  = mean(rb_depth_all),
      mean_rb_depth_a4s2 = mean(rb_depth_a4s2),
      breadth_a4s2_pct   = 100 * mean(rb_depth_a4s2 > 0),
      median_duplex_frac = median(duplex_frac, na.rm = TRUE),
      mean_lowmq_frac    = mean(lowmq_frac, na.rm = TRUE),
      # Evenness: restriction-digested libraries give spikier coverage than
      # nuclease/sonication, so CV and the within-2-fold fraction act as a
      # fragmentation fingerprint independent of the metadata label.
      cv_read_depth      = sd(read_depth_raw) / mean(read_depth_raw),
      frac_within_2x_med = mean(read_depth_raw >= med / 2 & read_depth_raw <= med * 2))
  }, by = sample]
  setorder(summ, -mean_rb_depth_a4s2)
  fwrite(summ, summ_out, sep = "\t")
}

# ---- draw ------------------------------------------------------------------
pdf(out_pdf, width = 7.8, height = 8.1)
par(mar = c(1, 1, 3.4, 1))
if (mode == "combined") {
  plot_combined()
} else {
  for (s in samples) plot_sample(s)
}
invisible(dev.off())
cat("wrote", out_pdf, "\n")
