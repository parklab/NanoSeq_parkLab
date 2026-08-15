#!/usr/bin/env Rscript
# Negative-binomial burden rates for mito NanoSeq. Reads the consolidated
# burden table, writes a machine-readable rates table.
#
#   mito_burden_nb.R <burden.tsv> <out.rates.tsv> [exclude1,exclude2,...]
#
# Why NB and not the Poisson CIs Sanger's post step reports: duplex burdens are
# overdispersed across libraries (donor-to-donor biology plus library-to-library
# variation on top of sampling noise). Poisson intervals assume sampling noise
# is the only source, so they are too narrow -- and cell-type conclusions here
# turn on whether intervals overlap.
#
# Three nested models per cohort, because on 2026-08-15 the unadjusted fit
# reported OL/Neuron = 0.30 (p = 0.002) that turned out to be a batch effect:
# every OL sat in one batch and 100 of 115 neurons in another, and neurons alone
# differed 1.77x between those batches. The paired within-donor estimate was
# 1.24. Always read `model` before quoting a ratio.
#
#   celltype              -- unadjusted; confounded by anything batch-structured
#   celltype + batch      -- fitted when >1 batch is present
#   celltype + donor      -- the paired contrast; only fitted when enough donors
#                            carry more than one cell type for it to be estimable
#
# A donor-adjusted fit on a cohort where donors are nearly one-per-library is
# saturated, not adjusted: theta explodes and the estimate is meaningless. Such
# fits are emitted with flag="degenerate" rather than silently reported.

suppressMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: mito_burden_nb.R <burden.tsv> <out.tsv> [exclude,...]")
inf <- args[1]; outf <- args[2]
excl <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else character(0)

d <- read.delim(inf, stringsAsFactors = FALSE)
d <- d[!(d$sample %in% excl), ]
d <- d[!is.na(d$dup_bases) & d$dup_bases > 0, ]
if (!nrow(d)) stop("no usable rows in ", inf)

res <- list()
add <- function(...) res[[length(res) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

fit <- function(sub, ycol, model, form, dep, npair) {
  sub$y <- sub[[ycol]]
  m <- try(suppressWarnings(glm.nb(form, data = sub)), silent = TRUE)
  if (inherits(m, "try-error")) return(invisible(NULL))
  rdf <- m$df.residual
  npar <- length(coef(m))
  # Degeneracy is a parameter-count problem, not a theta problem: a donor factor
  # with ~one level per library saturates the mean and the cell-type coefficient
  # stops being identified (hg38: 62 donor levels on 125 libraries).
  # A large theta on its own only says dispersion is Poisson-like, which is a
  # legitimate result -- flag it separately rather than discarding the fit.
  flag <- if (rdf < 5 || npar > nrow(sub) / 3) "degenerate"
          else if (m$theta > 1e4) "poisson_like" else "ok"
  cf <- coef(m); sm <- summary(m)$coefficients
  ci <- try(suppressWarnings(suppressMessages(confint(m))), silent = TRUE)
  if (inherits(ci, "try-error")) {
    se <- sm[, 2]; ci <- cbind(cf - 1.96 * se, cf + 1.96 * se)
  }
  base <- exp(cf[["(Intercept)"]])
  for (nm in grep("^celltype", names(cf), value = TRUE)) {
    add(deployment = dep, metric = ycol, model = model,
        celltype = sub("^celltype", "", nm), reference = "Neuron",
        rate = signif(base * exp(cf[[nm]]), 4),
        ref_rate = signif(base, 4),
        ratio = signif(exp(cf[[nm]]), 3),
        ratio_lo = signif(exp(ci[nm, 1]), 3), ratio_hi = signif(exp(ci[nm, 2]), 3),
        p = signif(sm[nm, 4], 3), theta = signif(m$theta, 4),
        n_lib = nrow(sub), n_donor = length(unique(sub$donor_id)),
        n_paired_donor = npair, resid_df = rdf, flag = flag)
  }
}

for (dep in sort(unique(d$deployment))) {
  sub <- d[d$deployment == dep, ]
  if (length(unique(sub$celltype)) < 2) next
  sub$celltype <- relevel(factor(sub$celltype), ref = "Neuron")
  # donors carrying >1 cell type: the only ones informing a within-donor contrast
  npair <- sum(tapply(sub$celltype, sub$donor_id,
                      function(x) length(unique(x))) > 1)
  nbatch <- length(unique(sub$batch))
  for (ycol in c("n_snv_mol", "n_snv_site")) {
    fit(sub, ycol, "celltype", y ~ celltype + offset(log(dup_bases)), dep, npair)
    if (nbatch > 1)
      fit(sub, ycol, "celltype+batch",
          y ~ celltype + batch + offset(log(dup_bases)), dep, npair)
    if (npair >= 2)
      fit(sub, ycol, "celltype+donor",
          y ~ celltype + donor_id + offset(log(dup_bases)), dep, npair)
  }
}

out <- if (length(res)) do.call(rbind, res) else
  data.frame(deployment = character(0))
write.table(out, outf, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d rate estimates (%d flagged degenerate)\n",
            outf, nrow(out), sum(out$flag == "degenerate", na.rm = TRUE)))
