# NF1 exon 23a inclusion (isoform II fraction) in GTEx whole blood and
# EBV-transformed lymphocytes, via Snaptron (snapcount).
#
# Run in the R console:   source(file.choose(), echo = TRUE)
# Uses the snapcount install you already have. The query takes a few minutes.

needed  <- c("snapcount", "SummarizedExperiment")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install first: BiocManager::install(c(\"",
                          paste(missing, collapse = "\", \""), "\"))")
suppressPackageStartupMessages({
  library(snapcount)
  library(SummarizedExperiment)
})

EXPECT_EXON <- c(31252938, 31253000)   # exon 23a (hg38), found in the SRA run
EXON_LEN    <- 63
MIN_READS   <- 10   # matches the thesis (Section 2.8)
TISSUES     <- c("Whole Blood", "Cells - EBV-transformed lymphocytes")

# 1. Query GTEx (separate object name so the earlier SRA `jx` isn't reused) -----------
if (!exists("jx_gtex")) {
  comps <- intersect(c("gtexv2", "gtex"), names(Compilation))
  if (!length(comps)) stop("No GTEx compilation in your snapcount. Available: ",
                           paste(names(Compilation), collapse = ", "))
  for (COMP_NAME in comps) {
    message("Trying compilation: ", COMP_NAME)
    jx_gtex <- tryCatch(
      query_jx(QueryBuilder(compilation = Compilation[[COMP_NAME]], regions = "NF1")),
      error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    if (!is.null(jx_gtex)) break
  }
  if (is.null(jx_gtex)) { rm(jx_gtex); stop("All GTEx compilations failed") }
  message("Using compilation: ", COMP_NAME)
}
message("Junctions: ", nrow(jx_gtex), "  Samples: ", ncol(jx_gtex))

# 2. Find the tissue column and keep the two tissues -----------------------------------
md <- as.data.frame(colData(jx_gtex))
tissue_col <- NULL
for (cn in names(md)) {
  v <- as.character(md[[cn]])
  if (any(v %in% TISSUES)) { tissue_col <- cn; break }
}
if (is.null(tissue_col)) stop("Couldn't find a column with GTEx tissue names (e.g. 'Whole Blood')")
message("Tissue column: ", tissue_col)
keep <- as.character(md[[tissue_col]]) %in% TISSUES
print(table(as.character(md[[tissue_col]])[keep]))

# 3. Find exon 23a (rank 63-bp cassettes by their weakest junction) ---------------------
gr  <- rowRanges(jx_gtex)
j   <- data.frame(idx = seq_along(gr), start = start(gr), end = end(gr))
tot <- Matrix::rowSums(assay(jx_gtex, 1)[, keep, drop = FALSE])

L  <- merge(j, j, by = "start", suffixes = c(".skip", ".left"));  L <- L[L$end.left < L$end.skip, ]
R  <- merge(j, j, by = "end",   suffixes = c(".skip", ".right")); R <- R[R$start.right > R$start.skip, ]
LR <- merge(L, R, by = "idx.skip")
LR <- LR[LR$start.right - LR$end.left - 1 == EXON_LEN, ]
if (!nrow(LR)) stop("No 63-bp cassette exon found in NF1")
LR$exon_start <- LR$end.left + 1
LR$exon_end   <- LR$start.right - 1
LR$cov_left   <- tot[LR$idx.left]; LR$cov_right <- tot[LR$idx.right]; LR$cov_skip <- tot[LR$idx.skip]
LR <- LR[order(-pmin(LR$cov_left, LR$cov_right, LR$cov_skip)), ]
print(head(LR[, c("exon_start", "exon_end", "cov_left", "cov_right", "cov_skip")], 5))
cx <- LR[1, ]
message("Selected cassette exon: chr17:", cx$exon_start, "-", cx$exon_end)
if (cx$exon_start != EXPECT_EXON[1] || cx$exon_end != EXPECT_EXON[2])
  warning("Selected exon differs from expected exon 23a coordinates - check before using")

# 4. Per-sample isoform II fraction ------------------------------------------------------
cnt  <- as.matrix(assay(jx_gtex, 1)[c(cx$idx.left, cx$idx.right, cx$idx.skip), keep, drop = FALSE])
incl <- (cnt[1, ] + cnt[2, ]) / 2
res  <- data.frame(sample = colnames(cnt), tissue = as.character(md[[tissue_col]])[keep],
                   left = cnt[1, ], right = cnt[2, ], skip = cnt[3, ],
                   informative = incl + cnt[3, ], isoformII_frac = incl / (incl + cnt[3, ]))

# 5. Summary per tissue: per-sample median/IQR (deep samples) + pooled with Wilson CI -----
wilson <- function(x, n) {
  z <- 1.96; p <- x / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}
summ <- do.call(rbind, lapply(split(res, res$tissue), function(g) {
  i <- sum((g$left + g$right) / 2); s <- sum(g$skip); n <- i + s; ci <- wilson(i, n)
  deep <- g$isoformII_frac[g$informative >= MIN_READS]
  q <- if (length(deep)) quantile(deep, c(.25, .5, .75)) else rep(NA, 3)
  data.frame(tissue = g$tissue[1], n_samples = nrow(g), n_deep = length(deep),
             total_reads = round(n), median_sample = round(q[2], 3),
             iqr_low = round(q[1], 3), iqr_high = round(q[3], 3),
             pooled = round(i / n, 3), ci_low = round(ci[1], 3), ci_high = round(ci[2], 3),
             row.names = NULL)
}))
print(summ, row.names = FALSE)

write.csv(res,  "nf1_exon23a_gtex_samples.csv", row.names = FALSE)
write.csv(summ, "nf1_exon23a_gtex_summary.csv", row.names = FALSE)
message("Wrote nf1_exon23a_gtex_samples.csv and nf1_exon23a_gtex_summary.csv in ", getwd())
