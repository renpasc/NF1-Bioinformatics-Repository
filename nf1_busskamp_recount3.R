# NF1 exon 31 (legacy 23a) PSI and splicing-regulator expression in a public
# NGN induction time course that starts from undifferentiated stem cells
# (Busskamp et al., 2014, GEO GSE60548), using recount3.
#
# Answers: does the isoform II -> isoform I switch happen in the first days of
# NGN-driven induction, and when do REST, CELF3-6 and ELAVL2-4 change?
#
# Run in a fresh R session:   source(file.choose(), echo = TRUE)
# Needs internet access to NCBI GEO and the recount3 servers (duffel.rail.bio).
# Outputs are written to the R working directory (check with getwd()).
#
# Method matches the thesis GTEx analyses:
#   PSI = mean(left inclusion, right inclusion) / (mean(left, right) + skipping)
#   TPM from recount3 gene counts (transform_counts), denominators over all genes.

GSE        <- "GSE60548"   # Busskamp et al. (2014)
PROJECT    <- NA           # set to the SRP accession by hand if automatic lookup fails
MIN_READS  <- 10           # minimum informative junction reads per sample (match Section 2.8)
EXPECT_EXON <- c(31252938, 31253000)   # exon 31, GRCh38
EXON_LEN   <- 63
GENES      <- c("NF1", "REST", "CELF3", "CELF4", "CELF5", "CELF6",
                "ELAVL2", "ELAVL3", "ELAVL4")
DAYS_FILE  <- "busskamp_sample_days.csv"   # optional override: columns run, day, group

# 0. Packages (installed automatically on a new computer) -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("recount3", "SummarizedExperiment", "GenomicRanges", "IRanges", "Matrix")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})
options(timeout = max(600, getOption("timeout")))

# 1. Find the SRA project for the GEO series ---------------------------------------------
if (is.na(PROJECT)) {
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSE,
                "&targ=self&form=text&view=brief")
  txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) character())
  srp <- unique(regmatches(txt, regexpr("SRP[0-9]+", txt)))
  if (!length(srp)) stop("Couldn't read the SRA project from GEO. Open ", GSE,
                         " on GEO, copy the SRP number under 'Relations', and set PROJECT.")
  PROJECT <- srp[1]
}
message("SRA project: ", PROJECT)

hp <- available_projects()
pinfo <- hp[hp$project == PROJECT & hp$project_type == "data_sources", ]
if (!nrow(pinfo)) stop(PROJECT, " is not in recount3. Use the nf-core route instead: ",
                       "download the FASTQs, run nf-core/rnaseq, then nf1_exon31_psi.py.")
print(pinfo)

# 2. Sample metadata and differentiation day ---------------------------------------------
# Defined as functions so the parsing can be checked without downloading anything.
parse_day <- function(x) {
  x <- tolower(x)
  pat <- "(?<![a-z])(?:day|d)[ _.:-]*([0-9]+(?:\\.[0-9]+)?)"
  has <- grepl(pat, x, perl = TRUE)
  d <- rep(NA_real_, length(x))
  m <- regmatches(x, regexpr(pat, x, perl = TRUE))
  d[has] <- as.numeric(sub(pat, "\\1", m, perl = TRUE))
  pat2 <- "([0-9]+(?:\\.[0-9]+)?) ?(?:days?|d)(?![a-z])"   # e.g. "2 days", "4d"
  has2 <- is.na(d) & grepl(pat2, x, perl = TRUE)
  m2 <- regmatches(x[has2], regexpr(pat2, x[has2], perl = TRUE))
  d[has2] <- as.numeric(sub(pat2, "\\1", m2, perl = TRUE))
  d[is.na(d) & grepl("undiff|uninduced|no dox|-dox|ipsc only", x)] <- 0
  d
}

sample_table <- function(md) {
  cols <- intersect(c("sra.sample_title", "sra.experiment_title", "sra.sample_attributes"),
                    names(md))
  label <- apply(md[, cols, drop = FALSE], 1, function(r) paste(r[!is.na(r)], collapse = " | "))
  data.frame(run = rownames(md), label = label, day = parse_day(label),
             group = "", stringsAsFactors = FALSE)
}

# 3. PSI from junction counts --------------------------------------------------------------
# j: data.frame(idx, start, end) of intron coordinates; tot: total reads per junction.
find_cassette <- function(j, tot, exon_len = EXON_LEN) {
  L  <- merge(j, j, by = "start", suffixes = c(".skip", ".left"));  L <- L[L$end.left < L$end.skip, ]
  R  <- merge(j, j, by = "end",   suffixes = c(".skip", ".right")); R <- R[R$start.right > R$start.skip, ]
  LR <- merge(L, R, by = "idx.skip")
  LR <- LR[LR$start.right - LR$end.left - 1 == exon_len, ]
  if (!nrow(LR)) return(NULL)
  LR$exon_start <- LR$end.left + 1
  LR$exon_end   <- LR$start.right - 1
  LR$cov_left   <- tot[LR$idx.left]; LR$cov_right <- tot[LR$idx.right]; LR$cov_skip <- tot[LR$idx.skip]
  LR$support    <- pmin(LR$cov_left, LR$cov_right, LR$cov_skip)
  LR[order(-LR$support), ]
}

psi_table <- function(cnt3, runs) {
  # cnt3: 3 x samples matrix (left inclusion, right inclusion, skipping)
  incl <- (cnt3[1, ] + cnt3[2, ]) / 2
  n    <- incl + cnt3[3, ]
  data.frame(run = runs, left = cnt3[1, ], right = cnt3[2, ], skip = cnt3[3, ],
             informative = n, psi_pct = ifelse(n > 0, 100 * incl / n, NA),
             stringsAsFactors = FALSE)
}

wilson <- function(x, n, z = 1.96) {
  p <- x / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}

summarise_by_day <- function(res, min_reads = MIN_READS) {
  ok <- res[res$informative >= min_reads & !is.na(res$day), ]
  if (!nrow(ok)) return(NULL)
  keys <- unique(ok[, c("group", "day")])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    g <- ok[ok$group == keys$group[k] & ok$day == keys$day[k], ]
    i <- sum((g$left + g$right) / 2); n <- sum(g$informative)
    ci <- wilson(i, n)
    data.frame(group = keys$group[k], day = keys$day[k], n_samples = nrow(g),
               mean_psi_pct = round(mean(g$psi_pct), 2),
               sd_psi_pct = if (nrow(g) > 1) round(sd(g$psi_pct), 2) else NA,
               pooled_psi_pct = round(100 * i / n, 2),
               ci95_low = round(100 * ci[1], 2), ci95_high = round(100 * ci[2], 2),
               informative_reads = round(n, 1))
  }))
  out[order(out$group, out$day), ]
}

# 4. TPM from recount3 gene counts ----------------------------------------------------------
tpm_for_genes <- function(cnt, len_bp, gname, genes = GENES) {
  len_kb <- len_bp / 1000
  ok <- !is.na(len_kb) & len_kb > 0
  rpk_total <- as.numeric(crossprod(1 / len_kb[ok], as.matrix(cnt[ok, , drop = FALSE])))
  rows <- which(gname %in% genes & ok)
  found <- unique(gname[rows])
  rows <- vapply(found, function(g) {
    r <- rows[gname[rows] == g]; r[which.max(rowSums(as.matrix(cnt[r, , drop = FALSE])))]
  }, integer(1))
  tpm <- sweep(as.matrix(cnt[rows, , drop = FALSE]) / len_kb[rows], 2, rpk_total, "/") * 1e6
  rownames(tpm) <- found
  tpm
}

if (!exists("SKIP_DOWNLOAD") || !SKIP_DOWNLOAD) {

  # 5. Junctions: exon 31 PSI per sample ---------------------------------------------------
  if (!exists("rse_jx")) rse_jx <- create_rse(pinfo, type = "jxn")
  message("Junctions: ", nrow(rse_jx), "  Samples: ", ncol(rse_jx))
  md <- as.data.frame(colData(rse_jx))
  st <- sample_table(md)

  if (file.exists(DAYS_FILE)) {
    over <- read.csv(DAYS_FILE, stringsAsFactors = FALSE)
    st$day   <- over$day[match(st$run, over$run)]
    st$group <- ifelse(is.na(over$group[match(st$run, over$run)]), "",
                       over$group[match(st$run, over$run)])
    message("Using days and groups from ", DAYS_FILE)
  } else {
    write.csv(st, DAYS_FILE, row.names = FALSE)
    message("\nWrote ", DAYS_FILE, ". Check the parsed 'day' column against the labels.",
            "\nAdd a 'group' (e.g., NGN or control) if the study has more than one condition,",
            "\nsave, and re-run the script to use your corrections.")
  }
  print(st[, c("run", "day", "label")], row.names = FALSE)

  nf1 <- rse_jx[IRanges::overlapsAny(rowRanges(rse_jx),
                                     GenomicRanges::GRanges("chr17:31090000-31390000"),
                                     type = "within"), ]
  gr  <- rowRanges(nf1)
  j   <- data.frame(idx = seq_along(gr), start = start(gr), end = end(gr))
  tot <- Matrix::rowSums(assay(nf1, 1))
  LR  <- find_cassette(j, tot)
  if (is.null(LR)) stop("No 63 bp cassette exon found in NF1 for this project")
  cx  <- LR[1, ]
  message("Selected cassette exon: chr17:", cx$exon_start, "-", cx$exon_end)
  if (cx$exon_start != EXPECT_EXON[1] || cx$exon_end != EXPECT_EXON[2])
    warning("Selected exon differs from exon 31 (chr17:31,252,938-31,253,000). Check before using.")

  cnt3 <- as.matrix(assay(nf1, 1)[c(cx$idx.left, cx$idx.right, cx$idx.skip), , drop = FALSE])
  psi  <- psi_table(cnt3, colnames(nf1))
  res  <- merge(st[, c("run", "group", "day", "label")], psi, by = "run")

  # 6. Genes: regulator TPM per sample --------------------------------------------------
  if (!exists("rse_gene")) rse_gene <- create_rse(pinfo, type = "gene")
  gcnt <- transform_counts(rse_gene)
  tpm  <- tpm_for_genes(gcnt, rowData(rse_gene)$bp_length, rowData(rse_gene)$gene_name)
  tpm_long <- data.frame(run = rep(colnames(tpm), each = nrow(tpm)),
                         gene = rep(rownames(tpm), times = ncol(tpm)),
                         tpm = as.vector(tpm), stringsAsFactors = FALSE)
  tpm_wide <- as.data.frame(t(round(tpm, 2))); tpm_wide$run <- rownames(tpm_wide)
  res <- merge(res, tpm_wide, by = "run", all.x = TRUE)
  res <- res[order(res$group, res$day, res$run), ]

  # 7. Summaries ---------------------------------------------------------------------------
  by_day <- summarise_by_day(res)
  tpm_long <- merge(tpm_long, st[, c("run", "group", "day")], by = "run")
  tpm_day <- aggregate(tpm ~ group + day + gene, data = tpm_long,
                       FUN = function(x) round(mean(x), 2))
  tpm_day <- reshape(tpm_day, idvar = c("group", "day"), timevar = "gene", direction = "wide")
  names(tpm_day) <- sub("^tpm\\.", "", names(tpm_day))

  message("\nExon 31 PSI by day (samples with >= ", MIN_READS, " informative reads):")
  print(by_day, row.names = FALSE)
  message("\nMean TPM by day:")
  print(tpm_day[order(tpm_day$group, tpm_day$day), ], row.names = FALSE)

  write.csv(res,     "busskamp_nf1_samples.csv",    row.names = FALSE)
  write.csv(by_day,  "busskamp_nf1_psi_by_day.csv", row.names = FALSE)
  write.csv(tpm_day, "busskamp_regulators_tpm_by_day.csv", row.names = FALSE)
  message("Wrote busskamp_nf1_samples.csv, busskamp_nf1_psi_by_day.csv and ",
          "busskamp_regulators_tpm_by_day.csv in ", getwd())

  # 8. Plots -------------------------------------------------------------------------------
  if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(by_day)) {
    library(ggplot2)
    ok <- res[res$informative >= MIN_READS & !is.na(res$day), ]
    p1 <- ggplot(ok, aes(day, psi_pct, colour = group)) +
      geom_point(size = 2, alpha = 0.7) +
      stat_summary(fun = mean, geom = "line") +
      labs(x = "Day of NGN induction", y = "Exon 31 PSI (isoform II %)", colour = NULL) +
      theme_classic(base_size = 12)
    ggsave("busskamp_nf1_psi.png", p1, width = 5, height = 4, dpi = 600)
    tl <- tpm_long[!is.na(tpm_long$day), ]
    tl$gene <- factor(tl$gene, levels = intersect(GENES, unique(tl$gene)))
    p2 <- ggplot(tl, aes(day, log2(tpm + 1), colour = group)) +
      geom_point(size = 1.2, alpha = 0.7) +
      stat_summary(fun = mean, geom = "line") +
      facet_wrap(~ gene, scales = "free_y") +
      labs(x = "Day of NGN induction", y = expression(log[2]*"(TPM + 1)"), colour = NULL) +
      theme_classic(base_size = 11) +
      theme(strip.text = element_text(face = "italic"), strip.background = element_blank())
    ggsave("busskamp_regulators.png", p2, width = 8, height = 6, dpi = 600)
    message("Wrote busskamp_nf1_psi.png and busskamp_regulators.png")
  }
}
