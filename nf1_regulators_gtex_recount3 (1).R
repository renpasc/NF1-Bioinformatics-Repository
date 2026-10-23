# Expression of candidate NF1 exon 23a splicing regulators (REST, CELF3-6,
# ELAVL2-4) in GTEx whole blood and EBV-transformed lymphocytes, via recount3.
# NF1 is included as a reference gene.
#
# Run in the R console:   source(file.choose(), echo = TRUE)
# The GTEx BLOOD project download is large and may take several minutes.
#
# Why recount3 rather than snapcount: exon 23a PSI needs junction counts, but
# gene expression needs gene-level counts and gene lengths to calculate TPM.
# recount3 provides both for the same GTEx samples.

needed  <- c("recount3", "SummarizedExperiment")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install first: BiocManager::install(c(\"",
                          paste(missing, collapse = "\", \""), "\"))")
suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})

GENES   <- c("NF1", "REST", "CELF3", "CELF4", "CELF5", "CELF6",
             "ELAVL2", "ELAVL3", "ELAVL4")
TISSUES <- c("Whole Blood", "Cells - EBV-transformed lymphocytes")
EXPRESSED_TPM <- 1   # threshold for calling a gene expressed in a sample

# 1. Load the GTEx BLOOD project (contains both tissues) --------------------------------
if (!exists("rse_blood")) {
  proj <- available_projects()
  p <- subset(proj, file_source == "gtex" & project == "BLOOD" &
                    project_type == "data_sources")
  if (nrow(p) != 1) stop("Couldn't find the GTEx BLOOD project in recount3")
  rse_blood <- create_rse(p, type = "gene")
}
message("Genes: ", nrow(rse_blood), "  Samples: ", ncol(rse_blood))

# 2. Find the tissue column and keep the two tissues -----------------------------------
md <- as.data.frame(colData(rse_blood))
tissue_col <- NULL
for (cn in names(md)) {
  v <- as.character(md[[cn]])
  if (any(v %in% TISSUES)) { tissue_col <- cn; break }
}
if (is.null(tissue_col)) stop("Couldn't find a column with GTEx tissue names (e.g. 'Whole Blood')")
message("Tissue column: ", tissue_col)
keep <- as.character(md[[tissue_col]]) %in% TISSUES
print(table(as.character(md[[tissue_col]])[keep]))

# 3. Read counts and TPM -----------------------------------------------------------------
# recount3 stores base-pair coverage; transform_counts() scales it to read counts.
# TPM denominators use all genes, so the full matrix is needed for that step only.
rse <- rse_blood[, keep]
cnt <- transform_counts(rse)
len_kb <- rowData(rse)$bp_length / 1000
ok <- !is.na(len_kb) & len_kb > 0
rpk_total <- as.numeric(crossprod(1 / len_kb[ok], cnt[ok, , drop = FALSE]))  # per-sample RPK sum

gname <- rowData(rse)$gene_name
rows  <- which(gname %in% GENES & ok)
found <- unique(gname[rows])
if (length(setdiff(GENES, found)))
  warning("Not found in annotation: ", paste(setdiff(GENES, found), collapse = ", "))
# If a symbol maps to more than one gene ID, keep the ID with the highest total count
rows <- vapply(found, function(g) {
  r <- rows[gname[rows] == g]; r[which.max(rowSums(cnt[r, , drop = FALSE]))]
}, integer(1))

sub_cnt <- cnt[rows, , drop = FALSE]
sub_tpm <- sweep(sub_cnt / len_kb[rows], 2, rpk_total, "/") * 1e6
rm(cnt); invisible(gc())

# 4. Long-format per-sample table --------------------------------------------------------
res <- data.frame(
  sample = rep(colnames(sub_tpm), each = length(found)),
  tissue = rep(as.character(colData(rse)[[tissue_col]]), each = length(found)),
  gene   = rep(found, times = ncol(sub_tpm)),
  count  = round(as.vector(sub_cnt)),
  tpm    = as.vector(sub_tpm)
)
res$log2_tpm <- log2(res$tpm + 1)
res$gene <- factor(res$gene, levels = intersect(GENES, found))

# 5. Summary per tissue and gene ---------------------------------------------------------
summ <- do.call(rbind, lapply(split(res, list(res$tissue, res$gene), drop = TRUE), function(g) {
  q <- quantile(g$tpm, c(.25, .5, .75))
  data.frame(tissue = g$tissue[1], gene = as.character(g$gene[1]), n_samples = nrow(g),
             median_tpm = round(q[2], 2), iqr_low = round(q[1], 2), iqr_high = round(q[3], 2),
             mean_tpm = round(mean(g$tpm), 2), sd_tpm = round(sd(g$tpm), 2),
             pct_expressed = round(100 * mean(g$tpm >= EXPRESSED_TPM), 1),
             row.names = NULL)
}))
summ <- summ[order(summ$tissue, match(summ$gene, GENES)), ]
print(summ, row.names = FALSE)

write.csv(res,  "nf1_regulators_gtex_samples.csv", row.names = FALSE)
write.csv(summ, "nf1_regulators_gtex_summary.csv", row.names = FALSE)
message("Wrote nf1_regulators_gtex_samples.csv and nf1_regulators_gtex_summary.csv in ", getwd())

# 6. Plot (optional, needs ggplot2) ------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pl <- ggplot(res, aes(gene, log2_tpm, fill = tissue)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
    scale_fill_manual(values = c("white", "grey50"), name = NULL,
                      labels = c("EBV-transformed lymphocytes", "Whole blood")) +
    labs(x = NULL, y = expression(log[2]*"(TPM + 1)")) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(face = "italic"), legend.position = "top")
  ggsave("nf1_regulators_gtex.png", pl, width = 7, height = 4.5, dpi = 600)
  message("Wrote nf1_regulators_gtex.png")
} else message("Install ggplot2 to also save a plot")
