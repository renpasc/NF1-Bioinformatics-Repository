# Expression of candidate NF1 exon 23a splicing regulators (REST, CELF3-6,
# ELAVL2-4) across NGN2 neuronal differentiation, from nf-core/rnaseq output
# (STAR-Salmon route). NF1 is included as a reference gene.
#
# Run in the R console:   source(file.choose(), echo = TRUE)
# You will be asked to choose salmon.merged.gene_tpm.tsv if it isn't found at
# the default path (<outdir>/star_salmon/salmon.merged.gene_tpm.tsv).

GENES    <- c("NF1", "REST", "CELF3", "CELF4", "CELF5", "CELF6",
              "ELAVL2", "ELAVL3", "ELAVL4")
TPM_FILE <- "star_salmon/salmon.merged.gene_tpm.tsv"
META_FILE <- "sample_metadata.csv"   # optional: columns sample, line, day
EXPRESSED_TPM <- 1

# 1. Read the merged TPM table -----------------------------------------------------------
if (!file.exists(TPM_FILE)) {
  message("Choose salmon.merged.gene_tpm.tsv")
  TPM_FILE <- file.choose()
}
tpm <- read.delim(TPM_FILE, check.names = FALSE)
if (!all(c("gene_id", "gene_name") %in% names(tpm)))
  stop("Expected gene_id and gene_name columns - is this the nf-core merged gene TPM file?")
samples <- setdiff(names(tpm), c("gene_id", "gene_name"))
message("Genes: ", nrow(tpm), "  Samples: ", length(samples))

sub <- tpm[tpm$gene_name %in% GENES, ]
if (length(setdiff(GENES, sub$gene_name)))
  warning("Not found: ", paste(setdiff(GENES, sub$gene_name), collapse = ", "))
# If a symbol maps to more than one gene ID, keep the ID with the highest total TPM
sub <- sub[order(-rowSums(sub[, samples, drop = FALSE])), ]
sub <- sub[!duplicated(sub$gene_name), ]

# 2. Sample metadata: line and differentiation day ---------------------------------------
# Uses sample_metadata.csv if present. Otherwise parses the day from names like
# PB001_1_D14 or PB10_15_day7, and takes the text before the day as the line.
if (file.exists(META_FILE)) {
  meta <- read.csv(META_FILE)
  if (!all(c("sample", "line", "day") %in% names(meta)))
    stop(META_FILE, " needs columns: sample, line, day")
} else {
  day_re <- "[_.-]?[Dd](?:ay|AY|Ay)?[_.-]?(\\d+)"
  day  <- as.numeric(sub(paste0(".*?", day_re, ".*"), "\\1", samples, perl = TRUE))
  line <- sub(paste0("^(.*?)", day_re, ".*$"), "\\1", samples, perl = TRUE)  # text before the day
  meta <- data.frame(sample = samples, line = line, day = day)
}
message("Check the parsed metadata (write ", META_FILE, " to override):")
print(meta, row.names = FALSE)
if (any(is.na(meta$day))) warning("Day missing for some samples - add them to ", META_FILE)

# 3. Long-format per-sample table --------------------------------------------------------
res <- data.frame(
  gene   = rep(sub$gene_name, times = length(samples)),
  sample = rep(samples, each = nrow(sub)),
  tpm    = as.vector(as.matrix(sub[, samples]))
)
res <- merge(res, meta, by = "sample")
res$log2_tpm <- log2(res$tpm + 1)
res$gene <- factor(res$gene, levels = intersect(GENES, sub$gene_name))

# 4. Summary per gene, line and day ------------------------------------------------------
summ <- aggregate(tpm ~ gene + line + day, data = res, FUN = function(x)
  c(n = length(x), mean = mean(x), sd = if (length(x) > 1) sd(x) else NA))
summ <- do.call(data.frame, summ)
names(summ) <- c("gene", "line", "day", "n", "mean_tpm", "sd_tpm")
summ$mean_tpm <- round(summ$mean_tpm, 2); summ$sd_tpm <- round(summ$sd_tpm, 2)
summ$expressed <- summ$mean_tpm >= EXPRESSED_TPM
summ <- summ[order(summ$gene, summ$line, summ$day), ]
print(summ, row.names = FALSE)

write.csv(res,  "nf1_regulators_neurons_samples.csv", row.names = FALSE)
write.csv(summ, "nf1_regulators_neurons_summary.csv", row.names = FALSE)
message("Wrote nf1_regulators_neurons_samples.csv and nf1_regulators_neurons_summary.csv in ", getwd())

# 5. Plot (optional, needs ggplot2) ------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pl <- ggplot(res, aes(day, log2_tpm, colour = line, group = line)) +
    geom_point(size = 1.2, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_colour_grey(start = 0, end = 0.6, name = NULL) +
    labs(x = "Day of differentiation", y = expression(log[2]*"(TPM + 1)")) +
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(face = "italic"), strip.background = element_blank(),
          legend.position = "top")
  ggsave("nf1_regulators_neurons.png", pl, width = 8, height = 6, dpi = 600)
  message("Wrote nf1_regulators_neurons.png")
} else message("Install ggplot2 to also save a plot")
