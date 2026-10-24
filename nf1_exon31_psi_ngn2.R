# NF1 exon 31 (legacy 23a) PSI in the control NGN2 RNA-seq, from STAR junction tables.
# Same formula as the GTEx and Busskamp analyses:
#   PSI = mean(left inclusion, right inclusion) / (mean(left, right) + skipping)
#
# Input: the nf-core/rnaseq STAR junction tables, one per sample (<sample>.SJ.out.tab).
#   nf-core/rnaseq writes them to <outdir>/star_salmon/log/. Subfolders are searched,
#   so pointing SJ_FOLDER at <outdir> is fine.
#
# Optional: sample_metadata.csv in the working directory, columns sample, line, day,
# to override the line and day parsed from file names (e.g. LINE1_D14).
#
# Run in the R console:   source(file.choose(), echo = TRUE)
# Outputs go to the working directory (getwd()).

SJ_FOLDER <- NA        # folder holding the SJ.out.tab files, or NA to choose interactively
MIN_READS <- 10        # minimum informative junction reads per sample (match Section 2.8)
EXON_START <- 31252938 # exon 31, GRCh38
EXON_END   <- 31253000
META_FILE  <- "sample_metadata.csv"
EXCLUDE    <- "NTC|no.?template"   # sample names matching this are skipped

for (p in c("ggplot2", "ragg", "systemfonts")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)

# 1. Find the SJ.out.tab files ------------------------------------------------------------
if (is.na(SJ_FOLDER) || !dir.exists(SJ_FOLDER)) {
  SJ_FOLDER <- if (.Platform$OS.type == "windows") utils::choose.dir(caption = "Choose the folder with the SJ.out.tab files") else
    dirname(file.choose())
}
message("Reading junction tables from: ", SJ_FOLDER)
files <- list.files(SJ_FOLDER, pattern = "SJ\\.out\\.tab$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No *.SJ.out.tab files under ", SJ_FOLDER,
                         ". In nf-core/rnaseq they are in star_salmon/log/.")
sample <- sub("\\.SJ\\.out\\.tab$", "", basename(files))
keep <- !grepl(EXCLUDE, sample, ignore.case = TRUE)
files <- files[keep]; sample <- sample[keep]
message("Found ", length(files), " junction tables")

# 2. Read the junctions near NF1 exon 31 --------------------------------------------------
# STAR columns: chrom, intron first base, intron last base (1-based), strand, motif,
# annotated, unique reads, multimapping reads, max overhang. Unique reads are used.
read_sj <- function(f) {
  x <- utils::read.delim(f, header = FALSE, colClasses = c("character", rep("numeric", 8)))
  x <- x[x$V1 %in% c("17", "chr17") & x$V2 > 31150000 & x$V2 < 31350000, ]
  data.frame(start = x$V2, end = x$V3, reads = x$V7)
}
sj <- lapply(files, read_sj)

# Learn the exon 30 donor and exon 32 acceptor from all samples pooled, so a sample
# with zero inclusion reads still gets PSI = 0.
pool  <- aggregate(reads ~ start + end, data = do.call(rbind, sj), FUN = sum)
left  <- pool[pool$end == EXON_START - 1, ]
right <- pool[pool$start == EXON_END + 1, ]
if (nrow(left) && nrow(right)) {
  lk <- left[which.max(left$reads), c("start", "end")]
  rk <- right[which.max(right$reads), c("start", "end")]
} else {
  span <- pool[pool$start < EXON_START & pool$end > EXON_END, ]
  if (!nrow(span)) stop("No junctions spanning NF1 exon 31. Check genome build (GRCh38) and chromosome names.")
  s <- span[which.max(span$reads), ]
  lk <- data.frame(start = s$start, end = EXON_START - 1)
  rk <- data.frame(start = EXON_END + 1, end = s$end)
}
sk <- data.frame(start = lk$start, end = rk$end)
message("Exon 31: chr17:", EXON_START, "-", EXON_END,
        " | exon 30 ends ", lk$start - 1, " | exon 32 starts ", rk$end + 1)

get <- function(x, k) { v <- x$reads[x$start == k$start & x$end == k$end]; if (length(v)) sum(v) else 0 }

# 3. Sample line and day ------------------------------------------------------------------
if (file.exists(META_FILE)) {
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
  line <- meta$line[match(sample, meta$sample)]
  day  <- meta$day[match(sample, meta$sample)]
} else {
  day_re <- "[_.-]?[Dd](?:ay|AY|Ay)?[_.-]?(\\d+)"
  day  <- suppressWarnings(as.numeric(sub(paste0(".*?", day_re, ".*"), "\\1", sample, perl = TRUE)))
  line <- sub(paste0("^(.*?)", day_re, ".*$"), "\\1", sample, perl = TRUE)
  day[!grepl(day_re, sample, perl = TRUE)] <- NA
}

# 4. PSI per sample -----------------------------------------------------------------------
res <- data.frame(sample = sample, line = line, day = day,
                  left = vapply(sj, get, numeric(1), k = lk),
                  right = vapply(sj, get, numeric(1), k = rk),
                  skip = vapply(sj, get, numeric(1), k = sk))
res$informative <- (res$left + res$right) / 2 + res$skip
res$psi_pct <- ifelse(res$informative > 0, 100 * ((res$left + res$right) / 2) / res$informative, NA)
res$passes_min_reads <- res$informative >= MIN_READS
res <- res[order(res$line, res$day, res$sample), ]
message("\nCheck the parsed line and day (write ", META_FILE, " to override):")
print(res, row.names = FALSE, digits = 3)

# 5. Summaries: per line and day, and pooled across lines per day ------------------------
wilson <- function(x, n, z = 1.96) {
  p <- x / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}
summarise <- function(g, line_label) {
  i <- sum((g$left + g$right) / 2); n <- sum(g$informative); ci <- wilson(i, n)
  data.frame(line = line_label, day = g$day[1], n_samples = nrow(g),
             mean_psi_pct = round(mean(g$psi_pct), 2),
             sd_psi_pct = if (nrow(g) > 1) round(sd(g$psi_pct), 2) else NA,
             pooled_psi_pct = round(100 * i / n, 2),
             ci95_low = round(100 * ci[1], 2), ci95_high = round(100 * ci[2], 2),
             informative_reads = n)
}
ok <- res[res$passes_min_reads & !is.na(res$day), ]
if (!nrow(ok)) stop("No samples with >= ", MIN_READS, " informative reads and a known day")
by_line <- do.call(rbind, lapply(split(ok, list(ok$line, ok$day), drop = TRUE),
                                 function(g) summarise(g, g$line[1])))
by_day  <- do.call(rbind, lapply(split(ok, ok$day), function(g) summarise(g, "All lines")))
by_line <- by_line[order(by_line$line, by_line$day), ]
message("\nExon 31 PSI by day (all lines pooled):")
print(by_day, row.names = FALSE)

write.csv(res,     "ngn2_exon31_psi_samples.csv", row.names = FALSE)
write.csv(by_line, "ngn2_exon31_psi_by_line.csv", row.names = FALSE)
write.csv(by_day,  "ngn2_exon31_psi_by_day.csv",  row.names = FALSE)

# 6. Plot in the style of Figure 3.16 (mean +/- SD per line) --------------------------------
sf <- systemfonts::system_fonts()
FONT <- if (any(sf$family == "Lato" & sf$italic)) "Lato" else ""
lines_ <- sort(unique(by_line$line))
cols   <- setNames(c("black", "grey50", "grey30", "grey70")[seq_along(lines_)], lines_)
shapes <- setNames(c(16, 22, 15, 21)[seq_along(lines_)], lines_)
fills  <- setNames(c("black", "white", "grey30", "white")[seq_along(lines_)], lines_)
by_line$sd0 <- ifelse(is.na(by_line$sd_psi_pct), 0, by_line$sd_psi_pct)
p <- ggplot(by_line, aes(day, mean_psi_pct, colour = line, shape = line, fill = line)) +
  geom_errorbar(data = by_line[by_line$sd0 > 0, ],
                aes(ymin = pmax(mean_psi_pct - sd0, 0), ymax = mean_psi_pct + sd0),
                width = 1, linewidth = 0.45) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2, stroke = 0.7) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_fill_manual(values = fills, name = NULL) +
  scale_x_continuous(breaks = sort(unique(by_line$day))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Day of differentiation", y = "Exon 31 PSI (isoform II %)") +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(axis.text = element_text(colour = "black"), legend.position = "top")
ggsave("ngn2_exon31_psi.png", p, width = 4.5, height = 4, dpi = 600,
       device = ragg::agg_png, bg = "white")
message("\nWrote ngn2_exon31_psi_samples.csv, ngn2_exon31_psi_by_line.csv, ",
        "ngn2_exon31_psi_by_day.csv and ngn2_exon31_psi.png in ", getwd())
