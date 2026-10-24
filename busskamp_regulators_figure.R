# Two-panel figure of NF1 and candidate exon 31 splicing regulators, styled to match
# the GTEx/NGN2 regulator figure (Figures 3.16 and 3.17).
#
# (A) Heatmap of log2(TPM + 1) with TPM printed in each cell:
#     GTEx (median TPM) | Busskamp et al. (2014) by day | NGN2 neurons by day
#     The GTEx and NGN2 blocks are added only if their summary CSVs are in the folder.
# (B) Busskamp time course per gene (mean +/- SD across replicates).
#
# Inputs (same folder = R working directory):
#   busskamp_nf1_samples.csv              from nf1_busskamp_recount3.R    (required)
#   nf1_regulators_gtex_summary.csv       from nf1_regulators_gtex_recount3.R  (optional)
#   nf1_regulators_neurons_summary.csv    from nf1_regulators_neurons_nfcore.R (optional)
#
# Run in the R console:   source(file.choose(), echo = TRUE)
# Writes two separate 600 dpi figures: busskamp_regulators_heatmap.png (A) and
# busskamp_regulators_timecourse.png (B).

for (p in c("ggplot2", "ragg", "systemfonts")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)

GENES      <- c("NF1", "REST", "CELF3", "CELF4", "CELF5", "CELF6", "ELAVL2", "ELAVL3", "ELAVL4")
FONT       <- "Lato"   # used if installed, otherwise the default sans font
HEADER_BK  <- "Busskamp et al. (2014)"
HEADER_NGN <- "NGN2 neurons"
USE_GTEX   <- TRUE     # add GTEx columns if nf1_regulators_gtex_summary.csv is present
USE_NGN2   <- TRUE     # add NGN2 data to A and B if nf1_regulators_neurons_summary.csv is present
X_SCALE    <- "pseudolog" # "pseudolog" spreads out Days 0-7, "linear" matches Figure 3.16

# Fall back to the default sans font unless this font has an italic face (gene names are italic)
sf <- systemfonts::system_fonts()
if (!any(sf$family == FONT & sf$italic)) FONT <- ""

# 0. Locate the optional GTEx and NGN2 summaries ---------------------------------------
# Looks in the working directory first, then searches OneDrive and your home folder.
# Or set the full paths by hand, e.g. GTEX_FILE <- "C:/Users/you/Documents/nf1_regulators_gtex_summary.csv"
GTEX_FILE <- NA
NGN2_FILE <- NA

find_input <- function(name, path) {
  if (!is.na(path) && file.exists(path)) return(path)
  if (file.exists(name)) return(normalizePath(name))
  roots <- unique(Filter(nzchar, c(Sys.getenv("OneDrive"), Sys.getenv("OneDriveConsumer"),
                                   path.expand("~"), dirname(getwd()))))
  for (r in roots[dir.exists(roots)]) {
    hit <- list.files(r, pattern = paste0("^", gsub(".", "\\.", name, fixed = TRUE), "$"),
                      recursive = TRUE, full.names = TRUE)
    if (length(hit)) return(hit[order(file.mtime(hit), decreasing = TRUE)][1])
  }
  NA
}
GTEX_FILE <- find_input("nf1_regulators_gtex_summary.csv", GTEX_FILE)
NGN2_FILE <- find_input("nf1_regulators_neurons_summary.csv", NGN2_FILE)
message("GTEx summary: ", ifelse(is.na(GTEX_FILE), "NOT FOUND (GTEx block omitted)", GTEX_FILE))
message("NGN2 summary: ", ifelse(is.na(NGN2_FILE), "NOT FOUND (NGN2 data omitted)", NGN2_FILE))

# 1. Busskamp per-sample TPM -> mean and SD per day -------------------------------------
d <- read.csv("busskamp_nf1_samples.csv", check.names = FALSE)
d <- d[!is.na(d$day), ]
genes <- intersect(GENES, names(d))
long <- do.call(rbind, lapply(genes, function(g) data.frame(gene = g, day = d$day, tpm = d[[g]])))
bk <- aggregate(tpm ~ gene + day, data = long,
                FUN = function(x) c(mean = mean(x), sd = if (length(x) > 1) sd(x) else 0))
bk <- do.call(data.frame, bk)
names(bk) <- c("gene", "day", "mean", "sd")

# 2. Heatmap table ---------------------------------------------------------------------
blocks <- list()
if (USE_GTEX && !is.na(GTEX_FILE)) {
  gt <- read.csv(GTEX_FILE)
  tis <- c("Cells - EBV-transformed lymphocytes" = "EBV\nlymphocytes", "Whole Blood" = "Whole\nblood")
  gt <- gt[gt$tissue %in% names(tis), ]
  blocks$GTEx <- data.frame(gene = gt$gene, col = tis[gt$tissue], value = gt$median_tpm,
                            order = match(gt$tissue, names(tis)))
}
blocks[[HEADER_BK]] <- data.frame(gene = bk$gene, col = paste("Day", bk$day),
                                  value = bk$mean, order = bk$day)
if (USE_NGN2 && !is.na(NGN2_FILE)) {
  ng <- read.csv(NGN2_FILE)
  ng <- aggregate(mean_tpm ~ gene + day, data = ng, FUN = mean)   # mean across lines
  blocks[[HEADER_NGN]] <- data.frame(gene = ng$gene, col = paste0("Day ", ng$day, " "),
                                     value = ng$mean_tpm, order = ng$day)
}
# Note: NGN2 columns carry a trailing space so "Day 7" etc. never clash with Busskamp days.

hm <- do.call(rbind, lapply(seq_along(blocks), function(i) cbind(blocks[[i]], block = i)))
hm <- hm[hm$gene %in% genes, ]
col_order <- unique(hm[order(hm$block, hm$order), "col"])
hm$x <- match(hm$col, col_order)
hm$y <- match(hm$gene, rev(genes))
hm$z <- log2(hm$value + 1)
hm$label <- ifelse(hm$value >= 10, sprintf("%.0f", hm$value),
                   ifelse(hm$value >= 0.1, sprintf("%.1f", hm$value), "<0.1"))

ends   <- tapply(hm$x, hm$block, max)
starts <- tapply(hm$x, hm$block, min)
seps   <- head(ends, -1) + 0.5
heads  <- data.frame(x = (starts + ends) / 2, label = names(blocks))
zmax   <- max(9.5, max(hm$z, na.rm = TRUE))
pt     <- ggplot2::.pt

pA <- ggplot(hm, aes(x, y)) +
  geom_tile(aes(fill = z)) +
  geom_text(aes(label = label, colour = z > 5.5), size = 9 / pt, family = FONT) +
  geom_vline(xintercept = seps, linewidth = 0.9, colour = "black") +
  annotate("text", x = heads$x, y = length(genes) + 0.95, label = heads$label,
           fontface = "bold", size = 11 / pt, family = FONT) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, zmax),
                       breaks = seq(0, 8, 2), name = expression(log[2]*"(TPM + 1)")) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_x_continuous(breaks = seq_along(col_order), labels = trimws(col_order),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(genes), labels = rev(genes), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_colourbar(title.position = "right", barheight = unit(0.62, "npc"),
                                barwidth = unit(0.35, "cm"), ticks = FALSE)) +
  theme_void(base_family = FONT) +
  theme(axis.text.x = element_text(size = 10, margin = margin(t = 4), lineheight = 0.9),
        axis.text.y = element_text(size = 11, face = "italic", hjust = 1, margin = margin(r = 4)),
        legend.title = element_text(size = 10, angle = 90, hjust = 0.5),
        legend.text = element_text(size = 9),
        plot.margin = margin(28, 6, 6, 6))

# 3. Time course: Busskamp (Days 0-4) and NGN2 lines (Days 7-35) on shared axes -------
tc <- data.frame(series = HEADER_BK, gene = bk$gene, day = bk$day, mean = bk$mean, sd = bk$sd)
if (USE_NGN2 && !is.na(NGN2_FILE)) {
  nl <- read.csv(NGN2_FILE)
  nl$sd_tpm[is.na(nl$sd_tpm)] <- 0
  tc <- rbind(tc, data.frame(series = sub("_NGN2$", "", nl$line), gene = nl$gene, day = nl$day,
                             mean = nl$mean_tpm, sd = nl$sd_tpm))
}
tc <- tc[tc$gene %in% genes, ]
tc$gene   <- factor(tc$gene, levels = genes)
series    <- c(HEADER_BK, sort(setdiff(unique(tc$series), HEADER_BK)))
tc$series <- factor(tc$series, levels = series)
# Styles: Busskamp in blue triangles, NGN2 lines as in Figure 3.16
cols   <- setNames(c("#2171B5", "black", "grey50", "grey30")[seq_along(series)], series)
shapes <- setNames(c(17, 16, 22, 15)[seq_along(series)], series)
fills  <- setNames(c("#2171B5", "black", "white", "grey30")[seq_along(series)], series)
days   <- sort(unique(tc$day))
# Tick marks chosen so labels don't overlap on either axis type
xbreaks <- if (X_SCALE == "pseudolog") c(0, 1, 4, 7, 14, 35) else
  if (max(days) > 7) c(0, 7, 14, 21, 28, 35) else days
xbreaks <- xbreaks[xbreaks <= max(days)]
xscale <- if (X_SCALE == "pseudolog") {
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = xbreaks)
} else {
  scale_x_continuous(breaks = xbreaks)
}
dodge <- if (X_SCALE == "pseudolog") 0.06 else 0.6   # errorbar cap width

pB <- ggplot(tc, aes(day, mean, colour = series, shape = series, fill = series, group = series)) +
  geom_errorbar(data = tc[tc$sd > 0, ], aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
                width = dodge, linewidth = 0.45) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8, stroke = 0.7) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_fill_manual(values = fills, name = NULL) +
  xscale +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  # scales = "free" keeps x-axis ticks on every panel, as in Figure 3.16
  facet_wrap(~ gene, ncol = 3, scales = "free") +
  labs(x = "Day", y = "TPM") +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 11),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 9),
        panel.spacing.y = unit(1.1, "lines"),
        legend.position = if (length(series) > 1) "top" else "none",
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, -4, 0))

# 4. Save the two panels as separate figures -------------------------------------------
n_cols <- length(col_order)
ggsave("busskamp_regulators_heatmap.png", pA, width = 2.6 + 0.62 * n_cols, height = 6.2,
       dpi = 600, device = ragg::agg_png, bg = "white")
ggsave("busskamp_regulators_timecourse.png", pB, width = 6.6, height = 6.4,
       dpi = 600, device = ragg::agg_png, bg = "white")
message("Wrote busskamp_regulators_heatmap.png and busskamp_regulators_timecourse.png in ", getwd())
