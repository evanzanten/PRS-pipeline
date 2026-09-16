#!/usr/bin/env Rscript
#
# Figures for the QC steps (protocol paper sections 3.1 to 3.4).
#
# Called by qc.sh at the point in the pipeline where each figure is produced:
#
#   Rscript figures.R missingness   # 3.1  variant and sample missingness
#   Rscript figures.R sex           # 3.2  X chromosome F statistic
#   Rscript figures.R pca_het       # 3.3  ancestry clusters and heterozygosity
#   Rscript figures.R kinship       # 3.4  KING kinship coefficients
#
# The pca_het step also does the clustering and writes the list of
# heterozygosity outliers that the pipeline removes.
#
# Plots carry no titles, captions or axis labels, so they can go straight into
# the paper with the caption written alongside them. To label the axes, drop the
# axis.title line out of qc_theme() and add labs(x = ..., y = ...) per plot.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

step <- commandArgs(trailingOnly = TRUE)[1]

QC  <- Sys.getenv("OUT", unset = ".")
FIG <- Sys.getenv("FIG", unset = "../figures")

# Categorical colours, in this fixed order, for the ancestry clusters
COLORS <- c("#2a78d6", "#eb6834", "#1baf7a")

qc_theme <- function() {
  theme_classic(base_size = 9) +
    theme(axis.title  = element_blank(),
          plot.title  = element_blank(),
          legend.position = "none",
          axis.line   = element_line(linewidth = 0.3),
          axis.ticks  = element_line(linewidth = 0.3))
}

save_fig <- function(plot, file, width = 4, height = 3) {
  ggsave(file.path(FIG, file), plot, width = width, height = height, dpi = 300)
}

# plink writes tab-separated tables whose header starts with '#', and some
# column names contain brackets, so keep the names exactly as they are.
read_plink <- function(file, sep = "\t") {
  read.table(file.path(QC, file), header = TRUE, sep = sep,
             comment.char = "", check.names = FALSE)
}

# Bin a vector and drop the empty bins. geom_histogram would keep them, and on a
# log axis a count of zero becomes -Inf, which ggplot warns about on every plot.
binned <- function(x, bins = 80) {
  h <- hist(x, breaks = seq(min(x), max(x), length.out = bins + 1), plot = FALSE)
  data.frame(x = h$mids, n = h$counts, width = diff(h$breaks)[1])[h$counts > 0, ]
}


###MISSINGNESS###
#Log counts, because nearly every variant and sample sits in the first bin and a
#linear axis would hide the tail that the threshold has to be chosen against.

if (step == "missingness") {

  for (f in list(c("01_missing_raw.vmiss", "missingness_variant.png"),
                 c("01_missing_raw.smiss", "missingness_sample.png"))) {
    d <- read_plink(f[1])
    p <- ggplot(binned(d$F_MISS), aes(x = x, y = n, width = width)) +
      geom_col(fill = COLORS[1]) +
      scale_y_log10(labels = label_comma()) +
      qc_theme()
    save_fig(p, f[2])
  }
}


###SEX###
#Histogram of the X chromosome F statistic, with the 0.2 and 0.8 cutoffs marked.

if (step == "sex") {

  d <- read_plink("04_sexcheck.sexcheck", sep = "")
  p <- ggplot(d, aes(x = F)) +
    geom_histogram(bins = 80, fill = COLORS[1]) +
    geom_vline(xintercept = c(0.2, 0.8), colour = "#8a8a8a",
               linewidth = 0.3, linetype = "dashed") +
    qc_theme()
  save_fig(p, "sex_fstat.png")
}


###PCA AND HETEROZYGOSITY###
#Cluster on PC1/PC2, then flag individuals more than 3 SD from the mean
#heterozygosity of their own cluster rather than of the whole cohort.

if (step == "pca_het") {

  K <- 3   # broad ancestry clusters

  pca <- read_plink("06_pca.eigenvec")
  het <- read_plink("06_het.het")
  d   <- merge(pca, het, by = c("#FID", "IID"))

  # Standardise so that PC1 does not dominate the distance purely by scale
  z <- scale(as.matrix(d[, c("PC1", "PC2")]))
  set.seed(1)
  d$cluster <- factor(kmeans(z, centers = K, nstart = 25)$cluster)

  d$het_rate <- (d[["OBS_CT"]] - d[["O(HOM)"]]) / d[["OBS_CT"]]
  mu <- tapply(d$het_rate, d$cluster, mean)
  sd <- tapply(d$het_rate, d$cluster, stats::sd)
  d$outlier <- abs(d$het_rate - mu[d$cluster]) > 3 * sd[d$cluster]

  write.table(d[d$outlier, c("#FID", "IID")],
              file.path(QC, "06_het_outliers.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  p <- ggplot(d, aes(x = PC1, y = PC2, fill = cluster)) +
    geom_point(shape = 21, size = 1.1, colour = "white", stroke = 0.15) +
    scale_fill_manual(values = COLORS) +
    qc_theme()
  save_fig(p, "pca_clusters.png", height = 3.4)

  # Outlines rather than filled bars: three translucent fills would overlap
  # into colours that look like clusters of their own.
  p <- ggplot(d, aes(x = het_rate, colour = cluster)) +
    geom_freqpoly(bins = 45, linewidth = 0.5) +
    scale_colour_manual(values = COLORS) +
    qc_theme()
  save_fig(p, "heterozygosity.png")

  cat("cluster sizes:", table(d$cluster), "\n")
  cat("heterozygosity outliers removed:", sum(d$outlier), "\n")
}


###KINSHIP###
#Counts per relationship class, from Table 2 of the paper.

if (step == "kinship") {

  d <- read_plink("08_king.kin0")

  labels <- c("<0.04", "0.04-0.09", "0.09-0.18", "0.18-0.35", ">0.35")
  d$class <- cut(d$KINSHIP, breaks = c(-Inf, 0.04, 0.09, 0.18, 0.35, Inf),
                 labels = labels)
  counts <- as.data.frame(table(d$class))
  names(counts) <- c("class", "n")

  # One ordered variable, so a single hue light to dark rather than five colours
  shades <- c("#bcd4f0", "#8fb8e6", "#5f9adc", "#2a78d6", "#1b4f8c")

  p <- ggplot(counts, aes(x = class, y = n, fill = class)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = shades) +
    scale_y_log10(labels = label_comma()) +   # the unrelated bin dwarfs the rest
    qc_theme() +
    theme(axis.text.x = element_text(size = 7))
  save_fig(p, "kinship.png")

  print(counts, row.names = FALSE)
}
