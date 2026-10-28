#!/usr/bin/env Rscript
#
# Figures for the QC steps (protocol paper sections 3.1 to 3.8).
#
# Called by qc.sh at the point in the pipeline where each figure is produced:
#
#   Rscript figures.R missingness   # 3.1  variant and sample missingness
#   Rscript figures.R sex           # 3.2  X chromosome F statistic
#   Rscript figures.R pca_het       # 3.3  ancestry clusters and heterozygosity
#   Rscript figures.R kinship       # 3.4  KING kinship coefficients
#   Rscript figures.R harmonization # 3.8  allele frequencies against 1000 Genomes
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
  # Only the high tail. Excess heterozygosity means the sample is a mixture of
  # two people's DNA, which is a technical failure. A low value means the
  # opposite, long runs of homozygosity, which in a cohort like ours is
  # consanguinity rather than a bad sample, so we keep those individuals.
  d$outlier <- (d$het_rate - mu[d$cluster]) > 3 * sd[d$cluster]

  write.table(d[d$outlier, c("#FID", "IID")],
              file.path(QC, "06_het_outliers.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  p <- ggplot(d, aes(x = PC1, y = PC2, fill = cluster)) +
    geom_point(shape = 21, size = 1.1, colour = "white", stroke = 0.15) +
    scale_fill_manual(values = COLORS) +
    qc_theme()
  save_fig(p, "pca_clusters.png", height = 3.4)

  # One histogram per cluster, stacked so they share the x axis and the cluster
  # means can be compared. Each panel has its own y axis, because the clusters
  # differ a lot in size. The solid line is the cluster mean and the dashed line
  # is 3 SD above it. Bars to the right of the dashed line are the outliers.
  lines <- data.frame(cluster = factor(names(mu), levels = levels(d$cluster)),
                      mean = mu, upper = mu + 3 * sd)

  p <- ggplot(d, aes(x = het_rate, fill = cluster)) +
    geom_histogram(bins = 50, colour = "white", linewidth = 0.1) +
    geom_vline(data = lines, aes(xintercept = mean),
               colour = "#3a3a3a", linewidth = 0.3) +
    geom_vline(data = lines, aes(xintercept = upper),
               colour = "#8a8a8a", linewidth = 0.3, linetype = "dashed") +
    facet_wrap(~ cluster, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = COLORS) +
    # Whole-number counts only, with the bars sitting on the x axis
    scale_y_continuous(breaks = breaks_pretty(n = 4),
                       expand = expansion(mult = c(0, 0.05))) +
    qc_theme() +
    theme(strip.text = element_blank(), strip.background = element_blank(),
          panel.spacing = unit(1.5, "lines"))
  save_fig(p, "heterozygosity.png", height = 4.5)

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


###HARMONIZATION###
#Allele frequencies in our cohort against 1000 Genomes, before and after harmonization.
#
#HRC-1000G-check-bim.pl writes FreqPlot-12_maf-1000G.txt with, for every variant it
#could match: the variant ID, the 1000 Genomes frequency of the allele in our A1
#column, the frequency of that allele in our data, the difference, and a code for how
#our alleles line up with 1000 Genomes (1 = our A1 is the 1000 Genomes reference
#allele, 2 = our A1 is the 1000 Genomes alternative allele, 3 and 4 = the same after a
#strand flip).
#
#The x axis is the same in both panels: the frequency of the ALT allele as 1000 Genomes
#defines it.
#  Before: y is the allele PLINK called ALT in our data (A1). Where that is the same
#  allele as in 1000 Genomes, points lie on the diagonal. Where PLINK had REF and ALT
#  swapped, our ALT is the other allele, so its frequency is 1 minus the reference one:
#  the anti-diagonal. Our ALT is almost always the rarer allele, so y stays below ~0.5.
#  After: y is the frequency in our data of the allele 1000 Genomes calls ALT, which is
#  what our ALT allele is once step 4 has set the reference alleles. Now every point is
#  on the diagonal.

if (step == "harmonization") {

  library(grid)
  library(gtable)

  d <- read.table(file.path(QC, "FreqPlot-12_maf-1000G.txt"),
                  col.names = c("id", "ref_af", "our_af", "diff", "code"))
  excluded <- readLines(file.path(QC, "Exclude-12_maf-1000G.txt"))

  # Variants the tool removes for a frequency difference > 0.2 are left out: they do not
  # make it into the harmonized data anyway
  d <- d[!(d$id %in% excluded), ]

  swapped <- d$code %in% c(1, 3)

  # 1000 Genomes frequency of its ALT allele, whatever our A1 is
  d$kg_alt_af <- ifelse(swapped, 1 - d$ref_af, d$ref_af)

  # Our frequency of that same allele
  d$our_af_after <- ifelse(swapped, 1 - d$our_af, d$our_af)

  d$group <- factor(ifelse(swapped, "swapped", "concordant"),
                    levels = c("concordant", "swapped"))

  panel <- function(y, ymax) {
    ggplot(d, aes(x = kg_alt_af, y = .data[[y]], colour = group)) +
      geom_point(size = 0.15, alpha = 0.4, stroke = 0) +
      scale_colour_manual(values = c(concordant = COLORS[1], swapped = COLORS[2])) +
      scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, ymax + 0.01), clip = "off") +
      qc_theme() +
      # room on the right so the last x axis label ("1.00") is not cut off
      theme(plot.margin = margin(5, 14, 5, 5))
  }

  # Before stops at 0.6: our ALT frequencies are almost all below 0.5 and the highest
  # is 0.53, so every variant is still shown. After runs to 1.
  before <- ggplotGrob(panel("our_af", 0.6))
  after  <- ggplotGrob(panel("our_af_after", 1))

  # Put the two panels side by side with their x axes on one line. Both get an empty
  # row above the panel: the before panel leaves it empty, the after panel stretches
  # over it. Heights are in proportion to the axis ranges (0.62 and 1.02 including the
  # small margin at each end) so both panels use the same scale as the x axis (1.02).
  add_top_row <- function(g, stretch) {
    row <- g$layout$t[g$layout$name == "panel"]
    g <- gtable_add_rows(g, unit(0.40 / 1.02, "null"), pos = row - 1)
    g$heights[row + 1] <- unit(0.62 / 1.02, "null")
    if (stretch) {
      # the panel and its y axis start in the new row instead of below it
      on_panel_row <- g$layout$t == row + 1 & g$layout$b == row + 1
      g$layout$t[on_panel_row] <- row
    }
    g
  }
  g <- cbind(add_top_row(before, FALSE), add_top_row(after, TRUE), size = "first")

  # Choose the figure size so one unit on the x axis is exactly as long as one unit
  # on the y axis, which keeps the diagonals at 45 degrees
  side <- 3   # inches per unit of allele frequency
  pdf(NULL)
  fixed_w <- sum(convertWidth(g$widths[unitType(g$widths) != "null"], "in", valueOnly = TRUE))
  fixed_h <- sum(convertHeight(g$heights[unitType(g$heights) != "null"], "in", valueOnly = TRUE))
  invisible(dev.off())

  # bg = "white": otherwise the empty space above the before panel is transparent,
  # which many viewers show as black
  ggsave(file.path(FIG, "harmonization_frequencies.png"), g,
         width = fixed_w + 2 * side, height = fixed_h + side, dpi = 300, bg = "white")

  print(table(d$group))
}
