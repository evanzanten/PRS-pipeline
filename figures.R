#!/usr/bin/env Rscript
#
# Figures for the QC steps (sections 1.1 to 1.8 of the README).
#
# Run from the README at the point in the pipeline where each figure is produced:
#
#   Rscript figures.R missingness   # 1.1  variant and sample missingness
#   Rscript figures.R sex           # 1.2  X chromosome F statistic
#   Rscript figures.R pca_het       # 1.3  ancestry clusters and heterozygosity
#   Rscript figures.R kinship       # 1.4  KING kinship coefficients
#   Rscript figures.R harmonization # 1.8  allele frequencies against 1000 Genomes
#   Rscript figures.R pca_reference # 1.9  PCA of our cohort together with 1000 Genomes
#   Rscript figures.R imputation    # 2    imputation quality, on its own and against
#                                   #      allele frequency
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

  for (f in list(c("missingness.vmiss", "missingness_variant.png"),
                 c("missingness.smiss", "missingness_sample.png"))) {
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

  d <- read_plink("sexcheck.sexcheck", sep = "")
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

  # The clustering itself is pca_clusters.R's job, and it has already run by this
  # point in the pipeline. We read the cluster it assigned to each person rather
  # than repeating the k-means here: two copies of the same clustering would
  # drift apart, and only one of them decides who gets removed.
  pca <- read_plink("rough_pca.eigenvec")
  het <- read_plink("full_cohort_het.het")
  cl  <- read_plink("clusters.txt")
  d   <- merge(merge(pca, het, by = c("#FID", "IID")), cl, by = c("#FID", "IID"))
  d$cluster <- factor(d$cluster)

  d$het_rate <- (d[["OBS_CT"]] - d[["O(HOM)"]]) / d[["OBS_CT"]]
  mu <- tapply(d$het_rate, d$cluster, mean)
  sd <- tapply(d$het_rate, d$cluster, stats::sd)
  # The same one-sided rule pca_clusters.R applies, recomputed only so the
  # dashed line can be drawn in the right place. Nothing is written here.
  d$outlier <- (d$het_rate - mu[d$cluster]) > 3 * sd[d$cluster]

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
  cat("heterozygosity outliers shown:", sum(d$outlier), "\n")
}


###KINSHIP###
#Counts per relationship class, from Table 2 of the paper.

if (step == "kinship") {

  d <- read_plink("kinship.kin0")

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
#HRC-1000G-check-bim.pl writes FreqPlot-MAF-1000G.txt with, for every variant it
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

  d <- read.table(file.path(QC, "FreqPlot-MAF-1000G.txt"),
                  col.names = c("id", "ref_af", "our_af", "diff", "code"))
  excluded <- readLines(file.path(QC, "Exclude-MAF-1000G.txt"))

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


###PCA WITH REFERENCE###
#PC1 against PC2 for our cohort merged with 1000 Genomes. The reference individuals are
#round and coloured by superpopulation, using the panel file from 1000 Genomes (KG_PANEL).
#Everyone not in the panel file is one of ours; they are drawn as translucent dark grey
#triangles on top, so the reference individuals underneath stay visible.
#
#This is the one plot with axis labels and text: the cluster names are what makes it
#readable, so PC1 and PC2 are named as well.

if (step == "pca_reference") {

  pca <- read_plink("pca_with_reference.eigenvec")

  # The panel file has the columns sample, pop, super_pop and gender
  panel <- read.table(Sys.getenv("KG_PANEL"), header = TRUE, sep = "\t", fill = TRUE)
  pca$group <- panel$super_pop[match(pca$IID, panel$sample)]
  pca$group[is.na(pca$group)] <- "cohort"

  # Fixed colour per superpopulation, checked for colour blindness
  SUPERPOP <- c(AFR = "#2a78d6", AMR = "#e34948", EAS = "#1baf7a",
                EUR = "#eda100", SAS = "#4a3aa7")

  reference <- pca[pca$group != "cohort", ]
  cohort    <- pca[pca$group == "cohort", ]

  # One label per superpopulation, at the middle of its cluster. The clusters sit close
  # together on the right, so each label is shifted by hand into free space; adjust these
  # two columns if your own clusters end up somewhere else.
  labels <- aggregate(cbind(PC1, PC2) ~ group, reference, mean)
  nudge <- data.frame(group = c("AFR",  "AMR",   "EAS",   "EUR",  "SAS"),
                      dx    = c( 0.000, -0.011, -0.007,  -0.004, -0.008),
                      dy    = c( 0.005, -0.005,  0.000,   0.005, -0.004))
  labels <- merge(labels, nudge, by = "group")

  p <- ggplot() +
    geom_point(data = reference, aes(x = PC1, y = PC2, colour = group),
               size = 0.7, alpha = 0.7, stroke = 0) +
    geom_point(data = cohort, aes(x = PC1, y = PC2),
               shape = 17, colour = "#3a3a3a", size = 0.9, alpha = 0.35, stroke = 0) +
    geom_text(data = labels, aes(x = PC1 + dx, y = PC2 + dy, label = group, colour = group),
              size = 3, show.legend = FALSE) +
    scale_colour_manual(values = SUPERPOP) +
    labs(x = "PC1", y = "PC2") +
    qc_theme() +
    theme(axis.title = element_text(size = 9))
  save_fig(p, "pca_reference.png", width = 4, height = 3.6)

  # Besides the plot, put a number on where our individuals fall: give each of them the
  # name of the cluster whose centre is closest. The centre of a cluster is the average
  # PC1 and PC2 of its individuals. This is a rough summary and not an ancestry
  # assignment: an admixed individual lies between clusters and still ends up in one.
  centres <- aggregate(cbind(PC1, PC2) ~ group, data = reference, FUN = mean)
  closest <- apply(cohort[, c("PC1", "PC2")], 1, function(person) {
    distance <- (centres$PC1 - person["PC1"])^2 + (centres$PC2 - person["PC2"])^2
    centres$group[which.min(distance)]
  })

  cat("individuals per group:\n"); print(table(pca$group))
  cat("\nour individuals, by closest cluster:\n"); print(table(closest))
  cat("\nas a percentage:\n"); print(round(100 * prop.table(table(closest))))
}


#Imputation quality. Beagle scores every variant it imputed with DR2, between 0 and 1,
#which is its own estimate of how well it managed. Two things are worth seeing: how the
#scores are distributed, and how they depend on allele frequency, since a rare variant
#sits on few haplotypes in the reference panel and is the hardest to impute.
#
#Read the scores out of the imputed VCFs first, before any filtering, which is one line
#per variant and no genotypes:
#
#  for chr in {1..22}; do
#    bcftools query -f '%INFO/DR2\t%INFO/AF\n' "$OUT/chr${chr}_imputed.vcf.gz"
#  done | awk -F'\t' '{n = split($1, dr2, ","); split($2, af, ",")
#                       for (i = 1; i <= n; i++) print dr2[i]"\t"af[i]}' \
#     > "$OUT/imputation_quality.txt"
#
#The awk is there because a variant with several alternate alleles carries one DR2 and
#one frequency per allele, comma separated, which would otherwise be read as text.

if (step == "imputation") {

  CUTOFF_DR2 <- 0.8
  CUTOFF_MAF <- 0.01

  d <- data.table::fread(file.path(QC, "imputation_quality.txt"),
                         col.names = c("dr2", "af"))
  d <- d[is.finite(d$dr2) & is.finite(d$af), ]

  # Frequency of the rarer allele, which is what the filter in section 2 works on.
  d$maf <- pmin(d$af, 1 - d$af)

  cat("variants:", nrow(d), "\n")
  cat("passing DR2 >=", CUTOFF_DR2, ":", sum(d$dr2 >= CUTOFF_DR2), "\n")
  cat("passing both filters:",
      sum(d$dr2 >= CUTOFF_DR2 & d$maf >= CUTOFF_MAF), "\n")

  h <- binned(d$dr2, bins = 100)
  p <- ggplot(h, aes(x = x, y = n)) +
    geom_col(width = h$width, fill = COLORS[1]) +
    geom_vline(xintercept = CUTOFF_DR2, linewidth = 0.3, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"),
                       expand = expansion(mult = c(0, 0.05))) +
    qc_theme()
  save_fig(p, "imputation_dr2.png")

  # Thirty million points is more than a figure can show or a reader can take in, so
  # count them into a grid and shade each cell by how many landed in it. Frequency goes
  # on a log axis because most variants are rare, and the axis is plotted as log10 with
  # relabelled breaks rather than with scale_x_log10, so that the cells stay evenly
  # spaced. DR2 is reported to two decimals, so its own resolution sets the row height;
  # cutting it finer would leave empty rows striped through the picture.
  # With 880 people the rarest frequency the data can really express is about 1 in
  # 1,760 haplotypes, and Beagle writes frequencies to four decimals, so below roughly
  # 0.001 the only values that occur are rounding artefacts sitting in isolated
  # columns. The axis therefore starts there.
  d <- d[d$maf >= 1e-3, ]
  x_bins <- 80
  x_edge <- seq(log10(1e-3), log10(0.5), length.out = x_bins + 1)
  x_mid  <- (x_edge[-1] + x_edge[-length(x_edge)]) / 2
  cell <- table(findInterval(log10(d$maf), x_edge, rightmost.closed = TRUE),
                round(d$dr2, 2))
  grid <- as.data.frame(cell, stringsAsFactors = FALSE)
  names(grid) <- c("xi", "dr2", "n")
  grid <- grid[grid$n > 0, ]
  grid$x <- x_mid[as.integer(grid$xi)]
  grid$y <- as.numeric(grid$dr2)

  p <- ggplot(grid, aes(x = x, y = y, fill = n)) +
    geom_tile(width = diff(x_edge)[1], height = 0.01) +
    geom_hline(yintercept = CUTOFF_DR2, linewidth = 0.3, linetype = "dashed") +
    geom_vline(xintercept = log10(CUTOFF_MAF), linewidth = 0.3, linetype = "dashed") +
    # Most cells hold a handful of variants and a few hold hundreds of thousands, so
    # shade by the logarithm or everything but the densest band looks empty.
    scale_fill_gradient(low = "#dce9f8", high = COLORS[1], trans = "log10") +
    scale_x_continuous(breaks = log10(c(1e-3, 1e-2, 1e-1, 0.5)),
                       labels = c("0.001", "0.01", "0.1", "0.5"), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    qc_theme()
  save_fig(p, "imputation_dr2_maf.png", width = 4.6)
}
