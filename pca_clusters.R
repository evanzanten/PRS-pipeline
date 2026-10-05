#!/usr/bin/env Rscript
#
# This is the only R this pipeline needs for the QC steps. It has two jobs, and you pick
# one on the command line:
#
#   Rscript pca_clusters.R clusters    # step 1.3
#   Rscript pca_clusters.R ancestry    # step 1.9
#
# The "clusters" job splits the cohort into broad ancestry clusters using the
# preliminary PCA, then flags heterozygosity outliers inside each cluster. It
# writes the two files later steps read:
#
#   clusters.txt       FID, IID and cluster number for everybody. Step 1.6
#                      reads this to test Hardy-Weinberg per cluster.
#   het_outliers.txt   The individuals to drop, as FID and IID with no header,
#                      which is the layout PLINK's --remove expects.
#
# The "ancestry" job takes the PCA that was run together with 1000 Genomes and
# tells you which superpopulation each of your individuals sits closest to. It
# only prints a table; nothing later in the pipeline reads it. You need it to
# choose an LD reference panel and a GWAS arm in sections 3.4 to 3.6.
#
# Only base R is used, so there are no packages to install.

# Where the PLINK output lives. Set OUT in your shell, or replace this with the
# literal path to your working directory.
OUT  <- Sys.getenv("OUT", unset = ".")

# For the ancestry job only: the 1000 Genomes panel file, which lists each
# reference sample and the superpopulation it belongs to.
KG_PANEL <- Sys.getenv("KG_PANEL", unset = "integrated_call_samples_v3.20130502.ALL.panel")

job <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(job)) stop("say which job to run: clusters or ancestry")

# PLINK writes tab separated tables whose header line starts with a '#', and
# some column names contain brackets, so read the names exactly as they are.
read_plink <- function(file) {
  read.table(file.path(OUT, file), header = TRUE, sep = "\t",
             comment.char = "", check.names = FALSE)
}


### STEP 1.3: ANCESTRY CLUSTERS AND HETEROZYGOSITY OUTLIERS ###

if (job == "clusters") {

  K <- 3   # how many clusters to split into, see the note at the bottom

  pca <- read_plink("rough_pca.eigenvec")
  het <- read_plink("full_cohort_het.het")

  # Both files have one row per individual. Merging on the sample ID rather than
  # trusting both files to be in the same order is safer.
  d <- merge(pca, het, by = c("#FID", "IID"))

  # Put PC1 and PC2 on the same scale first, so that neither one dominates the
  # distance purely because it covers a wider numeric range. PC1 always spreads
  # wider than PC2, so without this the clustering would almost ignore PC2.
  z <- scale(as.matrix(d[, c("PC1", "PC2")]))

  # k-means starts from randomly chosen centres, so fix the seed to get the same
  # clusters every time. nstart = 25 repeats it from 25 different starts and
  # keeps the tightest result, which stops one unlucky start deciding the split.
  set.seed(1)
  d$cluster <- kmeans(z, centers = K, nstart = 25)$cluster

  # Heterozygosity rate: the share of a person's called genotypes that are not
  # homozygous. PLINK gives the number of homozygous calls, O(HOM), and the
  # number of calls it looked at, OBS_CT.
  d$het_rate <- (d[["OBS_CT"]] - d[["O(HOM)"]]) / d[["OBS_CT"]]

  # Compare each person with the mean of their own cluster, not of the whole
  # cohort. Ancestry groups genuinely differ in heterozygosity, so a cohort-wide
  # mean would flag whole groups instead of the bad samples we are after.
  mu <- tapply(d$het_rate, d$cluster, mean)
  sd <- tapply(d$het_rate, d$cluster, stats::sd)
  d$outlier <- abs(d$het_rate - mu[d$cluster]) > 3 * sd[d$cluster]

  # clusters.txt keeps its header, because the awk command in step 1.6 skips the
  # first line and reads the cluster number out of the third column.
  write.table(d[, c("#FID", "IID", "cluster")],
              file.path(OUT, "clusters.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # het_outliers.txt gets no header, because PLINK's --remove reads every line
  # in the file as a sample.
  write.table(d[d$outlier, c("#FID", "IID")],
              file.path(OUT, "het_outliers.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat("individuals per cluster:", table(d$cluster), "\n")
  cat("heterozygosity outliers flagged:", sum(d$outlier), "\n")
  cat("mean heterozygosity per cluster:", round(mu, 4), "\n")
}


### STEP 1.9: WHICH SUPERPOPULATION EACH INDIVIDUAL SITS CLOSEST TO ###

if (job == "ancestry") {

  pca <- read_plink("pca_with_reference.eigenvec")

  # The panel file has the columns sample, pop, super_pop and gender. Anyone in
  # the PCA who is not in it is one of our own individuals.
  panel <- read.table(KG_PANEL, header = TRUE, sep = "\t", fill = TRUE)
  pca$group <- panel$super_pop[match(pca$IID, panel$sample)]
  pca$group[is.na(pca$group)] <- "cohort"

  reference <- pca[pca$group != "cohort", ]
  cohort    <- pca[pca$group == "cohort", ]

  # The centre of a superpopulation is just the average PC1 and PC2 of its
  # reference individuals.
  centres <- aggregate(cbind(PC1, PC2) ~ group, data = reference, FUN = mean)

  # Give each of our individuals the name of the nearest centre. This is a rough
  # summary, not an ancestry assignment: an admixed person sits between two
  # clusters and still gets put in one of them.
  closest <- apply(cohort[, c("PC1", "PC2")], 1, function(person) {
    distance <- (centres$PC1 - person["PC1"])^2 + (centres$PC2 - person["PC2"])^2
    centres$group[which.min(distance)]
  })

  cat("individuals in the PCA, per group:\n");   print(table(pca$group))
  cat("\nour own individuals, by closest cluster:\n"); print(table(closest))
  cat("\nthe same as a percentage:\n")
  print(round(100 * prop.table(table(closest)), 1))
}


# Two things to check against your own data:
#
# K = 3 suited our cohort, which spreads along a European-African axis with an
# admixed group in between. Look at your own PC1 against PC2 first and count the
# groups you actually see. If your cohort is one homogeneous group, set K <- 1
# and the heterozygosity filter becomes the ordinary cohort-wide 3 SD one.
#
# If your data has no family IDs, PLINK writes the header as '#IID' and there is
# no FID column. Replace c("#FID", "IID") with "#IID" in the merge and in both
# write.table calls.
