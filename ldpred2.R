#!/usr/bin/env Rscript
#
# LDpred2-auto: shrink every variant's effect size using how strongly variants
# correlate with each other, and write one weight per variant.
#
# The correlations come from the reference the LDpred2 authors publish, computed in
# UK Biobank Europeans over the HapMap3+ variants. It is one .rds file per chromosome
# plus map_hm3_plus.rds describing the variants. Our own 870 individuals are far too
# few to estimate the correlations from: every chain returns NA.
#
# Usage:  Rscript ldpred2.R          (submit it, it needs several hours and ~64 GB)
#
# Change the five settings below to point at your own files.

suppressPackageStartupMessages({library(bigsnpr); library(data.table)})

GWAS    <- "sumstats_eur.txt"      # the tidy file from section 3.1
REF_DIR <- "ldpred2_reference"     # where LD_with_blocks_chr*.rds and map_hm3_plus.rds live
TMP_DIR <- "ldpred2_tmp"           # scratch for the on-disk correlation matrix
OUTFILE <- "ldpred2_weights.txt"
NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# bigsnpr refuses to run if the BLAS library is multi-threaded as well
bigparallelr::set_blas_ncores(1)

## 1. the reference variants, and our summary statistics matched to them ##########
# snp_match looks for the columns chr, pos, a0, a1 and beta, so we rename ours to
# those. Note a1 is the effect allele and a0 the other one, which is the opposite of
# what those names suggest in some other tools.
#
# We keep the reference's own ld column: it is the LD score of each variant, which
# step 3 needs and which is already computed for us.

map <- readRDS(file.path(REF_DIR, "map_hm3_plus.rds"))
setDT(map)
cat("reference variants:", nrow(map), "\n")

ss <- fread(GWAS)
setnames(ss, c("SNP","CHR","BP","A1","A2","BETA","SE","P","N"),
             c("rsid","chr","pos","a1","a0","beta","beta_se","p","n_eff"))
info <- snp_match(ss, as.data.frame(map[, .(chr, pos, a0, a1, rsid, ld)]),
                  join_by_pos = FALSE)
setDT(info)
cat("variants matched to the reference:", nrow(info), "\n")

## 2. the correlation matrix ######################################################
# The reference is stored one chromosome at a time and is far too large to hold in
# memory, so we build it on disk (an SFBM) and add one chromosome at a time, keeping
# only the variants we have.
#
# snp_match reports _NUM_ID_, the row of the variant in the whole reference, while
# each chromosome's file is numbered from 1 within that chromosome. match() converts
# between the two.

dir.create(TMP_DIR, showWarnings = FALSE)
tmp <- file.path(TMP_DIR, "corr_ref")

for (ch in 1:22) {
  ind.chr  <- info$`_NUM_ID_`[info$chr == ch]
  if (!length(ind.chr)) next
  ind.ref  <- match(ind.chr, which(map$chr == ch))
  corr_ch  <- readRDS(file.path(REF_DIR, paste0("LD_with_blocks_chr", ch, ".rds")))[ind.ref, ind.ref]
  if (ch == 1) corr <- as_SFBM(corr_ch, tmp, compact = TRUE)
  else         corr$add_columns(corr_ch, nrow(corr))
  cat("chr", ch, "done\n")
}

## 3. heritability, then the model ################################################
# LD score regression reads the heritability off the relationship between a variant's
# chi-square and how much LD it sits in. ld_size is the size of the whole reference,
# not of our subset, because that is what the ld column was computed over.

ldsc <- snp_ldsc(info$ld, ld_size = nrow(map), chi2 = (info$beta / info$beta_se)^2,
                 sample_size = info$n_eff, blocks = NULL)
cat("LD score regression h2:", ldsc[["h2"]], "\n")

auto <- snp_ldpred2_auto(corr, as.data.frame(info), h2_init = max(ldsc[["h2"]], 0.001),
                         vec_p_init = seq_log(1e-4, 0.2, 30), ncores = NCORES,
                         allow_jump_sign = FALSE, shrink_corr = 0.95)

## 4. keep the chains that agree ##################################################
# auto runs 30 independent chains from different starting points. Chains that wander
# off return a heritability of NA, or one far from the rest, and are dropped. The
# weights are the average over what is left. If nothing is left the model has not
# converged and the result must not be used.

h2s <- sapply(auto, function(a) a$h2_est)
ps  <- sapply(auto, function(a) a$p_est)
cat("chains with a usable h2:", sum(is.finite(h2s)), "of", length(auto), "\n")

ok <- which(is.finite(h2s) & is.finite(ps) & h2s > 0)
keep <- if (length(ok) >= 2) {
  med <- median(h2s[ok]); ok[h2s[ok] > 0.7 * med & h2s[ok] < 1.4 * med]
} else ok
cat("chains kept:", length(keep), "| h2 range:",
    paste(round(range(h2s[keep]), 4), collapse = " - "), "\n")
if (!length(keep)) stop("LDpred2-auto did not converge")

# sapply gives a vector rather than a matrix when only one chain survives
bmat <- sapply(auto[keep], function(a) a$beta_est)
beta <- if (is.matrix(bmat)) rowMeans(bmat) else bmat

out <- data.frame(SNP = info$rsid, A1 = info$a1, BETA = beta)
out <- out[is.finite(out$BETA), ]
fwrite(out, OUTFILE, sep = "\t")
cat("weights written:", nrow(out), "\n")
