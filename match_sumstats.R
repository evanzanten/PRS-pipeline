#!/usr/bin/env Rscript
#
# Attach our variant names to a GWAS summary statistics file, and keep the variants
# that are in both.
#
# The GIGASTROKE files give a chromosome, a position and two alleles, but no variant
# name, while every PRS tool wants one. Our own imputed data does have names: the
# rsIDs of the reference panel, which Beagle attached during imputation. So we match
# the two files on chromosome, position and alleles, and take the name from our data.
#
# Usage:  Rscript match_sumstats.R
#
# Change the four settings below to point at your own files.

library(data.table)

PVAR    <- "imputed.pvar"                 # our imputed variants (from section 2.1)
GWAS    <- "GCST90104540_buildGRCh37.tsv.gz"
OUTFILE <- "sumstats_eur.txt"
N_EFF   <- 236506                         # effective N of this GWAS (see section 3.1)

# ---------------------------------------------------------------- our variants --
# The .pvar file starts with a block of ## header lines, which we skip. The columns
# we need are the chromosome, the position, the name, and the two alleles.
ours <- fread(cmd = paste("grep -v '^##'", PVAR))
setnames(ours, c("#CHROM", "POS", "ID", "REF", "ALT"),
               c("chr", "pos", "name", "ref", "alt"), skip_absent = TRUE)

# A few thousand variants have no rsID (they are written as a dot), because the
# reference panel has no name for them either. We give those a name of their own,
# so that every variant can be told apart from every other one.
ours[name == ".", name := paste(chr, pos, ref, alt, sep = ":")]

# One key per variant: chromosome, position and the two alleles sorted alphabetically.
# Sorting them means a variant is recognised whichever way round a file writes its
# alleles (A/G in one file and G/A in the other are the same variant).
ours[, key := paste(chr, pos, pmin(ref, alt), pmax(ref, alt), sep = ":")]

# ------------------------------------------------------------ the GWAS results --
gwas <- fread(GWAS)
gwas[, key := paste(chromosome, base_pair_location,
                    pmin(effect_allele, other_allele),
                    pmax(effect_allele, other_allele), sep = ":")]

# ------------------------------------------------------------------ the match --
m <- merge(gwas, ours[, .(key, name)], by = "key")

# The columns that the PRS tools expect. A1 is the effect allele: the allele that
# the effect size belongs to.
out <- m[, .(SNP = name, CHR = chromosome, BP = base_pair_location,
             A1 = effect_allele, A2 = other_allele,
             BETA = beta, SE = standard_error, P = p_value, N = N_EFF)]
fwrite(out, OUTFILE, sep = " ")

cat(nrow(gwas), "variants in the GWAS,", nrow(out), "of which are also in our data\n")
