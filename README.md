


##Let's first inspect the raw data, how many variants (and what type of variants) and samples do we have?

bcftools stats *samples* 

# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1019
SN      0       number of records:      864725
SN      0       number of no-ALTs:      144311
SN      0       number of SNPs: 720414
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0


###LENIENT VARIANT FILTER###
#We use PLINK to first do a lenient variant pass (missingness >80%) 
plink2 --bfile *samples* --geno 0.2 --out lenient_filter 

###SAMPLE MISSINGNESS FILTER###
plink2 --bfile lenient_filter --mind 0.02 --out sample_missingness 

###SEX CONCORDANCE###
plink2 --bfile sample_missingness --check-sex 

#Plot simple histogram of the sex (x axis the F statistic and the y axis count) !!NO TITLES OR TEXT JUST THE PLOT 

###PRELIMINARY PCA AND HETEROZYGOSITY###
#1. Compute a within-sample PCA to cluster individuals 
#2. For each cluster compute heterozygosity and remove individuals +- 3 SD away 

###KINSHIP ESTIMATION###
#1. Perform KING in plink2, make a simple histogram/barplot with count of kinship coefficient <0.04 (no relationship), 0.04-0.09 (third degree relatives), 0.09-0.18 (second degree relatives), 0.18-0.35 (first degree), and >0.35 (duplicate sample). !!NO TITLES/CAPTIONS/TEXT, JUST THE PLOT 

###STRICT VARIANT FILTER###
#We use PLINK to first do a strict variant pass (missingness >98%) 
plink2 --bfile sample_missingness --geno 0.02 --out strict_filter 

###HARDY WEINBERG###
plink2 --bfile strict_filter --hwe 1e-6 --out hwe_done

###MAF###
plink2 --bfile strict_filter --maf 0.01 --out maf_done