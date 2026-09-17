# **From genotype data to polygenic risk scores: a practical end-to-end guide for researchers without bioinformatics training**

**Author: van Zanten, E.S.**


**This GitHub is created alongside our paper on PRS computation for non-bioinformaticians. We go through the pipeline in the order described in the paper, starting with some preparation steps, followed by QC and PRS computation. Note that for generating the figures, we have a separate R script in this repository called figures.R.**
**Contents**
**1.1 Variant and sample missingness**  
**1.2 Sex concordance**   
**1.3 Preliminary PCA and heterozygosity**  
**1.4 Relatedness**  
**1.5 Strict variant filtering**   
**1.6 Hardy-Weinberg**  
**1.7 MAF**  
  

### Preparation 

1. **Software, data and programming languages**  
   Make sure you have the software needed for this pipeline downloaded, see [insert singularity container and point to table]. The coding itself is done in bash, and we make the figures in R [version]. For several of the QC and PRS steps, we need to download genetic data of a reference panel. To choose the right reference panel, we refer to our paper section 4.1.

2. **Input data**  
   This pipeline is built for SNP array data, although some of the steps are also relevant for whole exome sequencing (WES) and whole genome sequencing (WGS). For full WES and WGS pipelines, we recommend using [insert good QC papers]. In our case, we use SNP array data in a variant call format (VCF) file. Our data is derived from a stroke GWAS in 986 Brazilian individuals (513 cases, 473 controls).

3. **Configuration**  
   In many scripts, values are often set multiple times throughout the script, and it is often easier and cleaner to assign these values at the beginning of the script so that if you change the values, you only have to do that once at the beginning of the script. Throughout this pipeline, we assign multiple variables to values:  
   *  For reproducibility, we assign a seed. This means that in steps where randomization is applied, we will get the same results every time we run the script.  
   *  When working on computing clusters, it is advisable to set the number of CPU cores that PLINK may use, since if you do not do this, PLINK will attempt to use all available cores on a node, which may exceed your allocated resources. As a default, we will use 4 threads.  
   *  We also specify the output path that will be used to put all the results in, and the input VCF and metadata (for us, the phenotype textfile that states whether each individual is a case or a control, and what their age and sex is)  

```
SEED=1
THREADS=4
OUT=/path_to_output/
VCF=/path_to_vcf/genotypes.vcf.gz
PHENO=/path_to_metadata/phenotypes.txt
```


### Inspecting your data  

Let's start by inspecting our VCF file and the phenotype data, just to get an idea of what each file looks like. Since we specified the paths to those files above, we just have to reference to them in this step.  

1. First inspect the phenotype data. We do not want to print the entire file, just the first two lines of the beginning of the file is enough (head -n 2 does this)   
```
cat "$PHENO" | head -n 2
```  
|Sample ID|FID|IID|Status|Sex|Age
|---------|---|---|------|---|---|
|00000301 |00000301_2-375928.CEL|00000301_2-375928.CEL|Case|FEMALE|60|
|00000305 |00000305_2-362470.CEL|00000305_2-362470.CEL|Case|MALE  |53|


In our case, the FID (family ID) and IID (individual ID) are the same, since we do not have families in our cohort to our knowledge (spoiler: later in the pipeline, we will find out we do have relatives). The first two individuals are both cases, one is male (aged 53) and one is female (aged 60).  

2. Let's also inspect what the VCF looks like, and what types of variants our VCF contains.  
A VCF is often gzipped (compressed; ending in .gz) since it is very large, and it has a long header that contains specifics on the genotypes (e.g. build, how the VCF was derived). Since we just want to inspect the data itself, we therefore have to use a different command than for the phenotyping file, and BCFtools is designed to do this.  

```
bcftools view -H "$VCF" | head -n 1
```
This skips the header (-H), and shows the data for the first variant. For readability, we just paste the genotypes of the first 4 individuals:
|CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|           |
|-----|---|--|---|---|----|------|----|------|-----------|
1     |  86028  | AX-13216142   |  T  |     C    |   .   |    .     |  PR        |      GT  |    0/0     0/0     0/0     0/0 |             

Per column:  
1: chromosome  
86028: position on the chromosome  
AX-13216142: variant ID (in our case, the Axiom assay probe name, not an rsID)  
T: reference allele  
C: alternative allele  
.: quality score (empty, since SNP arrays don't give one)  
.: filter status (empty, no filters applied)  
PR: extra information: here a flag that the reference allele is provisional, we will check this later on in the pipeline  
GT: format of the sample columns: they hold genotypes (GT)  

Genotypes are then counted as follows:  
0/0: two reference alleles (T/T)  
0/1: one reference, one alternative allele (T/C)  
1/1: two alternative alleles (C/C)  
./.: missing, the array failed to call it  

Now we check what kinds of variants, and how many, we have.  
```
bcftools stats "$VCF"
```
This prints a lot of interesting data. For us, the first table is most relevant since it is a summary:  
|SN|[2]id|[3]key|[4]value|  
|--|-----|------|--------|
SN |     0  |     number of samples:   |   986  |
SN |    0    |   number of records:     | 864725 | 
SN  |    0    |   number of no-ALTs:     | 144311 | 
SN   |   0     |number of SNPs: |720414  |
SN     | 0     |  number of MNPs:| 0  |
SN      |0      | number of indels:|       0|  
SN      |0       |number of others: |      0 | 
SN      |0       |number of multiallelic sites:|   0|  
SN      |0       |number of multiallelic SNP sites:  |     0  |

So we have 986 individuals and 864,725 variants, all of them SNPs. The 144,311 "no-ALTs" are SNPs at which everybody in our cohort turned out to have the same genotype, so only one allele was ever seen.

3. We now rewrite the phenotype file into the two columns expected by PLINK. PLINK reads sex as 1 for male and 2 for female, but also accepts the words, so MALE and FEMALE can be passed through unchanged. Case/control status has to become 2 for a case and 1 for a control, which is the part people get backwards most often. Change the column numbers if your file is laid out differently: below, $2 is FID, $3 is IID, $4 is the status and $5 is the sex.
```
awk -F'\t' 'NR==1 {print "#FID\tIID\tSEX\tPHENO"; next}
            {p = ($4=="Case") ? 2 : ($4=="Control") ? 1 : "NA"
             print $2"\t"$3"\t"$5"\t"p}' "$PHENO" > sex_pheno.txt


Now we can use PLINK. Remember to specify the threads and the output directory we put in our configuration! 
plink2 --vcf "$VCF" --double-id \
       --update-sex sex_pheno.txt \
       --pheno sex_pheno.txt --pheno-name PHENO \
       --threads "$THREADS" --make-bed --out "$OUT/first_step"
Here, we first refer to our VCF again, then tell PLINK that the IID and FID are the same. We then update the sex and phenotype status and specify the number of threads PLINK may use. We then generate PLINK files, including a BED, BIM and FAM file (see [insert link] for details on those filetypes). Here, we refer to our output directory for the first time. 

```

What does the result look like?: PLINK tells you what it managed to attach.
For us:
```
1 binary phenotype loaded (513 cases, 473 controls).
--update-sex: 986 samples updated.
```

### 1.1 Missingness
Missingness is the fraction of genotypes that the array failed to call. It can be counted per variant (a probe that works badly in everyone) or per sample (an error in someones DNA that worked badly for every probe). Both variant and sample missingness need a threshold, and the (default) thresholds used in this pipeline are as follows: 
|Step|Threshold|Removes if missing in more than|Keeps call rate of at least|
|------|----|----|---|
|Variant filter (lenient)|0.2|20% of samples|80%|
|Sample filter|0.02|2% of variants|98%|
|Variant filter (strict)|0.02|2% of samples|98%|  

At this stage, we will first perform the lenient variant filtering and the sample filtering. The strict variant filter will be applied later in the pipeline. 
PLINK applies --mind before --geno when both are given in one command, and we want to do the lenient variant filtering before sample filtering (see paper for detailed explanation on this). 

```
plink2 --bfile "$OUT/first_step" --geno 0.2 --threads "$THREADS" --out "$OUT/lenient_var"
In our data, this excludes 1248 variants

plink2 --bfile "$OUT/lenient_var" --mind 0.02 --threads "$THREADS" --out "$OUT/sample_missingn"
In our data, this excludes 8 samples (individuals)

We are now left with 978 samples and 863477 variants. 

``` 


### 1.2 Sex concordance 
We now check that the sex recorded in the phenotype file matches the sex we can read off the genotypes. A mismatch usually means a sample was swapped somewhere between the clinic and the plate, and a swapped sample carries the wrong phenotype, so it has to go. The check works on the X chromosome. Males have one copy and so cannot be heterozygous on it; females have two and are heterozygous at many positions.
PLINK summarises this as an inbreeding coefficient F, which lands near 1 for males and near 0 for females. Below we call F < 0.2 female and F > 0.8 male, and treat anything in between as ambiguous.
Note that this step can only be done in PLINK 1.9 (not PLINK 2, where --check-sex does not exist). 
```
plink --bfile "$OUT/sample_missingn" --check-sex 0.2 0.8 --out "$OUT/sexcheck"

```
What does the result look like? 
One line per sample, with the reported sex (PEDSEX), the sex read from the genotypes (SNPSEX), and STATUS saying OK or
PROBLEM. For us the first line is:
|FID|                     IID |                    PEDSEX|  SNPSEX|  STATUS  | F|
|---|---|---|---|---|---|
00000301_2-375928.CEL |  00000301_2-375928.CEL  | 2   |    1   |    PROBLEM | 0.9987

Plot a histogram of F, with the two cutoffs marked. In a clean cohort you see two tight groups, one at each end, and almost nothing in the middle.  
<img src="figures/sex_fstat.png" alt="Sex check F statistic" width="400">

Remove the samples flagged PROBLEM. For us, this removed 50 samples, leaving 928. That is about 5%, which is more than you would like to see: 32 samples reported male came out female and 14 reported female came out male. A rate this symmetric points at labelling rather than at the DNA, so it is worth asking whoever prepared the plates before you accept it. We remove them either way, just to be absolutely sure our data is clean for the next steps. 
```
We first use awk to create a textfile with the samples that need to be removed. NR>1 is used to skip the first line (the header), then we filter for those rows where column 5 contains "PROBLEM" ($5=="PROBLEM") and print column 1 (FID) and column 2 (IID), separated by a tab ({print $1"\t"$2})
awk 'NR>1 && $5=="PROBLEM" {print $1"\t"$2}' "$OUT/sexcheck.sexcheck" > "$OUT/to_be_removed_sex.txt"

We then tell plink that it should remove those samples that we just generated the textfile for, from our data
plink2 --bfile "$OUT/sample_missingness" --remove "$OUT/to_be_removed_sex.txt" --threads $THREADS --make-bed --out "$OUT/sex_check_finished"

Since in our data, we removed 50 samples, we end up with 928 individuals. 
```

### 1.3 Preliminary PCA and heterozygosity 
Heterozygosity is the fraction of a person's genotypes that carry two different alleles. Too high suggests the sample is a mixture of two DNAs, and too low suggests inbreeding or a degraded sample. Either way it is a sign of a sample we should not trust. However, the 'normal' amount of heterozygosity differs enormously between ancestry groups, so in an admixed cohort like ours, comparing everybody to the mean heterozygosity rates of the entire cohort would flag whole ancestry groups instead of just bad samples. This is why, for cross-ancestry or admixed cohorts, it is advisable to perform a PCA to split the cohort into broad ancestry clusters, then within each cluster we remove individuals more than 3 standard deviations away from that cluster's mean heterozygosity. For relatively homogeneous cohorts such as European cohorts, you can skip the PCA step and directly compute the heterozygosity rates for the entire cohort. 

**Principal Component Analysis**  
We perform a PCA to do a rough separation of our cohort into clusters (see the paper to get a full description of how this step works). PCA needs variants that are roughly independent of one another, so we first prune for linkage disequilibrium. We also drop the long-range LD regions of Price et al. (2008), a short list of places in the genome where correlations stretch so far that they dominate the top components and swamp the ancestry signal we are after. Those coordinates are on build GRCh37, so check the build of your own data before using them. See the file high_ld_b37.txt to get a list of these regions. 

```
In the following command, we use plink to:  
- Prune (drop) one of every correlated variant pairs at r2 > 0.2, using a sliding window the size of 200 variants, and steps of 50 variants. We restrict our pruning to the autosomes and common (MAF > 0.05) variants  
- Remove the regions that we specified in the high_ld_b37.txt file, from the Price et al. (2008) paper.   

plink2 --bfile "$OUT/sex_check_finished" --autosome --maf 0.05 \
       --exclude bed1 high_ld_b37.txt \
       --indep-pairwise 200 50 0.2 --threads $THREADS --out pruned_longrange_ld

This left us with 109,774 variants out of 863,477, written to pruned_longrange_ld.prune.in.
```

We use only the first two principal components, since we want a coarse split into groups, not an ancestry assignment (we leave that for section 3.9). We use PLINK's randomised PCA algorithm (approx). PLINK recommends it only above 5,000 samples and we have fewer, but it is much faster than the exact algorithm, and the difference between the two is far smaller than the spread of the clusters, so it doesn't matter here. Because this algorithm is randomised, we use the seed set in the configuration step. We then calculate heterozygosity for every sample in one go. Each sample's heterozygosity depends only on its own genotypes, so it does not matter that we calculate it for the whole cohort at once. The clusters are used in the next step, where we compare each sample with the average of its own cluster rather than with the whole cohort.


```
plink2 --bfile "$OUT/sex_check_finished" --extract pruned_longrange_ld.prune.in \
       --pca 2 approx --seed $SEED --threads $THREADS --out rough_pca

plink2 --bfile "$OUT/sex_check_finished" --extract pruned_longrange_ld.prune.in --het --threads $THREADS --out full_cohort_het

#figures.R does the clustering, draws the PCA and the heterozygosity plots, and writes the list of outliers to het_outliers.txt.

#For us it found clusters of 626, 243 and 59 individuals, and flagged 14 outliers.

plink2 --bfile "$OUT/sex_check_finished" --remove het_outliers.txt --threads $THREADS --make-bed --out het_finished

#Leaving 914 samples.
```
  <table>
    <tr>
      <td><img src="figures/pca_clusters.png" alt="PCA clusters"
  width="100%"></td>
      <td><img src="figures/heterozygosity.png" alt="Heterozygosity rate per cluster. Solid line: cluster mean; dashed lines: ±3SD"
  width="100%"></td>
    </tr>
    <tr>
      <td align="center"><em>PC1 and PC2, coloured by cluster</em></td>
      <td align="center"><em>Heterozygosity rate per cluster</em></td>
    </tr>
  </table>



###KINSHIP ESTIMATION###
#Related individuals break the assumption that every sample is an independent
#observation, which both the association testing and the PRS evaluation rely on.
#Duplicates are worse still: the same person counted twice.
#
#We estimate kinship with the KING algorithm, which is robust to the population
#structure we just saw in the PCA. It runs on the same pruned variants.
#
#A kinship coefficient is roughly the probability that an allele drawn from one
#person and one from the other are inherited from the same ancestor. The
#conventional bands are:
#
#  < 0.04       unrelated
#  0.04 - 0.09  third-degree relatives (first cousins)
#  0.09 - 0.18  second-degree (half-siblings, grandparent-grandchild)
#  0.18 - 0.35  first-degree (parent-child, full siblings)
#  > 0.35       the same person twice

plink2 --bfile 07_het_ok --extract 06_prune.prune.in --make-king-table --threads $THREADS --out 08_king

#Plot the counts in each band. Note this reports every possible pair, so the
#unrelated bin is enormous and the plot uses a log scale.
Rscript figures.R kinship

#For us, out of the 417,241 pairs among 914 samples:
#  <0.04       417149
#  0.04-0.09       37
#  0.09-0.18       14
#  0.18-0.35       33
#  >0.35            8
#
#The 8 pairs above 0.35 are the four people who were genotyped twice, which is
#worth knowing about before you find them here.
#
#This is also where the previous step earns its place. Run on the data before
#the heterozygosity filter, the 0.04-0.09 band held 1,531 pairs instead of 37:
#a contaminated sample looks slightly related to everybody, so 14 bad samples
#produced around 1,500 spurious relationships between themselves and the rest.

#Remove one of each pair closer than second-degree. --king-cutoff works out for
#itself which member of each pair to drop so that the fewest samples are lost.
plink2 --bfile 07_het_ok --extract 06_prune.prune.in \
       --king-cutoff 0.09 --threads $THREADS --out 08_king_cutoff

plink2 --bfile 07_het_ok --remove 08_king_cutoff.king.cutoff.out.id \
       --threads $THREADS --make-bed --out 09_unrelated

#For us this removed 44 samples, leaving 870. Note that the file PLINK writes
#has a header line, so it has 45 lines for 44 samples - --remove knows this, but
#do not be caught out if you count them yourself.


###STRICT VARIANT FILTER###
#Now that the bad samples are gone, we repeat the variant pass at the threshold
#we actually want: a 98% call rate. Doing it now rather than at the start means
#a variant is not condemned on the basis of samples that have since been removed.

plink2 --bfile 09_unrelated --geno 0.02 --threads $THREADS --make-bed --out 10_strict_filter

#For us this removed 9,026 variants, leaving 854,451.


###HARDY WEINBERG###
#Hardy-Weinberg equilibrium is the genotype frequency you expect from the allele
#frequency if mating is random. A variant that departs from it sharply is
#usually a genotyping error, so it is a standard filter.
#
#This is a case/control cohort, so we test in the controls only. A real risk
#variant is expected to depart from equilibrium in cases, which is the whole
#point of it, and testing there would throw out exactly the variants we are
#looking for. We collect the variants that pass in controls, then apply that
#list to everyone.

awk '$6==1 {print $1"\t"$2}' 10_strict_filter.fam > 11_controls.txt

plink2 --bfile 10_strict_filter --keep 11_controls.txt \
       --hwe 1e-6 --write-snplist --threads $THREADS --out 11_hwe_pass

plink2 --bfile 10_strict_filter --extract 11_hwe_pass.snplist --threads $THREADS --make-bed --out 11_hwe

#For us this removed 409 variants, leaving 854,042.


###MAF###
#Finally we drop variants with a minor allele frequency below 1%. At this sample
#size a rarer variant is carried by too few people to estimate an effect for, and
#it is also where genotyping errors concentrate.
#
#This is also where the 144,311 no-ALT variants we saw at the very beginning
#leave the dataset, since a variant with only one allele has a frequency of zero.

plink2 --bfile 11_hwe --maf 0.01 --threads $THREADS --make-bed --out 12_maf

#For us this removed 403,302 variants, leaving 450,740.


###SUMMARY###
#Samples and variants left after each step.

printf '\n%-28s %10s %10s\n' step samples variants
for step in 00_all 00_raw 02_lenient_filter 03_sample_missingness 05_sex_ok \
            07_het_ok 09_unrelated 10_strict_filter 11_hwe 12_maf; do
    printf '%-28s %10s %10s\n' "$step" "$(wc -l < $step.fam)" "$(wc -l < $step.bim)"
done

#For us:
#  step                            samples   variants
#  00_all                             1019     864725
#  00_raw                              986     864725
#  02_lenient_filter                   986     863477
#  03_sample_missingness               978     863477
#  05_sex_ok                           928     863477
#  07_het_ok                           914     863477
#  09_unrelated                        870     863477
#  10_strict_filter                    870     854451
#  11_hwe                              870     854042
#  12_maf                              870     450740

echo
echo "QC'd data: 12_maf"
echo "Figures:   $FIG"
