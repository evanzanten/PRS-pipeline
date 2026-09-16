# Genotype quality control - protocol paper sections 3.1 to 3.7
#
# Runs the QC steps in the order described in the paper: variant and sample
# missingness (3.1), sex concordance (3.2), preliminary PCA and heterozygosity
# (3.3), relatedness (3.4), strict variant missingness (3.5), Hardy-Weinberg
# (3.6) and minor allele frequency (3.7).
#
# Usage:  bash qc.sh
#
# Keep qc.sh and figures.R in the same folder, and run it from that folder:
# the script writes its output into whatever directory you start it in.
#

#Stop at the first command that fails. Without this, one failed step is
#followed by a dozen confusing errors from every step that needed its output.
set -e


####PREPARATION####
#For reproducibility, we assign a seed. This means that in steps where randomization is applied, we will get the same results every time we run the script.
SEED=1

#How many processor cores PLINK may use. Left to itself PLINK takes every core
#it can see, which is fine on your own machine but a bad idea on a shared
#server: it slows the machine down for everyone else, and on a login node it
#will often be refused outright and stop with
#"libgomp: Thread creation failed". Four is plenty for a cohort of this size.
THREADS=4



###Reference Data###
#For 1000 genomes:
#VCF: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#Beagle VCF: https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/
#Pedigree file: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
#FASTA files: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
#Legend file (for Beagle): https://www.chg.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz
#Map file (for Beagle): https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip


###QUALITY CONTROL STEPS###
#Throughout this pipeline, we will be working with SNP array data derived from ~1020 Brazilian individuals, of which we have 513 cases (stroke), and 473 controls.

#We start with a VCF file containing genotypes of all individuals. We also have a phenotype textfile that states whether each individual is a case or a control, and what their age and sex is.

#First make a variable of these two files. These are the only two lines in this
#script that you need to change to run it on your own data:

VCF=genotypes.vcf.gz
PHENO=phenotypes.txt

#The phenotype file is expected to be tab-separated with a header line, and to
#have the sample identifier, the case/control status and the sex somewhere in
#it. Ours looks like this, and the column numbers used further down refer to it:

#  Sample ID   FID                     IID                     Status   Sex      Age
#  00000301    00000301_2-375928.CEL   00000301_2-375928.CEL   Case     FEMALE   60
#  00000305    00000305_2-362470.CEL   00000305_2-362470.CEL   Case     MALE     53

#Our plots are drawn by figures.R, which sits next to this script. It needs to
#know where to find the PLINK output and where to put the plots:

export OUT=.
export FIG=figures
mkdir -p "$FIG"


#Preliminary steps:

  #1. First, we inspect the raw data to see how many, and what types of variants we have.
  bcftools stats "$VCF" > cohort_stats_preQC.txt

  #What does the result look like?:
  cat cohort_stats_preQC.txt

  #For us, the first lines look like this:
  # SN    [2]id   [3]key  [4]value
#  SN      0       number of samples:      1019
#  SN      0       number of records:      864725
#  SN      0       number of no-ALTs:      144311
#  SN      0       number of SNPs: 720414
#  SN      0       number of MNPs: 0
#  SN      0       number of indels:       0
#  SN      0       number of others:       0
#  SN      0       number of multiallelic sites:   0
#  SN      0       number of multiallelic SNP sites:       0

  #So we have 1019 individuals and 864,725 variants, all of them SNPs. The
  #144,311 "no-ALTs" are SNPs at which everybody in our cohort turned out to
  #have the same genotype, so only one allele was ever seen. They carry no
  #information and will disappear at the allele frequency step (3.7).

  #2. We convert the VCF into PLINK's own format, which every step after this
  #one uses. At the same time we attach the reported sex and the case/control
  #status, because these live in the phenotype file and not in the VCF. We need
  #the sex for the concordance check (3.2) and the status for Hardy-Weinberg (3.6).

  #PLINK identifies a sample by two names, a family ID and an individual ID,
  #while a VCF has only one name per sample. --double-id copies the VCF name
  #into both, which is the simplest thing to do when you have no family structure.

  #First we rewrite the phenotype file into the two columns PLINK expects.
  #PLINK reads sex as 1 for male and 2 for female, but also accepts the words,
  #so MALE and FEMALE can be passed through unchanged. Case/control status has
  #to become 2 for a case and 1 for a control, which is the part people get
  #backwards most often. Change the column numbers if your file is laid out
  #differently: below, $2 is FID, $3 is IID, $4 is the status and $5 is the sex.

  awk -F'\t' 'NR==1 {print "#FID\tIID\tSEX\tPHENO"; next}
              {p = ($4=="Case") ? 2 : ($4=="Control") ? 1 : "NA"
               print $2"\t"$3"\t"$5"\t"p}' "$PHENO" > sex_pheno.txt

  plink2 --vcf "$VCF" --double-id \
         --update-sex sex_pheno.txt \
         --pheno sex_pheno.txt --pheno-name PHENO \
         --threads $THREADS --make-bed --out 00_all

  #What does the result look like?: PLINK tells you what it managed to attach.
  #For us:
#  1 binary phenotype loaded (513 cases, 473 controls).
#  --update-sex: 986 samples updated.

  #Note that 986 is less than the 1019 samples in the VCF. The 33 that were not
  #updated are the subject of the next step.


###REMOVE SAMPLES THAT ARE NOT STUDY PARTICIPANTS###
#Genotyping arrays are usually run with reference DNA on every plate as a
#quality control for the lab, and those wells are delivered to you along with
#your own samples. Ours are CEPH1463-02 and the GIAB sample NA24385.
#
#They have to go before anything else, because they are not people in the study
#and they break three later steps: many copies of one DNA look like a large set
#of duplicate pairs in the kinship step (3.4), they form their own tight
#clusters in the preliminary PCA (3.3), and they have no reported sex to check
#(3.2) or case/control status to use (3.6).
#
#The phenotype file is what defines the cohort, so we keep the samples that
#appear in it and report anything genotyped that does not.

awk 'NR>1 {print $1"\t"$2}' sex_pheno.txt > 00_participants.txt

#What is in the genotype data but not in the phenotype file?
awk 'NR==FNR {p[$2]; next} !($2 in p) {print $2}' \
    00_participants.txt 00_all.fam > 00_non_participants.txt

cat 00_non_participants.txt

#For us this lists 33 samples, all of them named after one of the two control
#DNAs. If anything appears here that you do not recognise, stop and find out
#what it is before going on: it means the genotypes and the phenotype file
#disagree about who is in the study.

plink2 --bfile 00_all --keep 00_participants.txt --threads $THREADS --make-bed --out 00_raw

#We now have 986 samples and 864,725 variants.


###MISSINGNESS STATISTICS###
#Missingness is the fraction of genotypes that the array failed to call. It can
#be counted per variant (a probe that works badly in everyone) or per sample (a
#DNA that worked badly at every probe), and both need a threshold.
#
#Rather than take a threshold on faith, we compute the two distributions first
#and look at them. --missing writes one file for each: .vmiss per variant and
#.smiss per sample.

plink2 --bfile 00_raw --missing --threads $THREADS --out 01_missing_raw

#Plot both distributions as histograms. The counts are on a log scale, because
#nearly every variant and sample sits in the first bin and a linear axis would
#hide the tail that the threshold has to be chosen against.
Rscript figures.R missingness


###LENIENT VARIANT FILTER###
#We now do a lenient variant pass, dropping variants missing in more than 20%
#of samples (i.e. keeping those with a call rate above 80%).
#
#Variants are filtered before samples so that a handful of bad probes does not
#make a good sample look bad. The strict variant pass comes later (3.5), once
#the bad samples are out of the way.

plink2 --bfile 00_raw --geno 0.2 --threads $THREADS --make-bed --out 02_lenient_filter

#For us this removed 1,248 variants, leaving 863,477.


###SAMPLE MISSINGNESS FILTER###
#Now the samples. We drop those missing more than 2% of their genotypes, i.e.
#keeping those with a call rate above 98%. A sample below that has usually
#failed in the lab rather than in one particular place in the genome.

plink2 --bfile 02_lenient_filter --mind 0.02 --threads $THREADS --make-bed --out 03_sample_missingness

#For us this removed 8 samples, leaving 978.


###SEX CONCORDANCE###
#We now check that the sex recorded in the phenotype file matches the sex we can
#read off the genotypes. A mismatch usually means a sample was swapped somewhere
#between the clinic and the plate, and a swapped sample carries the wrong
#phenotype, so it has to go.
#
#The check works on the X chromosome. Males have one copy and so cannot be
#heterozygous on it; females have two and are heterozygous at many positions.
#PLINK summarises this as an inbreeding coefficient F, which lands near 1 for
#males and near 0 for females. Below we call F < 0.2 female and F > 0.8 male,
#and treat anything in between as ambiguous.
#
#This is the one step that needs PLINK 1.9: --check-sex does not exist in PLINK 2.

plink --bfile 03_sample_missingness --check-sex 0.2 0.8 --out 04_sexcheck

#What does the result look like?: one line per sample, with the reported sex
#(PEDSEX), the sex read from the genotypes (SNPSEX), and STATUS saying OK or
#PROBLEM. For us the first line is:
#  FID                     IID                     PEDSEX  SNPSEX  STATUS   F
#  00000301_2-375928.CEL   00000301_2-375928.CEL   2       1       PROBLEM  0.9987

#Plot a histogram of F, with the two cutoffs marked. In a clean cohort you see
#two tight groups, one at each end, and almost nothing in the middle.
Rscript figures.R sex

#Remove the samples flagged PROBLEM. The $3!=0 guard skips samples with no
#reported sex, which cannot be checked against anything (there should be none
#left now that the control DNAs have been dropped).
awk 'NR>1 && $5=="PROBLEM" && $3!=0 {print $1"\t"$2}' 04_sexcheck.sexcheck > 04_sex_remove.txt

plink2 --bfile 03_sample_missingness --remove 04_sex_remove.txt --threads $THREADS --make-bed --out 05_sex_ok

#For us this removed 50 samples, leaving 928. That is about 5%, which is more
#than you would like to see: 32 samples reported male came out female and 14
#reported female came out male. A rate this symmetric points at labelling rather
#than at the DNA, so it is worth asking whoever prepared the plates before you
#accept it. We remove them either way, because we cannot tell which of the two
#records is the wrong one.


###PRELIMINARY PCA AND HETEROZYGOSITY###
#Heterozygosity is the fraction of a person's genotypes that carry two different
#alleles. Too high suggests the sample is a mixture of two DNAs; too low
#suggests inbreeding or a degraded sample. Either way it is a sign of a sample
#we should not trust.
#
#The catch is that the normal amount of heterozygosity differs between ancestry
#groups, so in a cohort like ours, which is admixed, comparing everybody to one
#cohort-wide average would flag whole ancestry groups instead of bad samples.
#So we do it in two parts:
#
#  1. Compute a PCA to split the cohort into broad ancestry clusters
#  2. Within each cluster, remove individuals more than 3 SD from that
#     cluster's own mean heterozygosity
#
#PCA needs variants that are roughly independent of one another, so we first
#prune for linkage disequilibrium. We also drop the long-range LD regions of
#Price et al. (2008), a short list of places in the genome where correlations
#stretch so far that they dominate the top components and swamp the ancestry
#signal we are after. Those coordinates are on build GRCh37, so check the build
#of your own data before using them.

cat > high_ld_b37.txt <<'REGIONS'
1	48000000	52000000	highLD1
2	86000000	100500000	highLD2
2	134500000	138000000	highLD3
2	183000000	190000000	highLD4
3	47500000	50000000	highLD5
3	83500000	87000000	highLD6
3	89000000	97500000	highLD7
5	44500000	50500000	highLD8
5	98000000	100500000	highLD9
5	129000000	132000000	highLD10
5	135500000	138500000	highLD11
6	25000000	35000000	highLD12
6	57000000	64000000	highLD13
6	140000000	142500000	highLD14
7	55000000	66000000	highLD15
8	7000000	13000000	highLD16
8	43000000	50000000	highLD17
8	112000000	115000000	highLD18
10	37000000	43000000	highLD19
11	46000000	57000000	highLD20
11	87500000	90500000	highLD21
12	33000000	40000000	highLD22
12	109500000	112000000	highLD23
20	32000000	34500000	highLD24
REGIONS

#--indep-pairwise 200 50 0.2 slides a 200-variant window along the genome in
#steps of 50, and within each window drops one of every pair correlated at
#r2 > 0.2. We also restrict to the autosomes and to common variants.

plink2 --bfile 05_sex_ok --autosome --maf 0.05 \
       --exclude bed1 high_ld_b37.txt \
       --indep-pairwise 200 50 0.2 --threads $THREADS --out 06_prune

#This left us with 109,774 variants out of 863,477, written to 06_prune.prune.in.

#Only the first two components are needed here. We want a coarse split into
#broad groups, not an ancestry assignment (that is section 3.9).
#
#'approx' is PLINK's randomised algorithm. PLINK will warn you that it is only
#recommended above 5000 samples, and we have fewer, but the exact algorithm
#takes over ten minutes on a cohort this size while 'approx' takes seconds, and
#the difference between them is far smaller than the width of the clusters we
#are about to draw. Being randomised is why we set a seed at the top.
plink2 --bfile 05_sex_ok --extract 06_prune.prune.in \
       --pca 2 approx --seed $SEED --threads $THREADS --out 06_pca

plink2 --bfile 05_sex_ok --extract 06_prune.prune.in --het --threads $THREADS --out 06_het

#figures.R does the clustering, draws the PCA and the heterozygosity plots, and
#writes the list of outliers to 06_het_outliers.txt.
Rscript figures.R pca_het

#For us it found clusters of 626, 243 and 59 individuals, and flagged 14 outliers.

plink2 --bfile 05_sex_ok --remove 06_het_outliers.txt --threads $THREADS --make-bed --out 07_het_ok

#Leaving 914 samples.


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
