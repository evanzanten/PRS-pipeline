 # Software

  Every tool the pipeline calls, at the version it was run with. These are the
  same versions the container installs: `container/pipeline.def` is the
  machine-readable form of this table, and the two must be changed together.

  | Tool | Version | What it does here | Sections |
  |------|---------|-------------------|----------|
  | PLINK | 1.90b7.11 (2023-12-11) | The three things PLINK 2 does not do: `--check-sex`, `--test-missing`, and `--flip` with `--freq --keep-allele-order` | 1.2, 1.5, 1.8 |
  | PLINK 2 | 2.00a3.7 (2022-10-24) | Everything else: filtering, frequencies, Hardy-Weinberg, kinship, PCA, format conversion, scoring | 1.1-3.8 |
  | bcftools / HTSlib | 1.18 | Reading and filtering VCF files, and indexing them | inspection, 2.x |
  | Beagle | 5.5 (27Feb25.75f) | Phasing and imputation in one step | 2.x |
  | bref3 | 28Jun21.220 | Converts a reference panel to Beagle's compressed haplotype format | 2.x |
  | UCSC liftOver | undated binary | Only if your data is on GRCh38 and you want GRCh37 | 1.8 |
  | Perl | 5.x, from the base image | Runs Rayner's strand-checking script, which you fetch yourself | 1.8 |
  | PRSice-2 | 2.3.5 (2021-09-20) | Clumping and thresholding, and a score per p-value threshold | 3.3 |
  | R: bigsnpr | from the 2024-01-15 snapshot | LDpred2-auto | 3.4 |
  | R: bigstatsr | from the same snapshot | The on-disk matrix code bigsnpr is built on | 3.4 |
  | Python: numpy, scipy, h5py | from the base image | Required by PRS-CS and PRS-CSx | 3.5, 3.6 |
  | R | 4.3.2 | The summary-statistics join, LDpred2, the validation, and every figure | throughout |
  | R: data.table | 1.14.10 | Reading and joining large variant tables | throughout |
  | R: ggplot2 | 3.4.4 | Figures | figures.R |
  | R: scales | 1.3.0 | Axis formatting | figures.R |

  R packages are installed from a dated CRAN snapshot (2024-01-15) so that a
  rebuild gets these versions rather than whatever is current. HTSlib and
  bcftools are built from source rather than from the distribution's packages,
  which are several versions behind and would put a different bcftools in the
  image from the one named here.

  ## Why two versions of PLINK

  PLINK 2 is a rewrite, not an update, and a few things in PLINK 1.9 were not
  carried across. This pipeline uses PLINK 2 everywhere it can, and PLINK 1.9
  only for `--check-sex` in 1.2, `--test-missing` in 1.5, and `--flip` and
  `--freq --keep-allele-order` in 1.8.

  There is a trap in moving between them. PLINK 1.9 writes the sex chromosomes as
  numbers, where 23 is X, 24 is Y, 25 is the pseudo-autosomal region and 26 is the
  mitochondrion, while PLINK 2 writes them as `X`, `Y`, `XY` and `MT`. A test
  written for one convention silently returns zero on the other.

  ## What is not in the container

  Two things, both because they are data rather than software, are tens of
  gigabytes, and change on their own schedule:

  | What | Size | Used by |
  |------|------|---------|
  | 1000 Genomes Phase 3 VCFs, GRCh37 | ~15 GB | 1.9 |
  | 1000 Genomes in Beagle bref3 format | ~3 GB | 2.x |
  | LDpred2 HapMap3+ LD reference | 14 GB, 29 GB unpacked | 3.4 |
  | PRS-CS 1000 Genomes LD blocks, European | 4.6 GB | 3.5, 3.6 |
  | PRS-CS 1000 Genomes LD blocks, African | 7.4 GB | 3.6 |
  | GWAS summary statistics, per ancestry | ~1 GB each | 3.1 onwards |
