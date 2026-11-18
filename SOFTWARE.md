 # Software

  Every tool the pipeline calls, at the version it was run with.
  `container/pipeline.def` is the machine-readable form of this table, and the
  two must be changed together.

  The container installs these same versions with two exceptions, both noted in
  the table. PLINK 2: the build the results were produced with is no longer
  distributed, so the image pins a current one. BridgePRS: it is not in the image
  at all, so section 3.7 needs it fetched and installed separately. This is why
  the archived image is the reproducible artifact.

  | Tool | Version | What it does here | Sections |
  |------|---------|-------------------|----------|
  | PLINK | 1.9.0-b.7.11 (2025-08-19) | The three things PLINK 2 does not do: `--check-sex`, `--test-missing`, and `--flip` with `--freq --keep-allele-order` | 1.2, 1.5, 1.8 |
  | PLINK 2 | 2.00a3.7 (2022-10-24) | Everything else: filtering, frequencies, Hardy-Weinberg, kinship, PCA, format conversion, scoring. This is the build every number in the README came from | 1.1-3.8 |
  | PLINK 2 | v2.0.0-b.1 (2026-10-01) | What the container installs. PLINK 2 removes old dated builds from its download server, so 2.00a3.7 can no longer be fetched and the image pins the current build instead | 1.1-3.8 |
  | bcftools / HTSlib | 1.18 | Reading and filtering VCF files, and indexing them | inspection, 2.x |
  | Beagle | 5.5 (27Feb25.75f) | Phasing and imputation in one step | 2.x |
  | bref3 | 28Jun21.220 | Converts a reference panel to Beagle's compressed haplotype format | 2.x |
  | UCSC liftOver | undated binary | Only if your data is on GRCh38 and you want GRCh37 | 1.8 |
  | Perl | 5.x, from the base image | Runs Rayner's strand-checking script, which you fetch yourself | 1.8 |
  | PRSice-2 | 2.3.5 (2021-09-20) | Clumping and thresholding, and a score per p-value threshold | 3.3 |
  | PRS-CS | development version of 2024-11-21 (commit 5330390), 11 commits after release v1.1.0 | Continuous-shrinkage weights from one GWAS | 3.5 |
  | PRS-CSx | development version of 2025-11-12 (commit 37eaeb3), 15 commits after release v1.1.0 | Continuous-shrinkage weights from several GWAS jointly | 3.6 |
  | BridgePRS | development version of 2026-02-23 (commit 136ca4c). The repository's only release, 0.1, dates from 2023-09-28, so there is no release number for the code we ran. **Not in the container** | Cross-ancestry weights, fitting the European arm first and using it as the prior for the African arm | 3.7 |
  | R: bigsnpr | 1.12.2, from the 2024-01-15 snapshot | LDpred2-auto, and the optional `snp_autoSVD` route for long-range LD | 1.3, 3.4 |
  | R: bigstatsr | 1.5.12, from the same snapshot | The on-disk matrix code bigsnpr is built on | 1.3, 3.4 |
  | Python: numpy, scipy, h5py | from the base image | Required by PRS-CS and PRS-CSx | 3.5, 3.6 |
  | R | 4.3.2 | The summary-statistics join, the ancestry clustering, LDpred2 and the validation | throughout |
  | R: data.table | 1.14.10 | Reading and joining large variant tables | throughout |
  | R: ggplot2 | 3.4.4 | Drawing the figures, with `grid` and `gtable` from the same install. `pca_clusters.R` and `match_sumstats.R` need none of these: they use base R and data.table only | `figures.R` |
  | R: scales | 1.3.0 | Axis formatting, same as above | `figures.R` |

  R packages are installed from a dated CRAN snapshot (2024-01-15) so that a
  rebuild gets these versions. HTSlib and bcftools are built from source,
  because the distribution's packages are several versions behind and would put
  a different bcftools in the image from the one named here.

  ## Why two versions of PLINK

  PLINK 2 not an update but more of a rewrite, and a few things in PLINK 1.9 were not
  carried across. This pipeline uses PLINK 2 everywhere it can, and PLINK 1.9
  only for `--check-sex` in 1.2, `--test-missing` in 1.5, and `--flip` and
  `--freq --keep-allele-order` in 1.8.

  There is a trap in moving between them. PLINK 1.9 writes the sex chromosomes as
  numbers, where 23 is X, 24 is Y, 25 is the pseudo-autosomal region and 26 is the
  mitochondrion, while PLINK 2 writes them as `X`, `Y`, `XY` and `MT`. A test
  written for one convention silently returns zero on the other.

  ## What is not in the container

  The data hereunder is not included since they are tens of GB and change often: 
  
  | What | Size | Used by |
  |------|------|---------|
  | 1000 Genomes Phase 3 VCFs, GRCh37 | ~15 GB | 1.9 |
  | 1000 Genomes in Beagle bref3 format | ~3 GB | 2.x |
  | LDpred2 HapMap3+ LD reference | 14 GB, 29 GB unpacked | 3.4 |
  | PRS-CS 1000 Genomes LD blocks, European | 4.6 GB | 3.5, 3.6 |
  | PRS-CS 1000 Genomes LD blocks, African | 7.4 GB | 3.6 |
  | GWAS summary statistics, per ancestry | ~1 GB each | 3.1 onwards |

  
  And one script: `HRC-1000G-check-bim.pl` (version 4.3.3 here) from Will
  Rayner's tools page. It is a single Perl file on a personal page whose download
  URL changes between versions, so we do not bundle it. Perl is in the image, so
  fetch the script and it runs.

