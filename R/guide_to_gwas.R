# A guide to genome-wide association analysis and post-analytic interrogation

# https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605

# Last updated 2022-01-17

# Links ----

# https://poissonisfish.com/2017/10/09/genome-wide-association-studies-in-r/
# https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
# https://github.com/monogenea/GWAStutorial
# https://cran.r-project.org/web/packages/statgenGWAS/vignettes/GWAS.html
# https://biometris.github.io/statgenGWAS/index.html
# https://www.nature.com/articles/s43586-021-00056-9
# http://ww7.stat-gen.org/
# https://www.mtholyoke.edu/courses/afoulkes/Data/GWAStutorial/

# Workflow ----


# Ten steps:
# 1) reading data into R to create an R object
# 2) SNP-level filtering (part 1)
# 3) sample-level filtering
# 4) SNP-level filtering (part 2)
# 5) principal component analysis (PCA)
# 6) imputation of non-typed genotypes
# 7) association analysis of typed SNPs
# 8) association analysis of imputed data
# 9) integration of imputed and typed SNP results
# 10) visualisation and quality control of association findings

# Data pre-processing
# .ped file: information on each study participant
# .map file: 1 row for each SNP with corresponding chromosome and coordinate
# .bim file: same information as the .map file + the two observed alleles at each SNP from the .ped file
# .bed file: binary version of the genotype data
# .fam file: participant identification information (same as .ped, without the genotype)
# clinical data file: .txt or .csv, clinical data on each study subject (covariates and phenotypes)

# Install packages ----

# Packages :
# snpStats: read in various formats of genotype data, carry out QC, imputation & association analysis
# SNPRelate: sample-level QC & computationally efficient principal component calculation
# ggplot2, LDheatmap, postgwas: data visualisation
# plyr: data manipulation
# doParallel: parallel processing

# Install required packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("snpStats")
BiocManager::install("SNPRelate")
BiocManager::install("rtracklayer")
BiocManager::install("biomaRt")

install.packages("tidyverse")
install.packages("GenABEL")
install.packages("LDheatmap")
install.packages("doParallel")
install.packages("coin")
install.packages("igraph")
install.packages("https://github.com/merns/postgwas/releases/download/1.11-2/postgwas_1.11-2.zip", repos=NULL)

library(devtools)
install.packages("http:::cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")

# Load packages ----

library(tidyverse)
library(snpStats)
library(SNPRelate)

# Set paths ----

bed_file <- "tutorial/tutorial_files/GWAStutorial.bed"
bim_file <- "tutorial/tutorial_files/GWAStutorial.bim"
fam_file <- "tutorial/tutorial_files/GWAStutorial.fam"
clinical_file <- "tutorial/tutorial_files/GWAStutorial_clinical.csv"

# Specify parameters to be used in the data processing and analysis
#data.dir <- "tutorial/tutorial_files/"
#output.dir <- "tutorial/tutorial_output/"

# Input files
#gwas.fn <- lapply(c(bed = "bed", bim = "bim", fam = "fam", gds = "gds"),
 #                 function(n) sprintf("%s/chr16_1000g_CEU.%s", data.dir, n))
#clinical.fn <- sprintf("%s/GWAStutorial_clinical.csv", data.dir)
#onethou.fn <- lapply(c(info = "info", ped = "ped"),
 #                    function(n) sprintf("%s/chr16_1000g_CEU/%s", data.dir, n))
#protein.coding.coords.fname <- sprintf("%s/ProCodgene_coords.csv", output.dir)

# Output files
#gwaa.fname <- sprintf("%s/GWAStutorialout.txt", output.dir)
#gwaa.unadj.fname <- sprintf("%s/GWAStutorialoutUnadj.txt", output.dir)
#impute.out.fname <- sprintf("%s/GWAStutorial_imputationOut.csv", output.dir)
#CETP.fname <- sprintf("%s/CETP_GWASout.csv", output.dir)

# Step 1 - reading data into R to create an R object ----

# 1.1 - read in PLINK files to create list
geno <- snpStats::read.plink(bed = bed_file,
                             bim = bim_file,
                             fam = fam_file,
                             na.strings = ("-9"))

# 1.2.1 - obtain the genotypes SnpMatrix object from generated list
genotypes <- geno$genotypes
print(genotypes)

# 1.2.2 - obtain the SNP information table from list
snps <- geno$map %>% 
  dplyr::as_tibble() %>% 
  dplyr::select(chr = chromosome,
                snp = snp.name,
                gen_dist = cM,
                position,
                allele_1 = allele.1,
                allele_2 = allele.2)

# Clean global environment
rm(geno, bed_file, bim_file, fam_file)

# 1.3 - read in clinical file
clinical <- readr::read_csv(clinical_file) %>% 
  dplyr::mutate(fam_id = as.character(FamID),
                cad = as.factor(CAD),
                sex = as.factor(sex)) %>% 
  dplyr::select(fam_id, cad, sex:ldl) %>% 
  tibble::column_to_rownames(var = "fam_id")

# 1.4 - subset genotypes for individuals with clinical data
genotypes <- genotypes[rownames(clinical), ]

# Clean global environment
rm(clinical_file)

# Step 2 - SNP-level filtering (part 1) ----

# Typical SNP filtering workflow :
# 1) large amount of missing data
# 2) low variability
# 3) sample-level filtering (Step 3)
# 4) possible genotyping errors

# SNP-level filtering :
# - call rate = for a given SNP, proportion of inviduals in the study for which the corresponding
# SNP information is not missing (95% -> less than 5% of missing data)
# - minor allele frequency (MAF) = a large degree of homogeneity at a given SNP
# across study participants generally results in inedequate power to infer a statistically
# significant relationship between the SNP and the trait under study (here, we remove SNPs
# for which the MAF is less than 1%)

# 2.1 - Create SNP summary statistics (MAF, call rate, ...)
snp_summary <- snpStats::col.summary(genotypes)
print(head(snp_summary))

# 2.2 - set thresholds
cr_threshold <- 0.95
maf_threshold <- 0.01

# 2.3 - filter on MAF and call rate & remove NAs
use <- with(snp_summary, (!is.na(MAF) & MAF > maf_threshold) & Call.rate >= cr_threshold)
use[is.na(use)] <- FALSE
ncol(genotypes) - sum(use)  # 203,287 SNPs will be removed
table(use)  # 658,186 SNPs will be kept / 203,287 SNPs will be removed

# 2.4 - subset genotypes and SNP summary data for SNPs that pas CR and MAF criteria
genotypes <- genotypes[, use]
snp_summary <- snp_summary[use, ]

# Step 3 - sample-level filtering ----

# Criteria for sample-level filtering :
# - missing data
# - sample contamination
# - correlation (for population-based investigations)
# - racial, ethnic, or gender ambiguity or discordance

# Call rate : we exclude individuals who are missing genotype data across
# more and a pre-defined percentage of the typed SNPs (= sample call rate)

# Heterozygosity = presence of each of the two alleles at a given SNP within an 
# individual. This is expected under HWE to occur with probability 2*p*(1-p) where
# p is the dominant allele frequency at that SNP (assuming bi-allelic SNP). Excess
# heterozygosity across typed SNPs within an individual may be an indication of poor sample
# quality, while deficient heterozygosity can indicate inbreeding or other substructure in
# that person. Samples with an inbreeding coefficient |F| = (1 - O/E) > 0.10 are removed,
# where O and E are respectively the observed and expected counts of heterozygous SNPs
# within an individual.

# 3.1.1 - Create sample statistics (Call rate, heterozygosity)
sample_summary <- snpStats::row.summary(genotypes)

# 3.1.2 - Add the F stat (inbreeding coefficient) to sample_summary
maf <- snp_summary$MAF
call_matrix <- !is.na(genotypes)
het_exp <- 