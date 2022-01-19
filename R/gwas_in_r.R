# Genome-wide association studies in R

# https://poissonisfish.com/2017/10/09/genome-wide-association-studies-in-r/


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
#library(SNPRelate)

# Set paths ----

india_files <- paste("data/public/Genomics/105Indian_2527458snps",
                     c(".bed", ".bim", ".fam"), sep = "")

malaysia_files <- paste("data/public/Genomics/108Malay_2527458snps",
                        c(".bed", ".bim", ".fam"), sep = "")

china_files <- paste("data/public/Genomics/110Chinese_2527458snps",
                        c(".bed", ".bim", ".fam"), sep = "")

# Read files ----

snps_india <- snpStats::read.plink(bed = india_files[1],
                                   bim = india_files[2],
                                   fam = india_files[3])

snps_malaysia <- snpStats::read.plink(bed = malaysia_files[1],
                                      bim = malaysia_files[2],
                                      fam = malaysia_files[3])

snps_china <- snpStats::read.plink(bed = china_files[1],
                                   bim = china_files[2],
                                   fam = china_files[3])

# Read files ----

china_snps <- snpStats::read.plink(bed = china_files[1],
                                   bim = china_files[2],
                                   fam = china_files[3])

india_snps <- snpStats::read.plink(bed = india_files[1],
                                   bim = india_files[2],
                                   fam = india_files[3])

malaysia_snps <- snpStats::read.plink(bed = malaysia_files[1],
                                      bim = malaysia_files[2],
                                      fam = malaysia_files[3])

rm(china_files, india_files, malaysia_files)
all.equal(snps_china$genotypes)

load("data/conversionTable.RData")

# Merge the three SNP datasets ----

snps <- snps_malaysia
snps$genotypes <- rbind(snps_malaysia$genotypes, snps_india$genotypes, snps_china$genotypes)
colnames(snps$map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
snps$fam <- rbind(snps_malaysia$fam, snps_india$fam, snps_china$fam)

# Rename SNPs present in the conversion table into rs IDs ----

mapped_snps <- intersect(snps$map$SNP, names(conversionTable))
new_ids <- conversionTable[match(snps$map$SNP[snps$map$SNP %in% mapped_snps],
                                 names(conversionTable))]
snps$map$SNP[rownames(snps$map) %in% mapped_snps] <- new_ids
