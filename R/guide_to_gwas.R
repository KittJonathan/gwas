# Genome-wide association studies in R

# Last updated 2022-01-14

# Links ----

# https://poissonisfish.com/2017/10/09/genome-wide-association-studies-in-r/
# https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
# https://github.com/monogenea/GWAStutorial
# https://cran.r-project.org/web/packages/statgenGWAS/vignettes/GWAS.html
# https://biometris.github.io/statgenGWAS/index.html
# https://www.nature.com/articles/s43586-021-00056-9
# http://ww7.stat-gen.org/
# https://www.mtholyoke.edu/courses/afoulkes/Data/GWAStutorial/

# A guide to genome-wide association analysis and post-analytic interrogation ----

# https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605

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

# Specify parameters to be used in the data processing and analysis
data.dir <- "tutorial/tutorial_files/"
output.dir <- "tutorial/tutorial_output/"

# Input files
gwas.fn <- lapply(c(bed = "bed", bim = "bim", fam = "fam", gds = "gds"),
                  function(n) sprintf("%s/chr16_1000g_CEU.%s", data.dir, n))
clinical.fn <- sprintf("%s/GWAStutorial_clinical.csv", data.dir)
onethou.fn <- lapply(c(info = "info", ped = "ped"),
                     function(n) sprintf("%s/chr16_1000g_CEU/%s", data.dir, n))
protein.coding.coords.fname <- sprintf("%s/ProCodgene_coords.csv", output.dir)

# Output files
gwaa.fname <- sprintf("%s/GWAStutorialout.txt", output.dir)
gwaa.unadj.fname <- sprintf("%s/GWAStutorialoutUnadj.txt", output.dir)
impute.out.fname <- sprintf("%s/GWAStutorial_imputationOut.csv", output.dir)
CETP.fname <- sprintf("%s/CETP_GWASout.csv", output.dir)

# Step 1 - Read and format data

library(snpStats)

geno <- snpStats::read.plink(bed = "tutorial/tutorial_files/GWAStutorial.bed",
                             bim = "tutorial/tutorial_files/GWAStutorial.bim",
                             fam = "tutorial/tutorial_files/GWAStutorial.fam",
                             na.strings = ("-9"))

genotype <- geno$genotypes
print(genotype)
