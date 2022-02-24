# Tutorials

# Last updated 2022-02-24

# Load packages ----

# Population structure ----

# Load data
data <- read_delim(file = "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt")

#Define work directory
setwd("C:/Users/plasserre/Documents/TP_LASSERRE-Z_2021/TD2_Structure/")

data=read.table("geno_filtered_maf005_na010_prunedLD090.txt",header=F)

# How many genotypes, how many markers ?
dim(data)

# How alleles are coded ?
table(data[,1], useNA = "always")