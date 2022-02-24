# Tutorials

# Last updated 2022-02-24

# Load packages ----

library(tidyverse)

# Population structure ----

# Load data
data <- vroom::vroom(file = "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt",
                     col_names = FALSE)

dim(data)  # 1865 genotypes, 12038 markers

distinct(data[, 1])
