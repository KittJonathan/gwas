# Tutorials

# Last updated 2022-02-24

# Load packages ----

library(tidyverse)
library(patchwork)

# Population structure ----

# Load data

geno <- vroom::vroom(file = "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt",
                     col_names = FALSE)

dim(geno)  # 1865 genotypes, 12038 markers

distinct(geno[, 1])

# Visualise data

h1 <- ggplot(data = geno, mapping = aes(x = X1)) +
  geom_histogram(fill = "lightblue")

h1

hist(data$X1, col = "lightblue", main = "Genotypes for 1st marker in total population", xlab = "Genotypes")

ggplot(data = data, mapping = aes(x = X3)) +
  geom_histogram()

## Visualize total population distribution for 3 firt markers
par(mfrow=c(2,2)) # put graphs on 1 row by 2 columns
hist(data$X1,col="lightblue", main = "Genotypes pour le mk 1 dans la population totale",
     xlab="Genotypes")
hist(data$X2,col="grey", main = "Genotypes pour le mk 2 dans la population totale",
     xlab="Genotypes")
hist(data$X3,col="plum", main = "Genotypes pour le mk 3 dans la population totale",
     xlab="Genotypes")
