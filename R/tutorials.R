# Tutorials

# Last updated 2022-02-24

# Load packages ----

library(tidyverse)
library(tidylog)
library(patchwork)

# Population structure ----

# Load data

geno <- vroom::vroom(file = "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt",
                     col_names = FALSE)

dim(geno)  # 1865 genotypes, 12038 markers

distinct(geno[, 1])

# Visualise total population distribution for the first 3 markers

d1 <- geno %>% 
  select(marker1 = X1, marker2 = X2, marker3 = X3) %>% 
  tibble::add_column(genotype = 1:nrow(.), .before = "marker1") %>% 
  pivot_longer(-genotype, names_to = "marker", values_to = "allele") %>% 
  filter(!is.na(allele)) %>% 
  count(marker, allele) %>% 
  mutate(allele = factor(allele))

d1 %>% 
  filter(marker == "marker1") %>% 
  ggplot(aes(x = allele, y = n)) +
  geom_col()

ggplot(data = d1, mapping = aes(x = marker, y = n, fill = allele)) +
  geom_bar(stat = "identity", position = position_dodge(),
           width = 0.5) +
  scale_fill_brewer(palette = "Blues") +
  ggtitle("Genotypes for first 3 markers in total population") +
  labs(x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
  
  pivot_longer(names_to = "")

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
