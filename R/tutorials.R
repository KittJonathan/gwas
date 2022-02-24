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

ggplot(data = d1, mapping = aes(x = marker, y = n, fill = allele)) +
  geom_bar(stat = "identity", position = position_dodge(),
           width = 0.5) +
  scale_fill_brewer(palette = "Blues") +
  ggtitle("Genotypes for first 3 markers in total population") +
  labs(x = "") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5))