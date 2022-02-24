# Tutorials

# Last updated 2022-02-24

# Load packages ----

library(tidyverse)
library(ape)
library(broom)
library(FactoMineR)


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

# PCA

pca_fit <- prcomp(geno)

pca_fit %>% 
  tidy(matrix = "eigenvalues") %>% 
  filter(PC %in% 1:20) %>% 
  ggplot(mapping = aes(x = PC, y = percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  ggtitle("% of variance explained") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

#################################  ACP #########################################

# How works function CA ? What are the outputs ?
acp<-CA(geno,graph = F)
acp

# Eigen values
kable<- (acp$eig)
## head of the table
kable[1:10,]

# eigen values barplot
barplot(acp$eig[1:20,1], ylab="Eigen values")

# % of explained variance barplot
barplot(acp$eig[1:20,2], ylab="% explianed variance")	

## Factorial plan graphs
# Genotypes graph

plot(acp$row$coord[,"Dim 1"],	# Dim 1 is X axe
     acp$row$coord[,"Dim 2"],	# Dim 2 is Y axe
     main= "Genetic diversity",	# title
     pch=16,				# symbol circle
     cex=.5,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Axe 1",
     ylab="Axe 2"
)
abline(h=0,v=0,lty=2)			# adding lines