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

d1 <- augment(pca_fit, data = geno)

ggplot(d1, mapping = aes(.fittedPC1, .fittedPC2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal()





#################################  SNMF ########################################

# open LEA package documentation webpage
browseVignettes("LEA")

## Method snmf (use a specific genotyping file previously formated)
obj_snmf<-snmf(input.file = "geno_filtered_maf005_na010_prunedLD090.geno",
               K = 1:10, repetitions = 1, entropy = T,  
               ploidy = 2, project = "new",CPU = 2)

# SNMF function options:
# repetitions: nb of wanted repetitions
# entropy: "cross-entropy" approach will be applied (TRUE) or not applied (FALSE)
# ploidy: genotypes are diploides (2) or haploides (1)
# project: results will be saved and stored in a new directory 


### IF DON'T WORK AND CAUSES ERROR ==> DOWNLOAD PROJECT RESULTS
obj_snmf = load.snmfProject("Structure_snmf_output/geno_filtered_maf005_na010_prunedLD090.snmfProject")

# What are the files that have been created in your work directory ?

# What is allocated in "obj_snmf" object ?
obj_snmf

# Visualization 
plot(obj_snmf)

# Q(): head of the admixture coefficients matrix
head(Q(obj_snmf, K = 3))
# % d'appartenance de chaque génotype aux 3 groupes


# Allocate the admixture coefficients matrix to the object qmatrix of each individual to each genetic group
qmatrix = Q(obj_snmf, K = 3)
dim(qmatrix)

# save qmatrix in a "txt" file that will be used in GWAS
write.table(qmatrix,"structure_K3.txt", row.names=F,col.names = F, quote=F)

my.colors <- c("tomato", "lightblue","gold")

bp=barchart(obj_snmf, K = 3, 
            border = NA, space = 0,
            col = my.colors,
            xlab = "Individuals",
            ylab = "Ancestry proportions",
            main = "Ancestry matrix")

axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)


# What is this loop for ? (must be feasible with dplyr package)
vec_gr<-NULL
for (i in 1:nrow(qmatrix)){
  vec_gr[i]=paste0("Group",which.max(qmatrix[i,]))
}