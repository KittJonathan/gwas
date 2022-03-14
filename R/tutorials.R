# Tutorials

# Last updated 2022-02-24

# Load packages ----

library(tidyverse)
library(ape)
library(broom)
library(FactoMineR)
library(ggtree)
library(LEA)


# Population structure ----

# Load data

geno <- vroom::vroom(file = "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt",
                     col_names = FALSE)

dim(geno)  # 1865 genotypes, 12038 markers

distinct(geno[, 1])  # allele coding

# Visualise total population distribution for the first 3 markers

d1 <- geno %>% 
  select(marker1 = X1, marker2 = X2, marker3 = X3) %>%  # select and rename first 3 markers 
  tibble::add_column(genotype = 1:nrow(.), .before = "marker1") %>%  # add row with genotype number 
  pivot_longer(-genotype, names_to = "marker", values_to = "allele") %>%  # long format 
  filter(!is.na(allele)) %>%  # remove missing data 
  count(marker, allele) %>%  # count number of alleles for each marker
  mutate(allele = factor(allele))  # allele as factor

ggplot(data = d1, mapping = aes(x = marker, y = n, fill = allele)) +  # initiate ggplot
  geom_bar(stat = "identity", position = position_dodge(),  # barplot
           width = 0.5) +
  scale_fill_brewer(palette = "Blues") +  # set colours
  ggtitle("Genotypes for first 3 markers in total population") +  # add title
  labs(x = "") +  # remove x-axis label
  theme_minimal() +  # minimal ggplot theme
  theme(panel.grid.major.x = element_blank(),  # remove major x-axis grid
        plot.title = element_text(hjust = 0.5))  # center title

# Tree

geno <- geno %>% 
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))  # replace NAs by MAFs

subset_mk <- geno %>% 
  select(sample(x = 1:ncol(.), size = 2500, replace = FALSE))  # sample 2500 markers

distances <- subset_mk %>% 
  as.matrix() %>% 
  dist()  # calculate distance matrix

arbol <- nj(distances)  # neighbour-joining tree estimation

ggtree(tr = arbol, layout = "ape")  # View tree

tree <- hclust(distances, method = "ward.D2")  # hierarchical clustering

ggtree(tr = tree, layout = "dendrogram")  # View hierarchical tree

plot.new()
plot(tree)
rect.hclust(tree = tree, k = 3)

tree <- HCPC(res = as.matrix(subset_mk))

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
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

d1 <- augment(pca_fit, data = geno)

ggplot(d1, mapping = aes(.fittedPC1, .fittedPC2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal()

# SNMF

snmf_object <- snmf(input.file = "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.geno",
                    K = 1:10, repetitions = 1, entropy = TRUE,
                    ploidy = 2, project = "new", CPU = 2)

snmf_object <- load.snmfProject("TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.snmfProject")

plot(snmf_object)

ce_values <- list()

for (i in 1:10) {
  ce_values[[i]] <- cross.entropy(snmf_object, K = i)
}

cross_entropies <- tibble(
  k = 1:10,
  ce = unlist(ce_values))

cross_entropies %>% 
  ggplot(aes(x = k, y = ce)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Number of ancestral populations",
       y = "Cross-entropy") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank())

head(Q(snmf_object, K = 3))

q_matrix <- Q(snmf_object, K = 3) %>% 
  as_tibble()

dim(q_matrix)

admixture <- q_matrix %>%
  tibble::rowid_to_column() %>% 
  rename(individual = rowid, group_1 = V1, group_2 = V2, group_3 = V3) %>% 
  pivot_longer(-individual, names_to = "group", values_to = "value") %>% 
  group_by(individual) %>% 
  mutate(likely_assignment = group[which.max(value)],
         assignment_prob = max(value)) %>% 
  arrange(likely_assignment, desc(assignment_prob)) %>% 
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))
  

admixture

ggplot(data = admixture) +
  geom_col(aes(x = individual, y = value, fill = group)) +
  scale_fill_manual(values = c("tomato", "lightblue", "gold"))

vroom::vroom_write(q_matrix, "TD_Structure_et_GWAS1/TD2_Structure/TD2_Structure/structure_k3.txt")

my.colors <- c("tomato", "lightblue","gold")

bp=barchart(snmf_object, K = 3, 
            border = NA, space = 0,
            col = my.colors,
            xlab = "Individuals",
            ylab = "Ancestry proportions",
            main = "Ancestry matrix")

axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)
