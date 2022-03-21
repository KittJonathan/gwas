# GWAS the Tidy way

# Last updated 2022-03-21

# Links ----

# https://luisdva.github.io/rstats/model-cluster-plots/

# Install CRAN packages ----

# install.packages("ape")
# install.packages("ggtree)
# install.packages("FactoMineR")
# install.packages("tidyverse")

# Install Bioconductor packages ----

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("LEA")

BiocManager::install("trio")

# Load packages ----

library(FactoMineR)
library(ape)
library(ggtree)
library(LEA)
library(tidyverse)
library(vroom)
library(broom)
library(trio)

# Import .ped and .fam files ----

test <- read.pedfile(file = "data/transfer_2986650_files_a7f88f43/geno.ped")
test <- snpStats::read.pedfile("data/transfer_2986650_files_a7f88f43/geno.ped")

map_file <- vroom("data/transfer_2986650_files_a7f88f43/geno.map", col_names = FALSE)
# X1 = chromosome
# X2 = SNP
# X3 = Genetic Distance
# X4 = Position

ped_file <- vroom::vroom("data/transfer_2986650_files_a7f88f43/geno.ped")

# Import genotyping data ----

genotyping_data <- vroom(file = "data/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt",
                         col_names = FALSE)

# Total population distribution for the first three markers ----

genotyping_data %>% 
  select(1:3, ) %>%   # keep data for the first three markers
  pivot_longer(everything(), names_to = "marker", values_to = "allele") %>%   # long format
  mutate(marker = str_remove(string = marker, pattern = "X")) %>%  # remove "X" in marker name
  filter(!is.na(allele)) %>%  # remove NAs
  mutate(marker = fct_inseq(f = factor(marker)),  # transform variables into factors
         allele = factor(allele, levels = c(0, 1))) %>% 
  count(marker, allele) %>%  # count number of individuals for each marker and allele
  ggplot() +  # initiate ggplot
  geom_bar(aes(x = paste0("Marker ", marker), y = n, fill = allele),  # set x, y and fill for geom_bar()
           stat = "identity", position = position_dodge(),  # set stat and position for geom_bar()
           width = 0.5) +  # reduce bar width  
  scale_fill_brewer(palette = "Blues") +  # change colours
  labs(title = "Allele distribution for the first three markers",  # add title to plot
       x = "", y = "") +  # remove x and y axis labels
  theme_minimal() +  # use minimal ggplot theme
  theme(panel.grid.major.x = element_blank(),  # remove major grids on x axis
        panel.grid.minor.y = element_blank())  # remove minor grids on y axis

# Replace NAs by MAF ----

genotyping_data <- genotyping_data %>% 
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))

# Compute genetic distances for a subset of 2,500 markers ----

distances <- genotyping_data %>% 
  select(sample(1:ncol(.), size = 2500, replace = FALSE)) %>%  # sample 2500 markers
  as.matrix() %>%  # transform into matrix
  dist()  # calculate distances

# Compute neighbor-joining tree estimation ----

nj_tree_estimation <- nj(distances)

# View tree ----

ggtree(tr = nj_tree_estimation,
      layout = "ape")

# Hierarchical clustering ----

hc_tree <- hclust(d = distances, method = "ward.D2")

# View hierarchical clusters ----

ggtree(tr = hc_tree, 
       layout = "dendrogram")

# Conduct Principal Component Analysis ----

pca_fit <- prcomp(genotyping_data)

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

d1 <- augment(pca_fit, data = genotyping_data)

ggplot(d1, mapping = aes(.fittedPC1, .fittedPC2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Genetic diversity", x = "PC1", y = "PC2") +
  theme_minimal()

# Estimate individual ancestry coefficients and ancestral allele frequencies ----

LEA::snmf(input.file = "data/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.geno",
          K = 1:10,
          repetitions = 1,
          entropy = TRUE,
          ploidy = 2,
          project = "new",
          CPU = 2)

# Import the Sparse Non-negative Matrix Factorization algorithm results ----

snmf <- load.snmfProject("data/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.snmfProject")

# Extract cross-entropy values ----

ce_list <- list()

for (i in 1:length(snmf@runs)) {
  
  ce_list[[i]] <- LEA::cross.entropy(object = snmf, K = i)
  
}

ce_values <- tibble(
  k = 1:length(snmf@runs),
  ce_value = unlist(ce_list))

ce_values %>% 
  ggplot(mapping = aes(x = k, y = ce_value)) +
  geom_point(colour = "#56B4E9") +
  geom_line(colour = "#56B4E9") +
  scale_x_continuous(breaks = 1:length(snmf@runs)) +
  ggtitle("Cross-entropy") +
  labs(x = "Number of ancestral populations",
       y = "Cross-entropy") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Extract and save the admixture coefficients table ----

q_matrix <- LEA::Q(object = snmf, K = 3) %>% 
  as_tibble()

vroom::vroom_write(q_matrix, "data/TD2_Structure/q_matrix_k3.txt")

# Assign groups to individuals ----

group_assignment <- q_matrix %>%
  tibble::rowid_to_column() %>% 
  pivot_longer(-rowid, names_to = "group", values_to = "prob") %>% 
  rename(individual = rowid) %>% 
  mutate(group = str_remove(group, "V")) %>% 
  mutate(individual = fct_inseq(factor(individual)),
         group = fct_inseq(factor(group))) %>% 
  group_by(individual) %>% 
  mutate(assigned_group = group[which.max(prob)],
         assigned_prob = max(prob)) %>% 
  arrange(assigned_group, desc(assigned_prob))

# Admixture plot ----

ggplot(data = group_assignment) +
  geom_col(aes(x = individual, y = prob, fill = group)) +
  scale_fill_manual(values = c("lightblue", "tomato", "gold"))
  # scale_fill_manual(values = c("tomato", "lightblue", "gold")) 
  

my.colors <- c("tomato", "lightblue","gold")

bp=barchart(snmf, K = 3, 
            border = NA, space = 0,
            col = my.colors,
            xlab = "Individuals",
            ylab = "Ancestry proportions",
            main = "Ancestry matrix")

axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

# Add snmf groups to PCA data ----

groups_summary <- group_assignment %>% 
  group_by(individual) %>% 
  filter(row_number() == 1) %>% 
  select(individual, assigned_group)

pca_groups <- d1 %>% 
  rename(individual = .rownames) %>% 
  left_join(groups_summary)

# Plot PCA results with SNMF groups ----

ggplot(pca_groups, mapping = aes(.fittedPC1, .fittedPC2)) +
  geom_point(aes(colour = assigned_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Genetic diversity", x = "PC1", y = "PC2") +
  theme_minimal()

# Add complementary informations about rice genotypes to study ----

info <- read_tsv("data/TD2_Structure/info_taxoPop_traits_rice.txt")

genotypes_list <- read_csv("data/TD2_Structure/genotype_list.txt") %>% 
  left_join(info) %>% 
  rowid_to_column() %>% 
  rename(individual = rowid) %>% 
  mutate(individual = factor(individual)) %>% 
  left_join(groups_summary)

# Calculate mean membership of each taxonomic population to each genetic group ----

tax_pop <- genotypes_list %>% 
  group_by(assigned_group, POPULATION) %>% 
  summarise(count = n()) %>% 
  mutate(mean_membership = count / sum(count)) %>% 
  filter(mean_membership == max(mean_membership))

for (gr in sort(unique(genotypes_list$assigned_group))){
  cat(gr)
  print(table(info[genotypes_list$assigned_group == gr,"POPULATION"])/length(genotypes_list[genotypes_list$assigned_group==gr,"POPULATION"]))
  cat("\n\n")
}

# Pericarp colour ----

genotypes_list %>% 
  filter(!is.na(pericarp_color)) %>% 
  group_by(assigned_group, pericarp_color) %>% 
  summarise(count = n()) %>% 
  mutate(mean_membership = count / sum(count))

# Culm angle distribution by genetic group ----

genotypes_list %>% 
  filter(culm_angle != -9) %>% 
  group_by(assigned_group, culm_angle) %>% 
  summarise(count = n()) %>% 
  mutate(mean_membership = count / sum(count))

# Generate boxplot ----

genotypes_list %>% 
  filter(culm_angle != -9) %>% 
  group_by(assigned_group) %>% 
  ggplot(aes(assigned_group, culm_angle)) +
  geom_boxplot()

# PCA plot with taxonomy information ----

pca_taxo <- d1 %>% 
  mutate(individual = factor(.rownames)) %>% 
  select(individual, .fittedPC1, .fittedPC2) %>% 
  left_join(genotypes_list) %>% 
  group_by(assigned_group) %>% 
  mutate(mean_x = mean(.fittedPC1),
         mean_y = mean(.fittedPC2)) %>% 
  ungroup()
  
pca_taxo

ggplot() +
  geom_point(data = pca_taxo, aes(x = .fittedPC1, y = .fittedPC2,
                                  colour = assigned_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Genetic diversity", x = "PC1", y = "PC2") +
  geom_text(data = pca_taxo %>% group_by(assigned_group) %>% filter(row_number() == 1),
            aes(x = mean_x, y = mean_y, label = POPULATION)) +
  theme_minimal()


# Load data for GWAS ----

genotyping_data <- vroom::vroom("data/TD3_GWAS1/GWAS1_2000SNP_01_allele_code.txt") %>% 
  mutate(IID = paste(FID, IID, sep = "_")) %>% 
  select(-FID, -(PAT:PHENOTYPE)) %>% 
  rename(GENOTYPE = IID)

info <- read_tsv("data/TD3_GWAS1/info_taxoPop_traits_rice.txt") %>% 
  select(-SUBPOPULATION) %>% 
  rename(pericarp_color_char = pericarp_color) %>% 
  rename(pericarp_color = pericarp_color_num)

structure <- read_tsv("data/TD3_GWAS1/matrice_covar_structure_K3.txt")

# Merge data ----

d1 <- info %>% 
  left_join(structure) %>% 
  left_join(genotyping_data)

rm(genotyping_data, info, structure)

# Set first genotype data column ----

first_SNP_col <- 9

# Naive association model ----

naive_res_list <- list()

for (i in first_SNP_col:ncol(d1)) {
  
  snp <- colnames(d1[i])
  
  d1_tmp <- d1 %>% 
    filter(culm_angle != -9 & pericarp_color != "-9") %>% 
    select(1:5, i) %>% 
    rename(SNP = all_of(snp))
  
  pericarp_model <- lm(formula = pericarp_color ~ SNP, data = d1_tmp)
  
  culm_model <- lm(formula = culm_angle ~ SNP, data = d1_tmp)
  
  pericarp_res <- pericarp_model %>% 
    tidy() %>% 
    filter(term == "SNP") %>% 
    mutate(term = snp) %>% 
    select(term, pericarp_pval = p.value)
  
  culm_res <- culm_model %>% 
    tidy() %>% 
    filter(term == "SNP") %>% 
    mutate(term = snp) %>% 
    select(term, culm_pval = p.value)
  
  naive_res_list[[i]] <- pericarp_res %>% 
    left_join(culm_res)
  
}

naive_res_table <- do.call("rbind", naive_res_list)


# Load SNP positions ----

pos <- read_tsv("data/TD3_GWAS1/GWAS1_2000SNP_physical_positions.txt") %>% 
  mutate(SNP = as.character(SNP))

# Format naive_res_table and add SNP positions ----

naive_res_table <- naive_res_table %>% 
  mutate(term = str_sub(string = term, start = 1, end = 9)) %>% 
  rename(SNP = term) %>% 
  mutate(pericarp_log10_pval = -log10(pericarp_pval),
         culm_log10_pval = -log10(culm_pval)) %>% 
  left_join(pos) %>% 
  select(SNP, Chromosome:Position, everything())

# Bonferroni correction ----

signif_threshold = -log10(0.001 / ncol(d1))

# Manhattan plot for pericarp color with naive model ----

ggplot(data = naive_res_table,
       aes(x = Position, y = pericarp_log10_pval)) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = seq(0, 30e6, 2e6), labels = seq(0, 30, 2)) +
  geom_hline(yintercept = signif_threshold, colour = "blue", linetype = "dashed") +
  labs(x = "chromosome 7 positions (Mb)",
     y = "pericarp log10 pval") +
  ggtitle("Manhattan plot for pericarp color with naive model")

# Manhattan plot for culm angle with naive model ----

ggplot(data = naive_res_table,
       aes(x = Position, y = culm_log10_pval)) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = seq(0, 30e6, 2e6), labels = seq(0, 30, 2)) +
  geom_hline(yintercept = signif_threshold, colour = "blue", linetype = "dashed") +
  labs(x = "chromosome 7 positions (Mb)",
       y = "culm angle log10 pval") +
  ggtitle("Manhattan plot for culm angle color with naive model")

# Structure corrected association model ----

structure_res_list <- list()

for (i in first_SNP_col:ncol(d1)) {
  
  snp <- colnames(d1[i])
  
  d1_tmp <- d1 %>% 
    filter(culm_angle != -9 & pericarp_color != "-9") %>% 
    select(1:8, i) %>% 
    rename(SNP = all_of(snp))
  
  pericarp_model <- lm(formula = pericarp_color ~ SNP + Group1 + Group2 + Group3, data = d1_tmp)
  
  culm_model <- lm(formula = culm_angle ~ SNP + Group1 + Group2 + Group3, data = d1_tmp)
  
  pericarp_res <- pericarp_model %>% 
    tidy() %>% 
    filter(term == "SNP") %>% 
    mutate(term = snp) %>% 
    select(term, pericarp_pval = p.value)
  
  culm_res <- culm_model %>% 
    tidy() %>% 
    filter(term == "SNP") %>% 
    mutate(term = snp) %>% 
    select(term, culm_pval = p.value)
  
  structure_res_list[[i]] <- pericarp_res %>% 
    left_join(culm_res)
  
}

structure_res_table <- do.call("rbind", structure_res_list)

# Format structure_res_table and add SNP positions ----

structure_res_table <- structure_res_table %>% 
  mutate(term = str_sub(string = term, start = 1, end = 9)) %>% 
  rename(SNP = term) %>% 
  mutate(pericarp_log10_pval = -log10(pericarp_pval),
         culm_log10_pval = -log10(culm_pval)) %>% 
  left_join(pos) %>% 
  select(SNP, Chromosome:Position, everything())

# Bonferroni correction ----

signif_threshold = -log10(0.001 / ncol(d1))

# Manhattan plot for pericarp color with structure corrected model ----

ggplot(data = structure_res_table,
       aes(x = Position, y = pericarp_log10_pval)) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = seq(0, 30e6, 2e6), labels = seq(0, 30, 2)) +
  geom_hline(yintercept = signif_threshold, colour = "blue", linetype = "dashed") +
  labs(x = "chromosome 7 positions (Mb)",
       y = "pericarp log10 pval") +
  ggtitle("Manhattan plot for pericarp color with structure corrected model")

# Manhattan plot for culm angle with structure corrected model ----

ggplot(data = structure_res_table,
       aes(x = Position, y = culm_log10_pval)) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = seq(0, 30e6, 2e6), labels = seq(0, 30, 2)) +
  geom_hline(yintercept = signif_threshold, colour = "blue", linetype = "dashed") +
  labs(x = "chromosome 7 positions (Mb)",
       y = "culm angle log10 pval") +
  ggtitle("Manhattan plot for culm angle color with structure corrected model")

# Select most associated SNP ----

asso_snp <- structure_res_table %>% 
  filter(pericarp_log10_pval == max(pericarp_log10_pval)) %>% 
  pull(SNP)

# Create contingency table ----

allelic_effect <- d1 %>% 
  select(pericarp_color_char, starts_with(asso_snp)) %>% 
  filter(!is.na(pericarp_color_char)) %>% 
  rename(snp_218399762_A=`218399762_A`)
  
# Generate plot ----

library(vcd)

mosaic(snp_218399762_A ~ pericarp_color_char, data = allelic_effect)
