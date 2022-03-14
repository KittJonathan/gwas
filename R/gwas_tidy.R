# GWAS the Tidy way

# Last updated 2022-03-14

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

# Load packages ----

library(FactoMineR)
library(ape)
library(ggtree)
library(LEA)
library(tidyverse)
library(vroom)

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




#################################  ACP #########################################

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


################### Add info given by snmf on PCA graph ############

# What is doing cbind function ?
data<-cbind(Group=vec_gr,data)

dim(data)
#Do dimensions have changed ?


plot( acp$row$coord[,"Dim 1"],	# Dim 1 is X axe 
      acp$row$coord[,"Dim 2"],	# Dim 2 is Y axe
      main= "Genetic diversity",	# title
      pch=16,				# circle symbole 
      cex=.5,				# half size symbole 
      asp=1,       # orthonormal basis
      xlab="Axe 1",
      ylab="Axe 2",
      col=c("red", "green", "blue")[as.numeric(as.factor(data[,"Group"]))] # colours by group
)	
abline(h=0,v=0,lty=2)			# adding axes


gr<- aggregate(acp$row$coord,	# factors coordinates table 
               FUN="mean",		# compute mean ...
               by=list(data[,"Group"]))	# ...per genetic group

points(x = gr$`Dim 1`,y=gr$`Dim 2`,cex=2,pch=20)  #add points
text(x = gr$`Dim 1`,y=gr$`Dim 2`,cex=1,labels = c("Group1","Group2","Group3"),adj = c(1.5,-2)) # add point labels



### combine complementary informations about rice genotypes to previous study

genotypes_list=read.table("genotype_list.txt",header=T) # genotypes list
genotypes_list=cbind(genotypes_list,Group=vec_gr) # join "group" information (cbind cause same order but be careful!)

info=read.table("info_taxoPop_traits_rice.txt", header=T) # load complementary info

info=merge(info,genotypes_list) #join two tables


### To which taxonomic population are matching genetic groups 1, 2 and 3 ? 
# calculate the mean membership of each taxonomic population to each genetic group
for (gr in sort(unique(info$Group))){
  cat(gr)
  print(table(info[info$Group==gr,"POPULATION"])/length(info[info$Group==gr,"POPULATION"]))
  cat("\n\n")
}
# with dplyr:
info%>%group_by(Group,POPULATION)%>%summarise(count=n())%>%mutate(mean_membership=count/sum(count))

taxPop=data.frame(info%>%group_by(Group,POPULATION)%>%summarise(count=n())%>%
                    mutate(mean_membership=count/sum(count))%>%filter(mean_membership==max(mean_membership)))

### What do you notice concerning pericarp colour ? Comment
info%>%filter(!is.na(pericarp_color))%>%
  group_by(Group,pericarp_color)%>%summarise(count=n())%>%mutate(mean_membership=count/sum(count))

### What do you notice concerning culm angle distribution by genetic group?
# Comment  (scale 1 erected to 9 no erected)
info%>%filter(culm_angle!=-9)%>%
  group_by(Group,culm_angle)%>%summarise(count=n())%>%mutate(mean_membership=count/sum(count))

## visualization with ggplot2
info%>%filter(culm_angle!=-9)%>%group_by(Group)%>%ggplot(aes(Group,culm_angle))+geom_boxplot()



########################## ACP graph with taxonomic info on groups ################
plot( acp$row$coord[,"Dim 1"],	# Dim 1 is X axe 
      acp$row$coord[,"Dim 2"],	# Dim 2 is Y axe
      main= "Genetic diversity",	# title
      pch=16,				# circle symbole 
      cex=.5,				# half size symbole 
      asp=1,       # orthonormal basis
      xlab="Axe 1",
      ylab="Axe 2",
      col=c("red", "green", "blue")[as.numeric(as.factor(data[,"Group"]))] # colours by group
)	
abline(h=0,v=0,lty=2)			# adding axes


gr<- aggregate(acp$row$coord,	# factors coordinates table 
               FUN="mean",		# compute mean ...
               by=list(data[,"Group"]))	# ...per genetic group

points(x = gr$`Dim 1`,y=gr$`Dim 2`,cex=2,pch=20)  #add points
text(x = gr$`Dim 1`,y=gr$`Dim 2`,cex=1,labels = taxPop$POPULATION,adj = c(1.5,-2)) # add point labels




# info_2=read.table("pheno_GWAS_fastlmm_2.txt",header=T)
# info_2=info_2%>%mutate(GENOTYPE=paste('0_', IID, sep=""))%>%
#   select(pericarp_color, culm_angle, GENOTYPE)
# 
# colnames(info_2)=c("pericarp_color_num","culm_angle","GENOTYPE")
# info=merge(info,info_2)
# write.table(info,"info_classification_traits_rice.txt", row.names=F, col.names=T, sep="\t", quote=F)

