# GWAS the Tidy way

# Last updated 2022-03-14

# Install CRAN packages ----

# install.packages("ape")
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
library(LEA)
library(tidyverse)
library(vroom)

# Import genotyping data ----

genotyping_data <- vroom(file = "data/TD2_Structure/geno_filtered_maf005_na010_prunedLD090.txt",
                         col_names = FALSE)

# Clean genotyping data ----

genotyping_data_long <- genotyping_data %>% 
  rowid_to_column() %>%  # add rowid
  pivot_longer(-rowid,  # transform table into long format
               names_to = "marker",
               values_to = "genotype") %>% 
  rename(individual = rowid) %>%  # rename first column
  mutate(marker = str_remove(string = marker,  # remove "X" from marker column
                             pattern = "X")) %>% 
  mutate(individual = fct_inseq(f = factor(individual)),  # transform variables into factors
         marker = fct_inseq(f = factor(marker)),
         genotype = factor(x = genotype,
                           levels = c(0, 1, NA)))

head(genotyping_data_long)
  

d1 <- geno %>% 
  select(marker1 = X1, marker2 = X2, marker3 = X3) %>%  # select and rename first 3 markers 
  tibble::add_column(genotype = 1:nrow(.), .before = "marker1") %>%  # add row with genotype number 
  pivot_longer(-genotype, names_to = "marker", values_to = "allele") %>%  # long format 
  filter(!is.na(allele)) %>%  # remove missing data 
  count(marker, allele) %>%  # count number of alleles for each marker
  mutate(allele = factor(allele))  # allele as factor

################################################

# How alleles are coded ?
table(data[,1], useNA = "always")


################################ Data Visualization ############################
## Visualize total population distribution for 3 firt markers
par(mfrow=c(2,2)) # put graphs on 1 row by 2 columns
hist(as.numeric(data[,1]),col="lightblue", main = "Genotypes pour le mk 1 dans la population totale",
     xlab="Genotypes")
hist(as.numeric(data[,2]),col="grey", main = "Genotypes pour le mk 2 dans la population totale",
     xlab="Genotypes")
hist(as.numeric(data[,3]),col="plum", main = "Genotypes pour le mk 3 dans la population totale",
     xlab="Genotypes")


################################  TREE #########################################

# Replace NA by minor allele frequency computed as follow:
for (c in 1:ncol(data)){
  data[is.na(data[,c]),c]=mean(x = (data[,c]),na.rm=T)
}

# How works functions nj and dist ?
# How are computed genetic distances

subset_mk=sample(x = 1:ncol(data),size = 2500,replace = F) # sample of 2500 SNP to lighten the data
distances=dist(as.matrix(data[,subset_mk]))
arbol <- nj(distances)

tree=hclust(distances,method = "ward.D2")
plot(tree)
rect.hclust(tree,k=3)
tree <- HCPC(res = as.matrix(data[, subset_mk]))

par(mfrow=c(1,1)) # put graphs on 1 row by 2 columns
# visualization of the tree
plot(arbol,"unrooted", cex=0.5)


#################################  ACP #########################################

# How works function CA ? What are the outputs ?
acp<-CA(data,graph = F)
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
