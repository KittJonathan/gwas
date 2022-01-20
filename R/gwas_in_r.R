# Genome-wide association studies in R

# https://poissonisfish.com/2017/10/09/genome-wide-association-studies-in-r/


# Workflow ----


# Ten steps:
# 1) reading data into R to create an R object
# 2) SNP-level filtering (part 1)
# 3) sample-level filtering
# 4) SNP-level filtering (part 2)
# 5) principal component analysis (PCA)
# 6) imputation of non-typed genotypes
# 7) association analysis of typed SNPs
# 8) association analysis of imputed data
# 9) integration of imputed and typed SNP results
# 10) visualisation and quality control of association findings

# Data pre-processing
# .ped file: information on each study participant
# .map file: 1 row for each SNP with corresponding chromosome and coordinate
# .bim file: same information as the .map file + the two observed alleles at each SNP from the .ped file
# .bed file: binary version of the genotype data
# .fam file: participant identification information (same as .ped, without the genotype)
# clinical data file: .txt or .csv, clinical data on each study subject (covariates and phenotypes)

# Install packages ----

# Packages :
# snpStats: read in various formats of genotype data, carry out QC, imputation & association analysis
# SNPRelate: sample-level QC & computationally efficient principal component calculation
# ggplot2, LDheatmap, postgwas: data visualisation
# plyr: data manipulation
# doParallel: parallel processing

# Install required packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("snpStats")
BiocManager::install("SNPRelate")
BiocManager::install("rtracklayer")
BiocManager::install("biomaRt")

install.packages("tidyverse")
install.packages("GenABEL")
install.packages("LDheatmap")
install.packages("doParallel")
install.packages("coin")
install.packages("igraph")
install.packages("https://github.com/merns/postgwas/releases/download/1.11-2/postgwas_1.11-2.zip", repos=NULL)

library(devtools)
install.packages("http:::cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")

# Load packages ----

library(tidyverse)
library(snpStats)
library(SNPRelate)

# Set paths ----

india_files <- paste("data/public/Genomics/105Indian_2527458snps",
                     c(".bed", ".bim", ".fam"), sep = "")

malaysia_files <- paste("data/public/Genomics/108Malay_2527458snps",
                        c(".bed", ".bim", ".fam"), sep = "")

china_files <- paste("data/public/Genomics/110Chinese_2527458snps",
                        c(".bed", ".bim", ".fam"), sep = "")

# Read files ----

load("data/public/conversionTable.RData")

snps_india <- snpStats::read.plink(bed = india_files[1],
                                   bim = india_files[2],
                                   fam = india_files[3])

snps_malaysia <- snpStats::read.plink(bed = malaysia_files[1],
                                      bim = malaysia_files[2],
                                      fam = malaysia_files[3])

snps_china <- snpStats::read.plink(bed = china_files[1],
                                   bim = china_files[2],
                                   fam = china_files[3])

# Read files ----

china_snps <- snpStats::read.plink(bed = china_files[1],
                                   bim = china_files[2],
                                   fam = china_files[3])

india_snps <- snpStats::read.plink(bed = india_files[1],
                                   bim = india_files[2],
                                   fam = india_files[3])

malaysia_snps <- snpStats::read.plink(bed = malaysia_files[1],
                                      bim = malaysia_files[2],
                                      fam = malaysia_files[3])

rm(china_files, india_files, malaysia_files)
all.equal(ncol(snps_china$genotypes), ncol(snps_india$genotypes), ncol(snps_malaysia$genotypes))

load("data/conversionTable.RData")

# Merge the three SNP datasets ----

snps <- snps_malaysia
snps$genotypes <- rbind(snps_malaysia$genotypes, snps_india$genotypes, snps_china$genotypes)
colnames(snps$map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
snps$fam <- rbind(snps_malaysia$fam, snps_india$fam, snps_china$fam)

# Rename SNPs present in the conversion table into rs IDs ----

mapped_snps <- intersect(snps$map$SNP, names(conversionTable))
new_ids <- conversionTable[match(snps$map$SNP[snps$map$SNP %in% mapped_snps],
                                 names(conversionTable))]
snps$map$SNP[rownames(snps$map) %in% mapped_snps] <- new_ids

# Import lipid datasets & match SNP-Lipidomics samples ----

malaysia_lipids <- read.delim("data/public/Lipidomic/117Malay_282lipids.txt", row.names = 1)
india_lipids <- read.delim("data/public/Lipidomic/120Indian_282lipids.txt", row.names = 1)
china_lipids <- read.delim("data/public/Lipidomic/122Chinese_282lipids.txt", row.names = 1)

all(Reduce(intersect, list(colnames(malaysia_lipids),
                           colnames(india_lipids),
                           colnames(china_lipids))) == colnames(malaysia_lipids))

lipids <- rbind(malaysia_lipids, india_lipids, china_lipids)

# Country ----

country <- sapply(list(snps_malaysia, snps_india, snps_china), function(k){
  nrow(k$genotypes)
})

origin <- data.frame(sample.id = rownames(snps$genotypes),
                     country = factor(rep(c("M", "I", "C"), country)))

matchingSamples <- intersect(rownames(lipids), rownames(snps$genotypes))

snps$genotypes <- snps$genotypes[matchingSamples, ]
lipids <- lipids[matchingSamples, ]
origin <- origin[match(matchingSamples, origin$sample.id), ]

# Combine SNPs and lipidomics ----

genData <- list(SNP = snps$genotypes,
                MAP = snps$map,
                LIP = lipids)

# Write processed omics and GDS ----

save(genData, origin, file = "output/PhenoGenoMap.RData")
snpStats::write.plink("output/convertGDS", snps = snps$genotypes)

# Clear memory
rm(list = ls())

# Pre-processing ----

# Load files
load("output/PhenoGenoMap.RData")

# Use SNP call rate of 100% & MAF of 0.1 (very stringent)
maf <- 0.1
callRate <- 1

SNPstats <- snpStats::col.summary(genData$SNP)

maf_call <- with(SNPstats, MAF > maf & Call.rate == callRate)
genData$SNP <- genData$SNP[, maf_call]
genData$MAP <- genData$MAP[maf_call, ]
SNPstats <- SNPstats[maf_call, ]

# Sample call rate & heterozygosity
callMat <- !is.na(genData$SNP)
Sampstats <- snpStats::row.summary(genData$SNP)

hetExp <- callMat %*% (2 * SNPstats$MAF * (1 - SNPstats$MAF))
hetObs <- with(Sampstats, Heterozygosity * (ncol(genData$SNP)) * Call.rate)
Sampstats$hetF <- 1 - (hetObs / hetExp)

# Use sample call rate of 100%, het threshold of 0.1 (veryf stringent)
het <- 0.1
het_call <- with(Sampstats, abs(hetF) < het & Call.rate == 1)
genData$SNP <- genData$SNP[het_call, ]
genData$LIP <- genData$LIP[het_call, ]

# LD and kinship coeff
ld <- 0.2
kin <- 0.1

SNPRelate::snpgdsBED2GDS(bed.fn = "output/convertGDS.bed",
                         bim.fn = "output/convertGDS.bim",
                         fam.fn = "output/convertGDS.fam",
                         out.gdsfn = "output/myGDS",
                         cvt.chr = "char")

genofile <- SNPRelate::snpgdsOpen("output/myGDS", readonly = FALSE)
gds.ids <- gdsfmt::read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
gdsfmt::add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)
geno.sample.ids <- rownames(genData$SNP)

# First filter for LD
snpSUB <- SNPRelate::snpgdsLDpruning(genofile, ld.threshold = ld,
                                     sample.id = geno.sample.ids,
                                     snp.id = colnames(genData$SNP))
snpset.ibd <- unlist(snpSUB, use.names = FALSE)

# And now filter for MoM
ibd <- SNPRelate::snpgdsIBDMoM(genofile, kinship = TRUE,
                               sample.id = geno.sample.ids,
                               snp.id = snpset.ibd,
                               num.thread = 1)
ibdcoef <- SNPRelate::snpgdsIBDSelection(ibd)
ibdcoef <- ibdcoef[ibdcoef$kinship >= kin, ]

# Filter samples out
related.samples <- NULL
while (nrow(ibdcoef) > 0) {
  # count the number of occurrences of each and take the top one
  sample.counts <- sort(table(c(ibdcoef$ID1, ibdcoef$ID2)), decreasing = T)
  rm.sample <- names(sample.counts)[1]
  cat("Removing sample", rm.sample, "too closely related to",
      sample.counts[1], "other samples.\n")
  
  # remove from ibdcoef and add to list
  ibdcoef <- ibdcoef[ibdcoef$ID1 != rm.sample & ibdcoef$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

genData$SNP <- genData$SNP[!(rownames(genData$SNP) %in% related.samples), ]
genData$LIP <- genData$LIP[!(rownames(genData$LIP) %in% related.samples), ]

# Principal Component Analysis ----

# PCA
set.seed(100)
pca <- SNPRelate::snpgdsPCA(genofile, sample.id = geno.sample.ids,
                            snp.id = snpset.ibd, num.thread = 1)
pcatab <- data.frame(sample.id = pca$sample.id,
                     PC1 = pca$eigenvect[, 1],
                     PC2 = pca$eigenvect[, 2],
                     stringsAsFactors = FALSE)

# Subset and/or reorder origin accordingly
origin <- origin[match(pca$sample.id, origin$sample.id), ]

pcaCol <- rep(rgb(0, 0, 0, 0.3), length(pca$sample.id))
pcaCol[origin$country == "I"] <- rgb(1, 0, 0, 0.3)
pcaCol[origin$country == "M"] <- rgb(0, 0.7, 0, 0.3)

png("output/PCApopulation.png", width = 500, height = 500)
plot(pcatab$PC1, pcatab$PC2, xlab = "PC1", ylab = "PC2", col = pcaCol, pch = 16)
abline(h = 0, v = 0, lty = 2, col = "grey")
legend("top", legend = c("Chinese", "Indian", "Malay"), col = 1:3, pch = 16, bty = "n")
dev.off()

# Genome-Wide Association ----

source("R/GWASfunction.R")

target <- "Cholesterol"

phenodata <- data.frame("id" = rownames(genData$LIP),
                        "phenotype" = scale(genData$LIP[, target]), stringsAsFactors = FALSE)

start <- Sys.time()
GWAA(genodata = genData$SNP,
     phenodata = phenodata,
     filename = paste(target, ".txt", sep = ""))
Sys.time() - start  # 51 mins

# Manhattan plot
GWASout <- read.table("Cholesterol.txt", header = TRUE,
                      colClasses = c("character", rep("numeric", 4)))
GWASout$type <- rep("typed", nrow(GWASout))
GWASout$Neg_logP <- -log10(GWASout$p.value)

GWASout <- merge(GWASout, genData$MAP[, c("SNP", "chr", "position")])
GWASout <- GWASout[order(GWASout$Neg_logP, decreasing = TRUE), ]

readr::write_tsv(GWASout, "output/gwas_out.txt")

png(paste(target, ".png", sep = ""), height = 500, width = 500)
GWAS_Manhattan(GWASout)
dev.off
