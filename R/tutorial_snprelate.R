# Tutorials for the R/Bioconductor package SNPRelate

# https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html

# Last updated 2022-03-16

# Install packages ----

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

# Load packages ----

library(gdsfmt)
library(SNPRelate)

# Preparing data ----

# Genomic Data Structure file format
snpgdsSummary(snpgdsExampleFileName())

# Open a GDS file
genofile <- snpgdsOpen(snpgdsExampleFileName())

# Get the attributes of chromosome coding
get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))
