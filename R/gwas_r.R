# A guide to genome-wide association analysis

# Last updated 2022-02-11

# https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605

# Load packages ----

library(rtracklayer)
library(snpStats)

# Read files ----

# Set paths
bed.file <- "tutorial/tutorial_files/GWAStutorial.bed"
bim.file <- "tutorial/tutorial_files/GWAStutorial.bim"
fam.file <- "tutorial/tutorial_files/GWAStutorial.fam"

# Read in PLINK files to create list
geno <- snpStats::read.plink(bed = bed.file,
                             bim = bim.file,
                             fam = fam.file,
                             na.strings = "-9")
