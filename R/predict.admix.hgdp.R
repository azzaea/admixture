# This script estimates admixture proportions in HGDP samples using
# the EM algorithm.
#
# Here, I initialize the admixture proportions to the ADMIXTURE
# estimates to show that the EM algorithm recovers a very similar
# solution when the model parameters are initialized in the same
# way. (Note that the solution is not exactly the same because the
# models are slightly different.) You can choose instead a random
# initialization by setting input Q = NULL in the call to
# admixture.em.
#
# INSTRUCTIONS FOR OBTAINING GENOTYPES FROM HUMAN GENOME DIVERSITY PANEL:
#
# Follow these instructions to download the HGDP data. (Note that
# these steps are nearly identical to the instructions provided as
# part of the TeraStructure repository on github.)
#
# 1. Download the data from the website of the Human Genome Diversity
# Project, http://www.hagsc.org/hgdp, and extract the files from the
# ZIP archive: unzip hgdp.zip
#
# 2. Get the sample information from Noah Rosenberg's Annals of Human
# Genetics (2006) article, removing the last line of this file:
#
#   wget http://rosenberglab.stanford.edu/data/rosenberg2006ahg/SampleInformation.txt
#   head -n -1 SampleInformation.txt > temp.txt
#   mv -f temp.txt SampleInformation.txt
#
# 3. Convert the genotypes to the transposed PLINK (.tped) format by
# running the Python script (see the "data" subdirectory) in the same
# directory as the HGDP data.
# 
#   python hgdp2tped
#
# 4. Convert data in transposed PLINK format to the "TRAW" format used
# in PLINK. Here, I filter out SNPs that are not on autosomal
# chromosomes, and I (greedily) remove SNPs so that any pair of SNPs
# within each overlapping 500-kb window has a maximum correlation
# (r^2) of 0.2. The list of 159,846 SNPs obtained after these
# filtering steps is provided in the "data" subdirectory.
#
#  /DNAData/bin/plink2 -tfile hgdp --recode A-transpose spacex \
#     --chr 1-22 --extract hgdp-markers-pruned.txt --out hgdp2
#  mv hgdp2.traw hgdp.traw
#
suppressPackageStartupMessages({
library(parallel)
library(data.table)
library(turboEM)
source("misc.R")
source("mcmc.R")
source("read.data.R")
source("admixture.R")
dyn.load("mcmc.so")
dyn.load("admixture.so")})

# SCRIPT PARAMETERS
# -----------------
K        <- 7      # Number of ancestral populations.
e        <- 0.001  # Probability of a genotype error.
seed     <- 1      # Specifies the sequence of pseudorandom numbers.
mc.cores <- 20     # Number of CPUs to use.

# The .traw file containing the genotypes of the HGDP samples.
traw.file <- "hgdp.traw"

# LOAD GENOTYPES
# --------------
# Read the genotype data from the .traw file.
cat("Loading genotype data from .traw file.\n")
geno <- read.traw.file(traw.file)$geno

# Initialize the random number generator.
set.seed(seed)

# COMPUTE MAXIMUM LIKELIHOOD ADMIXTURE ESTIMATES USING TURBOEM
# ------------------------------------------------------------
cat("Fitting admixture model to data.\n")
r <- system.time(out <-
       admixture.em(geno,K,e = e,tol = 0.01,mc.cores = mc.cores))
with(out,
     cat(sprintf("Turbo-EM completed after %d iterations and %0.1f min.\n",
                 length(loglikelihood),r["elapsed"]/60)))

# SAVE RESULTS TO FILE
# --------------------
save(list = c("K","e","seed","out","r"),
     file = "hgdp.out.RData")
