# This script estimates admixture proportions in HGDP samples using
# the EM algorithm.
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
library(parallel)
library(data.table)
source("misc.R")
source("mcmc.R")
source("read.data.R")
source("admixture.R")
dyn.load("mcmc.so")
dyn.load("admixture.so")

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

# LOAD ADMIXTURE OUTPUT
# ---------------------
# Load the admixture proportions computed using ADMIXTURE.
Q0           <- read.table("../data/hgdp.admixture.K=7.admix",sep = " ",
                           header = TRUE,stringsAsFactors = FALSE)
rownames(Q0) <- Q0$id
Q0           <- as.matrix(Q0[-1])

# COMPUTE MAXIMUM LIKELIHOOD ADMIXTURE ESTIMATES USING EM
# -------------------------------------------------------
cat("Estimating admixture proportions in HGDP samples.\n")
r <- system.time(out <- admixture.em(geno,K,e = e,Q = Q0,mc.cores = mc.cores))
cat(sprintf("Computation took %0.1f min.\n",r["elapsed"]/60))
rm(r)

# COMPARE RESULT OF RUNNING ADMIXTURE AND EM ALGORITHM
# ----------------------------------------------------
cat(sprintf("Largest difference in admixture proportions is %0.3f.\n",
            max(abs(Q0 - out$Q))))
cat("Distribution of maximum differences across all HGDP samples:\n")
abs(Q0 - out$Q)
