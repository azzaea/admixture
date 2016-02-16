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

# LOAD ADMIXTURE OUTPUT
# ---------------------
# Load the admixture proportions computed using ADMIXTURE.
Q0           <- read.table("../data/hgdp.admixture.K=7.admix",sep = " ",
                           header = TRUE,stringsAsFactors = FALSE)
rownames(Q0) <- Q0$id
Q0           <- as.matrix(Q0[-1])

# COMPUTE MAXIMUM LIKELIHOOD ADMIXTURE ESTIMATES USING TURBOEM
# ------------------------------------------------------------
# For a random initialization of the admixture proportions instead of
# initializing to the estimates from ADMIXTURE, set Q = NULL.
cat("Estimating admixture proportions in HGDP samples using Turbo-EM.\n")
cat("ALGORITHM: 50 ITERATIONS OF QN, FOLLOWED BY DECME\n")
out <- admixture.em(geno,K,e = e,mc.cores = mc.cores,trace = FALSE,
                    init.iter = 50,method = "decme",tol = 0.01)
with(out$turboem.refine,
     cat(sprintf(paste("Turbo-EM made %d M-step updates, completing",
                       "after %d iterations and %0.1f min.\n"),
                 fpeval,itr,runtime[,"elapsed"]/60)))

stop()

# COMPARE RESULT OF RUNNING ADMIXTURE AND EM ALGORITHM
# ----------------------------------------------------
cat(sprintf("Largest difference in admixture proportions is %0.3f.\n",
            max(abs(Q0 - out$Q))))
cat("Distribution of largest differences across all HGDP samples:\n")
print(round(quantile(apply(abs(Q0 - out$Q),1,max),seq(0,1,0.1)),digits = 3))
