# This script demonstrates application of the EM algorithm together
# with a previously computed set of allele frequencies (the "panel")
# to project a set of samples onto the K ancestral populations; that
# is, we compute admixture estimates conditioned on the allele
# frequencies for a set of K populations. This reproduces the -P
# option that was introduced in version 1.3.0 of the ADMIXTURE
# software.
#
# To run this script, you can find the allele frequencies file
# (1kg_hgdp.7.P) from the github repository we created for a recent
# workshop presented at Stanford:
#
#   http://github.com/Ancestry/cehg16-workshop
#
# To obtain the genotypes for the 100 test samples (file
# 1kg_hgdp_test.traw), download the files from the same github
# repository, including file 1kg_hgdp_test.bed. Once you have the .bed
# file, convert this to the TRAW format using the following PLINK
# command:
#
#   plink2 -bfile 1kg_hgdp_test --recode A-transpose spacex
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
e        <- 0.001  # Probability of a genotype error.
seed     <- 1      # Specifies the sequence of pseudorandom numbers.
mc.cores <- 20     # Number of CPUs to use.

# These variables specify the input files: the output from ADMIXTURE
# containing allele frequency estimates for each ancestral population,
# and the .traw file containing the genotype data.
freq.file <- "1kg_hgdp.7.P"
traw.file <- "1kg_hgdp_test.traw"

# LOAD GENOTYPES
# --------------
# Read the genotype data from the .traw file.
cat("Loading genotype data from .traw file.\n")
geno <- read.traw.file(traw.file)$geno

# LOAD ALLELE FREQUENCIES
# -----------------------
# Load the allele frequencies from the specified .P file. Since the
# encoding in the .traw file is the opposite of what is used by
# ADMIXTURE, I need to take the inverse frequencies.
cat("Loading allele frequencies from .P file.\n")
F           <- 1 - as.matrix(read.table(freq.file))
K           <- ncol(F)
colnames(F) <- paste0("K",1:K)

# Initialize the random number generator.
set.seed(seed)

# COMPUTE MAXIMUM LIKELIHOOD ADMIXTURE ESTIMATES USING TURBOEM
# ------------------------------------------------------------
cat("Fitting admixture model to data.\n")
r <- system.time(out <- admixture.em(geno,K,e = e,F = F,update.F = FALSE,
                                     init.iter = 10,mc.cores = mc.cores))
with(out,
     cat(sprintf("Turbo-EM completed after %d iterations and %0.1f min.\n",
                 length(loglikelihood),r["elapsed"]/60)))
