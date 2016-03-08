# TO DO: Explain what this script does, and give instructions on how
# to use it. Also, add note here about ADMIXTURE with .bed files.
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
