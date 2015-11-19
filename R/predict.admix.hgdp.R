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
out  <- read.traw.file(sprintf("%s/%s",data.dir,traw.file))
map  <- out$map
geno <- out$geno
rm(out)

# Initialize the random number generator.
set.seed(seed)

# COMPUTE L0-PENALIZED ADMIXTURE ESTIMATES USING EM
# -------------------------------------------------
cat("Estimating admixture proportions in HGDP samples.\n")
cat("EM algorithm setings:\n")
cat("  Ancestral populations     ",K,"\n")
cat("  Strength of L0-penalty    ",a,"\n")
cat("  Genotype error probability",e,"\n")
n <- nrow(panel)
r <- system.time(out <-
       admixture.em(geno,K,e = e,a = a,F = out.admix$F,Q = out.admix$Q,
                    exact.q = FALSE,T = T,tolerance = 1e-4,cg = TRUE,
                    mc.cores = mc.cores))
rownames(out$Q) <- panel$id
colnames(out$Q) <- paste0("K",1:K)
cat(sprintf("Computation took %0.1f min.\n",r["elapsed"]/60))
rm(n,r)

# WRITE ESTIMATED ADMIXTURE PROPORTIONS TO FILE
# ---------------------------------------------
cat("Saving estimated admixture proportions to ",admix.pred.file,".\n",sep="")
write.admix.file(sprintf("%s/%s",data.dir,admix.pred.file),out$Q)
