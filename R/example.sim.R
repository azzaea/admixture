# Compare estimation of admixture proportions, with and without the L0
# penalty term, in simulated genotype data. In this case, all the
# samples are unlabeled.
suppressPackageStartupMessages({
library(parallel)
library(turboEM)
library(bayesm)
source("misc.R")
source("mcmc.R")
source("admixture.R")
source("sim.data.R")
dyn.load("mcmc.so")
dyn.load("admixture.so")})

# SCRIPT PARAMETERS
# -----------------
p        <- 2000  # Number of (unlinked) markers.
n        <- 500   # Number of samples.
e        <- 0.01  # Probability of genotype error.
a        <- 1e-3  # Strength of L0-penalty.
seed     <- 1     # Specifies the sequence of pseudorandom numbers.
mc.cores <- 20    # Number of CPUs to use.

# Specifies the amount of genetic drift for each population. For
# details on the relationship between this parameter and population
# demography under the neutral Wright–Fisher diffusion model, see:
# Nicholson et al (2002), "Assessing population differentiation and
# isolation from single-nucleotide polymorphism data," Journal of the
# Royal Statistical Society, Series B 64(4): 695–715.
d <- rep(0.1,10)

# The sequence of inverse temperatures for simulated annealing such
# that the first temperature is 1 and the last temperature is near 0.
T <- 1/cooling.sched.geom(1e4,6e-4)

# Initialize the random number generator.
set.seed(seed)

# LOAD ALLELE FREQUENCY DATA
# --------------------------
# Load the allele frequencies estimated from the HGDP data.
cat("Loading allele frequency data.\n")
af <- read.table("../data/hgdp.af.txt",sep = " ",stringsAsFactors = FALSE,
                 header = TRUE)$f1

# Randomly select a subset of the allele frequencies >5%; these are
# the allele frequencies of all loci in the ancestral population (af =
# ancestral frequency).
af <- af[af > 0.05]
af <- af[sample(length(af),p)]

# GENERATE POPULATION ALLELE FREQUENCIES
# --------------------------------------
# The allele frequencies are stored as a p x K matrix, where p is the
# number of markers (loci) and K is the number of populations.
cat("Generating population allele frequencies.\n")
f <- sample.af(af,d)

# GENERATE DATA SET
# -----------------
# In this data set, all individuals are either single-origin (one
# contributing population), or admixed between exactly two
# populations.
cat("Generating data.\n")
K        <- length(d)
Q        <- sample.admix.2way(n,K)
sim.data <- list(Q = Q,geno = sample.genotype.matrix(f,Q))
rm(Q)

# Remove any markers that show no variation in the sample.
r             <- colMeans(sim.data$geno)/2
markers       <- which(r > 0 & r < 1)
sim.data$geno <- sim.data$geno[,markers]
rm(r,markers)

# COMPUTE ADMIXTURE ESTIMATES USING TURBOEM
# -----------------------------------------
cat("Fitting admixture model to data.\n")
r <- system.time(out.em <- admixture.em(sim.data$geno,K,e = e,
                                        mc.cores = mc.cores))
with(out.em,
     cat(sprintf("Turbo-EM completed after %d iterations and %0.1f min.\n",
                 length(loglikelihood),r["elapsed"]/60)))
rm(r)

# Reorder the columns (ancestral populations) so that they best match
# the ground-truth admixture proportions.
cols <- rep(NA,K)
for (i in 1:K)
  cols[i] <- which.min(colSums(abs(out.em$Q - sim.data$Q[,i])))
out.em$Q <- out.em$Q[,cols]
out.em$F <- out.em$F[,cols]
rm(cols,i)

# COMPUTE L0-PENALIZED ADMIXTURE ESTIMATES USING EM
# -------------------------------------------------
cat("Fitting L0-penalized admixture model to data.\n")
r <- system.time(out.sparse <-
       admixture.em(sim.data$geno,K,e = e,a = a,exact.q = FALSE,T = T,
                    F = out.em$F,Q = out.em$Q,tol = 1e-4,mc.cores = mc.cores))
with(out.sparse,
     cat(sprintf("Turbo-EM completed after %d iterations and %0.1f min.\n\n",
         length(loglikelihood),r["elapsed"]/60)))
rm(r)

# SUMMARIZE ACCURACY OF ADMIXTURE ESTIMATES
# -----------------------------------------
# Create a table summarizing the error in the estimated admixture
# proportions.
cat("Overlap between estimated and ground-truth admixture proportions:\n")
bins <- c(seq(0,0.8,0.1),0.85,0.9,0.95,1)
r    <- rbind(table(cut(rowSums(pmin(sim.data$Q,out.em$Q)),bins)),
              table(cut(rowSums(pmin(sim.data$Q,out.sparse$Q)),bins)))
rownames(r) <- c("ML","L0")
colnames(r) <- bins[-length(bins)]
print(r)
cat("\n")

# Print a table summarizing the number of contributing ancestral
# populations (>1%).
cat("Number of contributing ancestral populations (>1%):\n")
r <- rbind(summary(factor(rowSums(sim.data$Q > 0.01),1:K)),
           summary(factor(rowSums(out.em$Q > 0.01),1:K)),
           summary(factor(rowSums(out.sparse$Q > 0.01),1:K)))
rownames(r) <- c("true","ML","L0")
print(r)
