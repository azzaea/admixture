# Illustrates how to use the EM algorithm (admixture.em) to predict
# admixture proportions when we have a reference set of labeled,
# single-origin samples.
library(parallel)
source("misc.R")
source("mcmc.R")
source("admixture.R")
source("sim.data.R")
dyn.load("mcmc.so")
dyn.load("admixture.so")

# SCRIPT PARAMETERS
# -----------------
p        <- 200   # Number of (unlinked) markers.
K        <- 20    # Number of ancestral populations.
n        <- 100   # Number of single-origin training samples.
prop.na  <- 0.01  # Proportion of genotypes that are missing.
e        <- 0.01  # Probability of genotype error.
a        <- 3e-3  # L0-penalty strength.
seed     <- 1     # Specifies the sequence of pseudorandom numbers.
mc.cores <- 20    # Number of CPUs to use.

# A vector that specifies, for each test individual, the number of
# populations contributing to the individual's genome.
k.test <- rep(c(1,2,4),each = 100)

# The sequence of inverse temperatures for simulated annealing such
# that the first temperature is 1 and the last temperature is near 0.
T <- 1/cooling.sched.geom(1e3,6e-3)

# All the data files are written to this directory.
data.dir <- "."

# Initialize the random number generator.
set.seed(seed)

# GENERATE TRAINING AND TEST DATA
# -------------------------------
cat("Simulating data with the following settings:\n")
cat("  markers              ",p,"\n")
cat("  ancestral populations",K,"\n")
cat("  training samples     ",n,"\n")
cat("  test samples         ",length(k.test),"\n")
cat("Generating genotype data for single-origin training samples.\n")

# Generate, for each ancestral population, the total proportion of
# chromosomes that are attributed to that population.
p.deme <- runif(K,0.2,0.8)
p.deme <- p.deme / sum(p.deme)

# Generate, for each ancestral population, the allele frequencies at
# at all the markers, assuming the markers are unlinked.
f <- runif(K*p)
f <- matrix(f,p,K)

# The admixture proportions are stored as an n x k matrix, and the
# genotypes are stored as an n x p matrix, where n is the number of
# samples, p is the number of markers, and k is the number of
# ancestral populations.
n.test     <- length(k.test)
q.train    <- matrix(0,n,K)
q.test     <- matrix(0,n.test,K)
geno.train <- matrix(0,n,p)
geno.test  <- matrix(0,n.test,p)
  
# Generate the single-origin training samples.
for (i in 1:n) {

  # Randomly sample the ancestral population of origin, then sample
  # the genotypes given the population-specific allele frequencies and
  # the ancestral population of origin.
  k              <- sample(1:K,size = 1,prob = p.deme)
  q.train[i,k]   <- 1
  geno.train[i,] <- sample.genotypes(f,q.train[i,])
}

# Generate the test samples.
cat("Generating genotype data for test samples.\n")
for (i in 1:n.test) {

  # Randomly sample the admixture proportions,
  k             <- sample(1:K,size = k.test[i],prob = p.deme)
  q.test[i,k]   <- 1/length(k)
  geno.test[i,] <- sample.genotypes(f,q.test[i,])
}
rm(i,k)

# Randomly set a small proportion of the genotypes to missing.
cat("Setting ",100*prop.na,"% of the genotypes to NA.\n",sep="")
geno.train[runif(n*p) < prop.na]     <- NA
geno.test[runif(n.test*p) < prop.na] <- NA

# COMPUTE MAXIMUM-LIKELIHOOD ADMIXTURE ESTIMATES USING EM
# -------------------------------------------------------
cat("Computing maximum-likelihood admixture proportion estimates.\n")
X <- rbind(geno.train,geno.test)
z <- c(q.train %*% 1:K,rep(NA,n.test))
r <- system.time(out.em <-
       admixture.em(X,K,z,e = e,cg = TRUE,mc.cores = mc.cores))
cat(sprintf("Computation took %0.1f min.\n",r["elapsed"]/60))
rm(r)

# COMPUTE L0-PENALIZED ADMIXTURE ESTIMATES USING EM
# -------------------------------------------------
cat("Computing L0-penalized admixture proportion estimates.\n")
r <- system.time(out.sparse <-
       admixture.em(X,K,e = e,a = a,F = out.em$F,Q = out.em$Q,cg = TRUE,
                    exact.q = FALSE,T = T,mc.cores = mc.cores))
cat(sprintf("Computation took %0.1f min.\n",r["elapsed"]/60))
rm(r)

# ASSESS ACCURACY IN TEST SAMPLES
# -------------------------------
cat("Overlap between estimated and ground-truth admixture proportions:\n")
bins <- c(seq(0,0.8,0.1),0.85,0.9,0.95,1)
r <- rbind(table(cut(rowSums(pmin(q.test,out.em$Q[-(1:n),])),bins)),
           table(cut(rowSums(pmin(q.test,out.sparse$Q[-(1:n),])),bins)))
rownames(r) <- c("ML","L0")
colnames(r) <- bins[-length(bins)]
print(r)
cat("\n")

cat("Number of contributing ancestral populations (>1%):\n")
r <- table(factor(rowSums(q.test > 0.01)),
           factor(rowSums(out.em$Q[-(1:n),] > 0.01)))
names(dimnames(r)) <- c("true","ML")
print(r)
cat("\n")
r <- table(factor(rowSums(q.test > 0.01)),
           factor(rowSums(out.sparse$Q[-(1:n),] > 0.01)))
names(dimnames(r)) <- c("true","L0")
print(r)
