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
nrep     <- 5     # Number of markers in each LD block.
K        <- 20    # Number of ancestral populations.
n        <- 100   # Number of single-origin training samples.
prop.na  <- 0.01  # Proportion of genotypes that are missing.
e        <- 0.01  # Probability of genotype error.
a        <- 3e-3  # L0-penalty strength.
seed     <- 1     # Specifies the sequence of pseudorandom numbers.
mc.cores <- 2     # Number of CPUs to use.

# A vector that specifies, for each test individual, the number of
# populations contributing to the individual's genome.
k.test <- rep(c(1,2,4),each = 100)

# The sequence of inverse temperatures for simulated annealing such
# that the first temperature is 1 and the last temperature is near 0.
T <- 1/cooling.sched.geom(1e4,6e-4)

# All the data files are written to this directory, and ADMIXTURE
# output files will be created in this directory.
data.dir <- "."

# Initialize the random number generator.
set.seed(seed)

# GENERATE TRAINING AND TEST DATA
# -------------------------------
cat("Simulating data with the following settings:\n")
cat("  markers              ",p,"x",nrep,"=",p*nrep,"\n")
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
p          <- p * nrep
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
  geno.train[i,] <- rep(sample.genotypes(f,q.train[i,]),each = nrep)
}

# Generate the test samples.
cat("Generating genotype data for test samples.\n")
for (i in 1:n.test) {

  # Randomly sample the admixture proportions,
  k             <- sample(1:K,size = k.test[i],prob = p.deme)
  q.test[i,k]   <- 1/length(k)
  geno.test[i,] <- rep(sample.genotypes(f,q.test[i,]),each = nrep)
}
rm(i,k)

# SAVE TRAINING AND TEST DATA TO FILE
# -----------------------------------
srcdir <- getwd()
setwd(data.dir)

# Save the genotypes of the training and test samples to a text file
# in the format used by EIGENSTRAT (a .geno file). It is important to
# reorder the training samples according to the population labels so
# that the admixture proportions are outputted in the same order as
# they are here.
cat("Writing genotype data to .geno file.\n")
z    <- c(q.train %*% 1:K)
rows <- order(z)
write.geno.file("sim.geno",rbind(geno.train[rows,],geno.test))

# To run ADMIXTURE in "supervised" mode, I need to create an
# additional .pop file.
cat("Writing population labels for training samples to .pop file.\n")
write.pop.file("sim.pop",c(z[rows],rep(NA,n.test)))
rm(z,rows)

# RUN ADMIXTURE
# -------------
# Run ADMIXTURE in "supervised" mode.
cat("Running ADMIXTURE.\n")
r <- system.time(system(sprintf(paste("%s --supervised --seed=1 -j%d",
                                      "sim.geno %d > sim.out"),
                                admix.exec,mc.cores,K)))
cat(sprintf("Computation took %0.1f s.\n",r["elapsed"]))

# Load the estimated allele frequencies and admixture proportions.
out.admix <- list(F = read.admixture.P.file(paste0("sim.",K,".P")),
                  Q = read.admixture.Q.file(paste0("sim.",K,".Q")))

# Restore the working directory.
setwd(srcdir)
rm(r)

# Randomly set a small proportion of the genotypes to missing.
cat("Setting ",100*prop.na,"% of the genotypes to NA.\n",sep="")
geno.train[runif(n*p) < prop.na]     <- NA
geno.test[runif(n.test*p) < prop.na] <- NA

# COMPUTE L0-PENALIZED ADMIXTURE ESTIMATES USING EM
# -------------------------------------------------
cat("Estimating allele frequencies and admixture proportions",
    "in unlabeled examples.\n")
cat("EM algorithm setings:\n")
cat("  Strength of L0-penalty    ",a,"\n")
cat("  Genotype error probability",e,"\n")
z <- c(q.train %*% 1:K,rep(NA,n.test))
r <- system.time(out <-
       admixture.em(rbind(geno.train,geno.test),K,z,e,a = a,F = out.admix$F,
                    Q = out.admix$Q,exact.q = FALSE,T = T,tolerance = 1e-4,
                    cg = TRUE,mc.cores = mc.cores))
cat(sprintf("Computation took %0.1f min.\n",r["elapsed"]/60))
rm(z,r)

# ASSESS ACCURACY IN TEST SAMPLES
# -------------------------------
cat("Overlap between estimated and ground-truth admixture proportions:\n")
bins <- c(seq(0,0.8,0.1),0.85,0.9,0.95,1)
r <- rbind(table(cut(rowSums(pmin(q.test,out.admix$Q[-(1:n),])),bins)),
           table(cut(rowSums(pmin(q.test,out$Q[-(1:n),])),bins)))
rownames(r) <- c("ADMIXTURE","EM + L0")
colnames(r) <- bins[-length(bins)]
print(r)
cat("\n")

cat("Number of contributing ancestral populations:\n")
r <- table(factor(rowSums(q.test > 0.001)),
           factor(rowSums(out.admix$Q[-(1:n),] > 0.001)))
names(dimnames(r)) <- c("true","ADMIXTURE")
print(r)
cat("\n")
r <- table(factor(rowSums(q.test > 0.001)),
           factor(rowSums(out$Q[-(1:n),] > 0.001)))
names(dimnames(r)) <- c("true","EM + L0")
print(r)
