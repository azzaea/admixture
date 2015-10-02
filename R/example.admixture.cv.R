# Find the best choice for the L0-penalty strength using
# cross-validation, in which we measure the error rate in the
# estimated of the left-out genotypes.
library(parallel)
library(Matrix)
source("misc.R")
source("mcmc.R")
source("admixture.R")
source("read.data.R")
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
n.cv     <- 8     # Number of trials for each choice of L0-penalty.
seed     <- 1     # Specifies the sequence of pseudorandom numbers.
mc.cores <- 10    # Number of CPUs to use.

# These variables specify: for each test individual, the number of
# populations contributing to the individual's genome (k.test);
# candidate values for L0-penalty strength (a); the sequence of
# inverse temperatures for simulated annealing such that the first
# temperature is 1 and the last temperature is near 0 (T).
k.test <- rep(c(1,2,4),each = 100)
a      <- 10^seq(-3,-1,0.25)
T      <- 1/cooling.sched.geom(1e4,6e-4)

# Pathname of the ADMIXTURE executable.
admix.exec <- "/DNAData2/bin/admixture"

# All the data files are written to this directory, and ADMIXTURE
# output files will be created in this directory.
data.dir <- "/DNAData2/tmp/peter/admixture"

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

# Define a function that takes as input a p x k matrix of
# population-specific allele frequencies (f), and a vector of
# admixture proportions of length k (q), where p is the number of
# markers and k is the number of ancestral populations, and outputs a
# vector of genotypes sampled i.i.d. according to the specified the
# allele frequencies and admixture poportions.
sample.genotypes <- function (f, q) {

  # Get the number of markers (p) and the number of ancestral
  # populations (K).
  p <- nrow(f)
  K <- ncol(f)
  
  # Independently for each marker, and for each of the two allele
  # copies, draw the population of origin.
  z1 <- sample(1:K,size = p,prob = q,replace = TRUE)
  z2 <- sample(1:K,size = p,prob = q,replace = TRUE)

  # For each marker, and for each of the two allele copies, draw the
  # allele, and output the (unphased) genotype, represented as an
  # allele count.
  return((runif(p) < f[cbind(1:p,z1)]) +
         (runif(p) < f[cbind(1:p,z2)]))
}

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

# FIT L0-PENALTY STRENGTH USING CROSS-VALIDATION
# ----------------------------------------------
# Initialize storage for the cross-validation statistics.
cat("Fitting L0-penalty strength using cross-validation.\n\n")
na  <- length(a)
err <- matrix(0,n.cv,na)

# Repeat for each cross-validation trial.
for (iter in 1:n.cv) {

  # Choose which genotypes to set to missing.
  nc   <- n + n.test
  r    <- sample(nc*p,floor(nc*p*prop.na))
  A    <- Matrix(FALSE,nc,p,sparse = TRUE)
  A[r] <- TRUE
  
  # Repeat for each candidate value of the L0-penalty strength.
  for (i in 1:na) {
    cat(sprintf("TRIAL=%d, LOG10a=%0.2f\n",iter,log10(a[i])))

    # Randomly set a proportion of the genotypes to missing.
    geno <- rbind(geno.train,geno.test)
    geno[which(A)] <- NA
    
    # Estimate the allele frequencies and admixture proportions in the
    # unlabeled examples.
    z <- c(q.train %*% 1:K,rep(NA,n.test))
    r <- system.time(out <-
           admixture.em(geno,K,z,e,a = a[i],F = out.admix$F,
                        Q = out.admix$Q,exact.q = FALSE,T = T,
                        tolerance = 1e-4,mc.cores = mc.cores))
    cat(sprintf("Computation took %0.1f min.\n",r["elapsed"]/60))

    # Calculate the expected absolute difference averaged over all
    # missing genotypes.
    geno        <- rbind(geno.train,geno.test)
    err[iter,i] <- calc.geno.error(geno,A,out$F,out$Q,e)
    cat(sprintf("Average error in missing genotypes is %0.4f.\n",err[iter,i]))
    cat("\n")
  }
}
rm(na,nc,iter,i,z,r,geno,A,out)

# SUMMARIZE RESULTS OF CROSS-VALIDATION
# -------------------------------------
cat("Average error in missing genotypes:\n")
r <- apply(err,2,function (x) c(quantile(x,0.1),mean(x),quantile(x,0.9)))
r <- cbind(log10(a),t(r))
colnames(r) <- c("log10a","Pr>.1","mean","Pr>.9")
print(round(r,digits = 3))
