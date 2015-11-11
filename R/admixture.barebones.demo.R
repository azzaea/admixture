# TO DO: Explain here what this script does.
library(bayesm)
source("misc.R")
source("sim.data.R")
source("admixture.barebones.R")

# These script parameters specify the number of unlinked genetic
# markers (p), the number of samples (n), the probability of a
# genotype error (e), and the amount of genetic drift for each of the
# K populations (d).
p <- 100
n <- 300
e <- 0.01
d <- rep(0.1,3)

# Initialize the random number generator.
set.seed(1)

# Generate the population allele frequencies. They are stored as a p x
# K matrix, where p is the number of markers (loci) and K is the
# number of populations.
af <- runif(p,min = 0.1,max = 0.9)
f  <- sample.af(af,d)
rm(af)

# Generate the genotype data. All individuals are either single-origin
# (one contributing population), or admixed between exactly two
# populations.
K    <- length(d)
Q    <- sample.admix.2way(n,K)
geno <- sample.genotype.matrix(f,Q)
rm(Q)

# Compute admixture estimates using the EM algorithm.
out <- admixture.em.barebones(geno,K)

# Reorder the columns so that the estimates best match the
# ground-truth admixture proportions.
