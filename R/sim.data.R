# SUMMARY
# -------
# This file contains functions for simulating data. Here is an
# overview of the functions defined in this file:
#
#   beta.mean.var(mu,s)
#   sample.af(af,d)
#   sample.admix.2way(n,K)
#   sample.genotypes(f,q)
#   sample.genotype.matrix(f,Q)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Return the parameters of the Beta distribution that yield a random
# variable with mean mu and variance s.
beta.mean.var <- function (mu, s) {
  a <- mu^2*((1-mu)/s - 1/mu)
  return(pmax(0,c(a,a*(1/mu - 1))))
}

# ----------------------------------------------------------------------
# Return a matrix of randomly generated population-specific allele
# frequencies f[i,k], in which f[i,k] is drawn from the Beta
# distribution with mean af[i] and variance d[k]*af[i]*(1-af[i]).
sample.af <- function (af, d) {

  # Get the number of populations.
  K <- length(d)

  # The allele frequencies are stored as a p x K matrix, where p is the
  # number of markers (loci) and K is the number of populations.
  f <- matrix(0,p,K)

  # Repeat for each locus, and for each population.
  for (i in 1:p)
    for (k in 1:K) {

      # Draw the population allele frequency from the Beta probability
      # distribution with mean mu = af and variance s = d*af*(1-af).
      x <- beta.mean.var(af[i],d[k]*af[i]*(1-af[i]))
      if (x[1] == 0 | x[2] == 0)
        stop("Invalid Beta distribution")
      f[i,k] <- rbeta(1,x[1],x[2])
    }

  return(f)
}

# ----------------------------------------------------------------------
# Return an n x K matrix of randomly generated admixture proportions
# for single-origin and 2-way admixed individuals.
sample.admix.2way <- function (n, K) {

  # The admixture proportions are stored as an n x K matrix, where n
  # is the number of samples and K is the number of ancestral
  # populations.
  Q <- matrix(0,n,K)

  # Repeat for each sample.
  for (i in 1:n) {

    # With the flip of a coin, simulate a single-origin individual or
    # a two-way admixed individual, in which contributing populations
    # are drawn uniformly at random. For the admixed individuals, the
    # contributing populations are drawn uniformly at random.
    if (runif(1) < 0.5)
      Q[i,sample(K,1)] <- 1
    else {
      x <- runif(1)
      Q[i,sample(K,2)] <- c(x,1-x)
    }
  }

  return(Q)
}

# ----------------------------------------------------------------------
# Define a function that takes as input a p x k matrix of population
# allele frequencies (f), and a vector of admixture proportions of
# length k (q), where p is the number of markers and k is the number
# of ancestral populations, and outputs a vector of genotypes sampled
# i.i.d. according to the specified the allele frequencies and
# admixture poportions.
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

# ----------------------------------------------------------------------
# Return an n x p matrix of genotypes drawn at random according to the
# population allele frequencies (f) and admixture proportions (Q) by
# calling function sample.genotypes.
sample.genotype.matrix <- function (f, Q) {

  # Get the number of samples (n) and the number of markers (p).
  n <- nrow(Q)
  p <- nrow(f)
  
  # Initialize the return value.
  geno <- matrix(0,n,p)

  # Repeat for each sample.
  for (i in 1:n)
    geno[i,] <- sample.genotypes(f,Q[i,])

  return(geno)
}
