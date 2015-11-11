# SUMMARY
# -------
# blah, blah, blah. Here is an overview of the functions defined in
# this file:
#
#   <function does here>
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns the probability of the unphased genotype (represented by the
# allele count) given the phased genotype (u1, u2). The last input
# argument (e) specifies the probability of a genotype error. Set
# input u = u1 + u2, where (u1,u2) is the phased genotype. This
# function does not allow for the inputs to have missing values (NA).
genoprob <- function (x, u, e)
  (1-2*e) * (x == u) + e * (x != u)

# ----------------------------------------------------------------------
# TO DO: Describe this function here.
admixture.Estep.slow <- function (X, F, Q, e) {

  # This is a small constant added to counts to ensure that they are
  # always greater than zero.
  eps <- 1e-6
  
  # Get the number of samples (n), the number of markers (p), and the
  # number of ancestral populations (K).
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(F)
  
  # Initialize the counts.
  m  <- matrix(eps,n,K)
  n0 <- matrix(eps,p,K)
  n1 <- matrix(eps,p,K)
  
  # Initialize storage for the posterior probabilities of the hidden
  # (phased) genotypes x population indicators for a single sample.
  r <- array(0,dim = c(p,4,K,K))
    
  # Update the expected allele counts and population counts. Repeat
  # for each sample.
  for (i in 1:n) {
    colnames(r) <- c("00","01","10","11")

    # Repeat for each combination of ancestral populations.
    for (j in 1:K)
      for (k in 1:K) {

        # Compute the posterior probabilities for all possible hidden
        # (phased) genotype configurations at each locus, given the
        # assigment (j,k) to the population-of-origin indicators. Note
        # that genotype configurations (0,1) and (1,0) do *not* have
        # the same posterior probability here (except when j and k are
        # the same).
        r[,"00",j,k] <- genoprob(X[i,],0,e) * (1 - F[,j]) * (1 - F[,k])
        r[,"01",j,k] <- genoprob(X[i,],1,e) * (1 - F[,j]) * F[,k]
        r[,"10",j,k] <- genoprob(X[i,],1,e) * F[,j]       * (1 - F[,k])
        r[,"11",j,k] <- genoprob(X[i,],2,e) * F[,j]       * F[,k]
        r[,,j,k]     <- r[,,j,k] * Q[i,j] * Q[i,k]
      }

    # Normalize the posterior probabilities.
    dim(r) <- c(p,4*K^2)
    r      <- r/rowSums(r)
      
    # Add the posterior probabilities to the sufficient statistics.
    dim(r)      <- c(p,4,K,K)
    colnames(r) <- c("00","01","10","11")
    m[i,]       <- m[i,] + apply(r,3,sum) + apply(r,4,sum)
    for (k in 1:K) {
      n0[,k] <- n0[,k] + rowSums(drop(r[,"00",k,])) +
                         rowSums(drop(r[,"01",k,])) +
                         rowSums(drop(r[,"00",,k])) +
                         rowSums(drop(r[,"10",,k]))
      n1[,k] <- n1[,k] + rowSums(drop(r[,"10",k,])) +
                         rowSums(drop(r[,"11",k,])) +
                         rowSums(drop(r[,"01",,k])) +
                         rowSums(drop(r[,"11",,k]))
    }
  }

  # Return a list containing the expected allele counts (n0 and n1)
  # and the expected population counts (m).
  return(list(n0 = n0,n1 = n1,m = m))
}

# ----------------------------------------------------------------------
# TO DO: Describe function here.
admixture.em.barebones <-
  function (X, K, e = 0.001, tolerance = 1e-4, max.iter = 1e3) {

  # Get the number of samples (n) and the number of markers (p).
  n <- nrow(X)
  p <- ncol(X)

  # Randomly initialize the p x K matrix of allele frequencies, and
  # initialize the n x K matrix of admixture proportions.
  F <- matrix(runif(p*K),p,K)
  Q <- matrix(1/K,n,K)
  
  # Repeat until convergence criterion is met.
  for (iter in 1:max.iter) {

    # Save the current parameter estimates.
    F0 <- F
    Q0 <- Q
    
    # E-STEP
    # Compute the expected allele counts and the expected population
    # counts.
    out <- admixture.Estep.slow(X,F,Q,e)
    
    # M-STEP
    # Adjust the allele frequencies and admixture proportions using
    # the standard M-step updates.
    F <- with(out,n1/(n0 + n1))
    Q <- with(out,m/rowSums(m))
  
    # CHECK CONVERGENCE
    # Print the status of the algorithm and check convergence.
    # Convergence is reached when the maximum absolute difference
    # between the parameters at two successive iterations is less than
    # the specified tolerance.
    err <- list(f = max(abs(F - F0)),
                q = max(abs(Q - Q0)))
    caterase(sprintf("iter=%d dF=%0.1e dQ=%0.1e",iter,max(err$f),max(err$q)))
    if (max(err$f,err$q) < tolerance)
      break
  }
  cat("\n")
  
  # Return a list containing the estimated allele frequencies and
  # admixture proportions.
  return(list(F = F,Q = Q))
}
