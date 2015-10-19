# SUMMARY
# -------
# This file contains various functions implementing estimation of
# admixture proportions from genotype data. Here is an overview of the
# functions defined in this file:
#
#   genoprob.given.q.fast(f,q,e)
#   calc.geno.error(X,A,F,Q,e)
#   update.q.sparse.exact(M,a)
#   update.q.sparse.approx(M,x0,a,T,u)
#   update.q.sparse.approx.mc(M,x0,a,T,u,mc.cores)
#   admixture.labeled.Estep.fast(X,F,z,e)
#   admixture.labeled.Estep.mc(X,F,z,e,mc.cores)
#   admixture.Estep.fast(X,F,Q,n0,n1,e)
#   admixture.Estep.mc(X,F,Q,n0,n1,e,mc.cores)
#   admixture.em(X,n0,n1,e,a,F,Q,tolerance,max.iter,exact.q,mc.cores,verbose,T)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns the probability of the genotype (i.e. allele count) given
# the vector of population-specific allele frequencies (f), the
# admixture proportions (q) and the genotype error probability (e).
# The return value is a vector of length 3 containing the probability
# that the allele count is 0, 1 and 2.
#
# This function calls "genoprob_given_q_Call", a
# function compiled from C code, using the .Call interface. For more
# details on how to load the C function into R, see the comments
# accompanying function admixture.labeled.Estep.fast.
genoprob.given.q.fast <- function (f, q, e) {

  # Initialize the return value.
  px <- rep(0,3)
  
  # Execute the C routine using the .Call interface. The only
  # component that changes is px.
  out <- .Call("genoprob_given_q_Call",
               f  = as.double(f),  # Allele frequencies.
               q  = as.double(q),  # Admixture proportions.
               e  = as.double(e),  # Genotype error probability.
               px = px)            # Genotype probabilities.
  return(px)
}

# ----------------------------------------------------------------------
# The inputs are as follows: X, the n x p genotype matrix, where n is
# the number of samples and p is the number of markers; A, a n x p
# logical matrix, that indicates which of the genotypes to use in
# calculating the genotype error (TRUE means that the corresponding
# genotype is included in the error calculation; FALSE, means that the
# corresponding genotype is not included in the error calculation); F,
# the matrix of population-specific allele frequencies; Q, the matrix
# of estimated admixture proportions; and e, the probability of a
# genotyping error. For more details about inputs F and Q, see
# function admixture.em.
#
# The return value is the expected absolute difference between the
# ground-truth and estimated genotypes given the allele frequencies
# and admixture proportions, averaged over all genotypes (i,j) such
# that A[i,j] = 1.
calc.geno.error <- function (X, A, F, Q, e) {

  # Get the number of samples.
  n <- nrow(X)

  # Initialize the return value.
  err <- 0
  
  # Repeat for each sample.
  for (i in 1:n) {
    
    # Get the loci corresponding to unobserved genotypes.
    markers <- which(A[i,])

    # Repeat for each unobserved genotype.
    for (j in markers) {

      # Compute the genotype probabilities given the estimated
      # admixture proportions and allele frequencies.
      px <- genoprob.given.q.fast(F[j,],Q[i,],e)
        
      # Compute the expected absolute difference between the estimated
      # and ground-truth genotype.
      err <- err + sum(px * abs(X[i,j] - 0:2))
    }
  }
  
  # Return the average error.
  return(err/sum(A))
}

# ----------------------------------------------------------------------
# Computes L0-penalized admixture proportion estimates given the
# expected population counts. Input M is an n x k matrix of expected
# population counts, where n is the number of samples and k is the
# number of ancestral populations. Input a specifies the strength of
# the L0-penalty.
update.q.sparse.exact <- function (M, a) {

  # Get the number of samples (n) and the number of ancestral populations (K).
  n <- nrow(M)
  K <- ncol(M)
  
  # Initialize the return value.
  Q <- matrix(0,n,K)
  
  # Generate all possible ways we can select K coordinates
  # (i.e. admixture proportions) that are nonzero, except for the case
  # when no coordinates are nonzero.
  X <- as.matrix(expand.grid(rep(list(0:1),K))[-1,])
  dimnames(X) <- NULL
  
  # Repeat for each sample.
  for (i in 1:n) {

    # Find the combination of nonzero admixture proportions which
    # yields the largest L0-penalized log-likelihood.
    x <- X[which.max(apply(X,1,function (x) sparse.multinom.logp(M[i,],x,a))),]

    # Store the final L0-penalized estimate.
    Q[i,] <- sparse.multinom.ml(M[i,],x)
  }

  # Return the updated admixture proportions.
  return(Q)
}

# ----------------------------------------------------------------------
# This is the same as function update.q.sparse.exact, except that it
# computes an approximate solution using simulated annealing. The two
# additional inputs are the inverse temperature schedule (T) and the
# sequence of pseudorandom numbers used to simulate the Markov chain
# (u). For more details on the last input argument, see function
# mh.sparse.multinom.fast.
update.q.sparse.approx <- function (M, x0, a, T, u) {

  # Add this number to the admixture proportions to ensure that we
  # never try compute log(0).
  e <- 1e-6

  # Get the number of samples (n) and the number of ancestral populations (K).
  n <- nrow(M)
  K <- ncol(M)
  
  # Initialize the return value.
  Q <- matrix(0,n,K)
  
  # Repeat for each sample.
  for (i in 1:n) {
    
    # Find the combination of nonzero admixture proportions which
    # yields the largest L0-penalized log-likelihood.
    x <- mh.sparse.multinom.fast(as.double(x0[i,]),M[i,],a,T,u)

    # Store the final L0-penalized estimate.
    Q[i,] <- sparse.multinom.ml(M[i,],x)
  }
  
  # Return the updated admixture proportions.
  return(Q)
}

# ----------------------------------------------------------------------
# This is the multicore variant of update.q.sparse.approx.
update.q.sparse.approx.mc <- function (M, x0, a, T, u, mc.cores = 2) {

  # Assign each sample to a CPU, and compute admixture proportion
  # estimates for each set of samples.
  rows <- distribute(1:nrow(M),mc.cores)
  out  <- mclapply(rows,function(i) update.q.sparse.approx(M[i,],x0[i,],a,T,u),
                   mc.cores = mc.cores)

  # Aggregate the outputs from the individual CPUs.
  Q <- do.call(rbind,out)
  Q[unlist(rows),] <- Q
  return(Q)
}

# ----------------------------------------------------------------------
# Compute the expected allele counts "n0" and "n1", the sufficient
# statistics for updating the population-specific allele frequencies.
# These sufficient statistics are stored in two p x k matrices, where
# p is the number of markers and k is the number of ancestral
# populations.
#
# This function calls "admixture_labeled_Estep_Call", a function
# compiled from C code, using the .Call interface. To load the C
# function into R, first build the "shared object" (.so) file using
# command
#
#   R CMD SHLIB admixture.c
#
# Next, load the shared objects into R using the R function dyn.load:
#
#   dyn.load("admixture.so")
#
admixture.labeled.Estep.fast <- function (X, F, z, e) {

  # Get the number of markers (p) and the number of ancestral
  # populations (K).
  p <- ncol(X)
  K <- ncol(F)
  
  # Check that genotype and allele frequency matrices are in double
  # precision.
  if (!is.double(X))
    stop("Input matrix X must be in double precision.")
  if (!is.double(F))
    stop("Input matrix F must be in double precision.")

  # Initialize the expected allele counts.
  n0 <- matrix(0,p,K)
  n1 <- matrix(0,p,K)
  
  # Execute the C routine using the .Call interface, and return the
  # expected allele counts in a list object. The main reason for using
  # .Call interface is that there is less of a constraint on the size
  # of the input matrices. The only components that change are n0 and
  # n1. Note that I need to subtract 1 from the indices because R
  # vectors start at 1, and C arrays start at 0.
  out <- .Call("admixture_labeled_Estep_Call",
               X  = X,              # Genotype matrix.
               F  = F,              # Allele frequency matrix.
               z  = as.double(z-1), # Ancestral population labels.
               e  = as.double(e),   # Genotype error probability.
               n0 = n0,             # Expected counts of "0" allele.
               n1 = n1)             # Expected counts of "1" allele.
  return(list(n0 = n0,n1 = n1))
}

# ----------------------------------------------------------------------
# This is the multicore variant of admixture.labeled.Estep.fast.
admixture.labeled.Estep.mc <- function (X, F, z, e, mc.cores = 2) {
  
  # Assign each sample to a CPU, and compute the expected allele
  # counts for each set of samples.
  out <- mclapply(distribute(1:nrow(X),mc.cores),
                  function (i) admixture.labeled.Estep.fast(X[i,],F,z[i],e),
                  mc.cores = mc.cores)
  
  # Aggregate the outputs from the individual CPUs.
  return(list(n0 = Reduce('+',lapply(out,function (x) x$n0)),
              n1 = Reduce('+',lapply(out,function (x) x$n1))))
}

# ----------------------------------------------------------------------
# Compute the expected allele counts "n0" and "n1" and the expected
# population counts "m", the sufficient statistics for updating the
# population-specific allele frequencies and the admixture
# proportions, respectively. Statistics n0 and n1 are stored in two p
# x k matrices, where p is the number of markers and k is the number
# of ancestral populations. The m counts are stored in an n x k
# matrix, where n is the number of samples.
#
# The expected allele counts are also input arguments. This allows the
# user to specify "prior counts" based on other information (e.g. a
# reference panel). If no prior counts are available, set all the
# entries of these two matrices to zero.
#
# This function calls "admixture_unlabeled_Estep_Call", a function
# compiled from C code, using the .Call interface. For more details on
# how to load the C function into R, see the comments accompanying
# function admixture.labeled.Estep.fast.
admixture.Estep.fast <- function (X, F, Q, n0, n1, e) {

  # Get the number of samples (n) and the number of ancestral
  # populations (K).
  n <- nrow(Q)
  K <- ncol(Q)

  # Check that genotype matrix, allele frequency matrix and admixture
  # proportions matrix are all in double precision.
  if (!is.double(X))
    stop("Input matrix X must be in double precision.")
  if (!is.double(F))
    stop("Input matrix F must be in double precision.")
  if (!is.double(Q))
    stop("Input matrix Q must be in double precision.")

  # Set the storage mode of inputs n0 and n1 to double precision.
  storage.mode(n0) <- "double"
  storage.mode(n1) <- "double"  
    
  # Initialize the expected population counts.
  m <- matrix(0,n,K)

  # Initialize storage for the posterior probabilities of the hidden
  # (phased) genotypes x population indicators.
  r00 <- matrix(0,K,K)
  r01 <- matrix(0,K,K)
  r10 <- matrix(0,K,K)
  r11 <- matrix(0,K,K)

  # Execute the C routine using the .Call interface, and return the
  # sufficient statistics in a list object. The main reason for using
  # .Call interface is that there is less of a constraint on the size
  # of the input matrices. The only components that change are m, n0
  # and n1. Note that I need to subtract 1 from the indices because R
  # vectors start at 1, and C arrays start at 0.
  out <- .Call("admixture_unlabeled_Estep_Call",
               X   = X,             # Genotype matrix.
               F   = F,             # Allele frequency matrix.
               Q   = Q,             # Admixture proportions matrix.
               e   = as.double(e),  # Genotype error probability.
               m   = m,             # Expected population counts.
               n0  = n0,            # Expected counts of "0" allele.
               n1  = n1,            # Expected counts of "1" allele.
               r00 = r00,           # Posterior probabilities.
               r01 = r01,
               r10 = r10,
               r11 = r11)
  return(list(m = m,n0 = n0,n1 = n1))
}

# ----------------------------------------------------------------------
# This is the multicore variant of admixture.Estep.fast.
admixture.Estep.mc <- function (X, F, Q, n0, n1, e, mc.cores = 2) {

  # Get the number of markers (p) and the number of ancestral
  # populations (K).
  p <- nrow(F)
  K <- ncol(F)
  
  # Assign each sample to a CPU, and compute the expected allele
  # counts for each set of samples.
  N    <- matrix(0,p,K)
  rows <- distribute(1:nrow(X),mc.cores)
  out  <- mclapply(rows,function (i) admixture.Estep.fast(X[i,],F,Q[i,],N,N,e),
                   mc.cores = mc.cores)

  # Aggregate the outputs from the individual CPUs. For the expected
  # allele counts (n0 and n1), also add the prior expected counts to
  # obtain the final result.
  n0 <- n0 + Reduce('+',lapply(out,function (x) x$n0))
  n1 <- n1 + Reduce('+',lapply(out,function (x) x$n1))
  m  <- do.call(rbind,lapply(out,function (x) x$m))
  m[unlist(rows),] <- m
  return(list(m = m,n0 = n0,n1 = n1))
}

# ----------------------------------------------------------------------
# Estimate population-specific allele frequencies and admixture
# proportions in unlabeled samples from genotypes. The non-optional
# inputs are as follows:
#
#   X   n x p genotype matrix, where n is the number of
#       samples and p is the number of markers;
#
#   K   number of ancestral populations;
#
#   z   vector giving the population of origin for each of the samples
#       (an integer between 1 and K), or NA if the sample is unlabeled.
#       If set to NULL, or not specified, all samples are unlabeled.
#
# The return value is a list with two list elements: F, the p x k
# matrix of population-specific allele frequency estimates; and Q, the
# n x k matrix of estimated admixture proportions, in which each row
# of Q sums to 1. For labeled samples, the admixture proportions are
# Q[i,k] = 1 when z[i] = k, otherwise all the other entries are
# exactly zero.
#
# There are two variations to the M-step update for the Q matrix. When
# the number of ancestral populations is small (k < 20), it is
# feasible to compute the L0-penalized estimate exactly by
# exhaustively calculating the posterior probability for each possible
# choice of the nonzero admixture proportions. Setting exact.q = TRUE
# will activate this option. However, for larger k, it is not feasible
# to compute the exact solution because the number of ways of choosing
# nonzero admixture proportions is too large. Instead, setting exact.q
# = FALSE computes an approximate solution using a simulated annealing
# algorithm. In this case, it is necessary to set input T. For an
# explanation of input T, see function update.q.sparse.approx.
#
# The cg parameter specifies the M-step update for the F matrix. When
# cg = FALSE, the binomial success rates are updated using the
# standard M-step solution that is derived by finding the roots of the
# partial derivatives of the expected complete log-likelihood. When cg
# = TRUE, the standard M-step update is adjusted using the conjugate
# gradient algorithm (specifically, using the Hestenes-Stiefel update
# formula). In some cases, I've found that the conjugate gradient
# upgrade can lead to improvements in the convergence rate of the EM
# iterates.
#
admixture.em <- function (X, K, z = NULL, e = 0.001, a = 0, F = NULL, Q = NULL,
                          tolerance = 1e-4, max.iter = 1000, exact.q = FALSE,
                          cg = FALSE,mc.cores = 1, verbose = TRUE, T = 1) {

  # This is a small constant added to allele frequency estimates to
  # ensure that the frequency estimates are always greater than zero.
  eps <- 0.001  

  # Get the number of samples (n) and the number of markers (p).
  n <- nrow(X)
  p <- ncol(X)

  # Following the "elastic net" approach, I multiply the penalty
  # strength by the sample size, which is equivalent to dividing the
  # log-likelihood by the sample size. In this case, the sample size
  # is 2 x n x p, where n is the number of individuals and p is the
  # number of markers. The 2 here accounts for the fact that there are
  # two allele copies per locus.
  a <- 2*n*p*a

  # If the population-of-origin labels are unspecified, set z to a
  # vector in which all the entries are missing.
  if (is.null(z))
    z <- rep(NA,n)
  
  # Get the set of samples that are labeled (i) and unlabeled (j).
  i <- which(!is.na(z))
  j <- which(is.na(z))
  z <- z[i]

  # Initialize the p x K matrix of allele frequencies.
  if (is.null(F))
    F <- matrix(runif(p*K),p,K)

  # Initialize the n x K matrix of admixture proportions. Acquire the
  # ground-truth admixture proportions for the labeled samples (i),
  # and set the admixture proportions for the unlabeled samples (j) to
  # be uniform.
  if (is.null(Q)) {
    Q             <- matrix(0,n,K)
    Q[cbind(i,z)] <- 1
    Q[j,]         <- 1/K
  }
  
  # For the approximate M-step update, generate a sequence of
  # pseudorandom numbers of length 4*n, where n is the length of the
  # inverse temperature schedule.
  if (a > 0 & !exact.q)
    u <- runif(4*length(T))
  
  # Repeat until convergence criterion is met.
  if (verbose)
    cat("iter delta-F delta-Q -beta-\n")
  g <- 0
  d <- 0
  for (iter in 1:max.iter) {

    # Save the current parameter estimates.
    F0 <- F
    Q0 <- Q

    # Save the gradient and search direction from the previous
    # iteration.
    g0 <- g
    d0 <- d

    # E-STEP
    # Compute the expected allele counts from the labeled, single-origin
    # samples. I add a small constant to n0 and n1 so that the counts are
    # never exactly zero.  
    if (length(i) == 0) {
      n0 <- matrix(eps,p,K)
      n1 <- matrix(eps,p,K)
    } else {
      if (mc.cores == 1)
        out <- admixture.labeled.Estep.fast(X[i,],F,z,e)
      else
        out <- admixture.labeled.Estep.mc(X[i,],F,z,e,mc.cores)
      n0 <- out$n0 + eps
      n1 <- out$n1 + eps
      rm(out)
    }

    # Compute the expected allele counts and the expected population
    # counts in the unlabeled samples only.
    if (mc.cores == 1)
      out <- admixture.Estep.fast(X[j,],F,Q[j,],n0,n1,e)
    else
      out <- admixture.Estep.mc(X[j,],F,Q[j,],n0,n1,e,mc.cores)
    M  <- out$m
    n0 <- out$n0
    n1 <- out$n1
    rm(out)
    
    # M-STEP
    if (cg) {

      # For all iterations after the first, compute the conjugate
      # gradient update for the binomial success rates. This update is
      # based on the Hestenes-Stiefel formula as described in Nocedal
      # and Wright (2006). Another option that doesn't seem to work as
      # well, but is included here for posterity, is the Polak-Ribiere
      # formula:
      #
      #   b = max(0,sum(g*(g - g0))/sum(g0^2))
      #
      # Note #1: The gradient is not the "true" gradient, but rather
      # the search direction implicit in the standard M-step
      # update. Note #2: I project the iterate onto [0,1] to ensure
      # that each of the parameters represents a proportion.
      g <- n1/(n0 + n1) - F
      if (iter == 1)
        b <- 0
      else
        b <- min(1,max(-1,sum(g*(g - g0))/sum((g - g0)*d0)))
      d   <- g + b*d0
      F   <- F + d
      F[] <- pmax(F,0.0001)
      F[] <- pmin(F,0.9999)
    }
    else {

      # Adjust the allele frequencies using the standard M-step
      # update.
      b <- NA
      F <- n1/(n0 + n1)
    }

    # Update the admixture proportions in the unlabeled samples.
    if (a == 0)
      Q[j,] <- M/rowSums(M)
    else if (exact.q)
      Q[j,] <- update.q.sparse.exact(M,a)
    else if (mc.cores == 1)
      Q[j,] <- update.q.sparse.approx(M,Q0[j,] > 0,a,T,u)
    else
      Q[j,] <- update.q.sparse.approx.mc(M,Q0[j,] > 0,a,T,u,mc.cores)
  
    # CHECK CONVERGENCE
    # Print the status of the algorithm and check convergence.
    # Convergence is reached when the maximum absolute difference
    # between the parameters at two successive iterations is less than
    # the specified tolerance.
    err <- list(f = max(abs(F - F0)),
                q = max(abs(Q - Q0)))
    if (verbose)
      caterase(sprintf("%4d %0.1e %0.1e %06.3f",iter,err$f,err$q,b))
    if (max(err$f,err$q) < tolerance)
      break
  }
  if (verbose)
    cat("\n")

  # Return a list containing the estimated allele frequencies and
  # admixture proportions.
  return(list(F = F,Q = Q))
}
