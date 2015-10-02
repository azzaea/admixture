# SUMMARY
# -------

# This file contains various functions to implement the
# Metropolis-Hastings algorithm for sampling from a binary state
# space. Here is an overview of the functions defined in this file:
#
#   sample.set(x)
#   sparse.multinom.ml(n,x)
#   sparse.multinom.logp(n,x,a)
#   cooling.sched.geom(n,p)
#   mh.bin(x0,f,T,samples.out,verbose)
#   mh.sparse.multinom.fast(x0,counts,a,T,u,samples.out)
# 
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns a single element uniformly at random from the set x.
sample.set <- function (x) {
  if (length(x) == 1)
    return(x)
  else
    return(sample(x,1))
}

# ----------------------------------------------------------------------
# Get the maximum-likelihood (ML) estimate of the nonzero multinomial
# probabilities given the category counts (n) and binary indicator for
# the nonzero coordinates (x). In the special case when the counts
# corresponding to the nonzero coordinates are all zero, I set the ML
# estimate to the uniform probability table.
sparse.multinom.ml <- function (n, x) {
  x <- as.logical(x)
  q <- rep(0,length(x))
  if (all(n[x] == 0))
    q[x] <- 1/sum(x)
  else
    q[x] <- n[x]/sum(n[x])
  return(q)
}

# ----------------------------------------------------------------------
# Computes the L0-penalized log-likelihood given the category counts
# (n) and a binary indicator for the nonzero coordinates (x).
sparse.multinom.logp <- function (n, x, a) {

  # Small number to ensure we never compute log(0).
  e <- 1e-6  

  # Get the maximum-likelihood (ML) estimate of the nonzero
  # multinomial probabilities.
  q <- sparse.multinom.ml(n,x)

  # Return the log-likelihood plus L0-penalty term.
  return(sum(n*log(q + e)) - a*sum(x))
}

# ----------------------------------------------------------------------
# Output a geometric sequence that acts as a cooling schedule for
# simulated annealing. The return value T is a vector of length n such
# that T[1] = 1, and T[n] is 0 or close to 0. The rate at which the
# temperatures approach zero is determined by p, which is analogous to
# the probability of success in a geometric distribution. Smaller
# values of p yield sequences that converge to zero more slowly.
cooling.sched.geom <- function (n, p)
  (1 - p)^(0:(n-1))

# ----------------------------------------------------------------------
# This function simulates a non-homogeneous Markov chain for a
# discrete (binary) state space. This is implemented using the
# Metropolis-Hastings algorithm, following Sec. 3.2 of Andrieu et al,
# 2003, "An Introduction to MCMC for Machine Learning," with birth,
# death and swap moves. States x in which sum(x) = 0 are assumed to
# have zero probability.
#
# The non-optional inputs are: x0, the initial state of the Markov
# chain (x0); f, a function f(x) that returns the log-probability of
# binary configuration x, up to a normalizing constant; and T, a
# sequence of inverse temperatures for each iteration of the Markov
# chain. Note that length(x0) must be greater than 2.
#
# If samples.out = FALSE, the return value is an estimate of the
# maximum-probability configuration. If samples.out = TRUE, the return
# value is a list containing the maximum-probability configuration
# (x.map) and a matrix storing the state at each iteration of the
# Markov chain.
#
# Some additional notes about the temperature schedule: if T =
# rep(1,n), this is equivalent to simulating a homogenous Markov chain
# (i.e. running MCMC) for n iterations; if the temperatures decrease
# slowly to 0, so that the inverse temperatures approach Inf, the
# return value is argmax f(x).
mh.bin <- function (x0, f, T, samples.out = FALSE, verbose = TRUE) {
  
  # Initialize the Markov chain.
  x <- x0

  # Get the number of coordinates (n) and the length of the Markov
  # chain to simulate (numiter).
  n       <- length(x)
  numiter <- length(T)
  
  # Keep track of the best solution we have found so far.
  out <- list(x = x,f = f(x))

  # If requested, store the Monte Carlo samples (that is, the states
  # of the Markov chain).
  if (samples.out)
    X <- matrix(0,n,numiter)

  # Repeat for each temperature in the cooling schedule.
  for (iter in 1:numiter) {
    if (verbose)
      caterase(sprintf("%d",iter))

    # Get the number of nonzero coordinates.
    n1   <- sum(x)
    xnew <- x
    if (n1 == n) {

      # DEATH MOVE (n1 = n)
      # When none of the coordinates are zero, the only possible move
      # we can make is a death move, in which we randomly set one of
      # the nonzero coordinates to zero.
      i       <- sample.set(1:n)
      xnew[i] <- 0
      r       <- n/3
    } else if (n1 == 1) {

      # When only one coordinate is nonzero, the only two possible
      # moves we can make are birth and swap.
      if (runif(1) < 0.5) {
        
        # BIRTH MOVE (n1 = 1)
        # Randomly activate one of the zero coordinates.
        i       <- sample.set(which(x == 0))
        xnew[i] <- 1
        r       <- (n-1)/3
      } else {
        
        # SWAP MOVE (n1 = 1)
        # Switch the state of the only nonzero coordinate with one of
        # the zeroed coordinates.
        i       <- which(x == 1)
        j       <- sample.set(which(x == 0))
        xnew[i] <- 0
        xnew[j] <- 1
        r       <- 1
      }
    } else {

      # If we have ended up here, then there is at least one zero
      # coordinate and at least one nonzero coordinate. Randomly
      # execute a birth, death or swap move with equal probability.
      u <- runif(1)
      if (u < 1/3) {

        # BIRTH MOVE (n1 > 1)
        # Randomly activate one of the zero coordinates.
        i       <- sample.set(which(x == 0))
        xnew[i] <- 1
        if (n1 == n-1)
          r <- 3/n
        else
          r <- (n-n1)/(n1+1)
      } else if (u > 2/3) {

        # DEATH MOVE (n1 < n)
        # Randomly set one of the nonzero coordinates to zero.
        i       <- sample.set(which(x == 1))
        xnew[i] <- 0
        if (n1 == 2)
          r <- 3/(n-1)
        else
          r <- n1/(n-n1+1)
      } else {

        # SWAP MOVE (n1 > 1)
        # Switch the state of nonzero coordinate i with zeroed
        # coordinate j.
        i       <- sample.set(which(x == 1))
        j       <- sample.set(which(x == 0))
        xnew[i] <- 0
        xnew[j] <- 1
        r       <- 1
      }
    }
  
    # Compute the acceptance probability, and move to the new state
    # according to this probability.
    fnew <- f(xnew)
    a    <- min(1,r * exp(T[iter] * (fnew - f(x))))
    if (runif(1) < a) { 
      x <- xnew
      
      # Update the solution to the maximization problem.
      if (fnew > out$f)
        out <- list(x = x,f = fnew)
    }

    # If requested, store the state of the Markov chain.
    if (samples.out)
      X[,iter] <- x
  }
  if (verbose)
    cat("\n")

  # Return the best solution we have found and, if requested, the
  # Monte Carlo samples.
  if (samples.out)
    return(list(x.map = out$x,samples = X))
  else
    return(out$x)
}

# ----------------------------------------------------------------------
# This function simulates a non-homogeneous Markov chain for a
# discrete (binary) state space with posterior distribution given by
# log-posterior density function sparse.multinom.logp(counts,x,a). This
# should produce the same or similar output to mh.bin(x0,f,T) when
#
#   f = function (x) sparse.multinom.logp(counts,x,a).
#
# Input u must be a vector of length 4*n, where n = length(T),
# containing uniform random draws from the unit interval [0,1]. The
# easiest solution is to set u = runif(4*n).
#
# This function calls "mh_sparse_multinom_Call", a function compiled
# from C code, using the .Call interface. For more details on how to
# load the C function into R, see the comments accompanying function
# admixture.labeled.Estep.fast in file admixture.R.
mh.sparse.multinom.fast <-
  function (x0, counts, a, T, u, samples.out = FALSE) {

  # Small number to ensure we never compute log(0).
  e <- 1e-6  

  # Get the number of coordinates (n) and the length of the Markov
  # chain to simulate (numiter).
  n       <- length(x0)
  numiter <- length(T)

  # Initialize the Markov chain (x), and keep track of the best
  # solution that we have found so far (xmap). Note that I call "c" to
  # force R to make a copy of x.
  x    <- as.double(x0)
  xmap <- c(x)

  # If requested, store the the states of the Markov chain.
  if (samples.out)
    samples <- matrix(0,n,numiter)
  else
    samples <- 0

  # Execute the C routine using the .Call interface. The main reason
  # for using .Call interface is that there is less of a constraint on
  # the size of the input matrices. (In particular, the vector of
  # inverse temperatures could be very long.) The only components that
  # change are x, xnew, xmap and samples.
  out <- .Call("mh_sparse_multinom_Call",
               x           = x,                 # State of Markov chain.
               xnew        = rep(0,n),          # Proposed state.
               xmap        = xmap,              # Best solution so far.
               T           = as.double(T),      # Inverse temperatures.
               counts      = as.double(counts), # Expected category counts.
               a           = as.double(a),      # Strength of L0-penalty.
               e           = e,                 # Small number.
               u           = as.double(u),      # Uniform random deviates.
               samples     = samples,           # All Markov chain states.
               samples.out = as.integer(samples.out))

  # Return the best solution we have found so far and, if requested,
  # the Monte Carlo samples.
  if (samples.out)
    return(list(x.map = xmap,samples = samples))
  else
    return(xmap)
}
