# SUMMARY
# -------
# This file contains some function definitions that didn't fit
# anywhere else. Here is an overview of the functions defined in this
# file:
#
#   distribute(x,k)
#   caterase(s)
#   beta.mean.var(mu,s)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Distribute the elements of x evenly (or as evenly as possible) into
# k list elements.
distribute <- function (x, k)
  split(x,rep(1:k,length.out = length(x)))

# ----------------------------------------------------------------------
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s)
    cat(s,rep("\b",nchar(s)),sep = "")

# ----------------------------------------------------------------------
# Return the parameters of the Beta distribution that yield a random
# variable with mean mu and variance s.
beta.mean.var <- function (mu, s) {
  a <- mu^2*((1-mu)/s - 1/mu)
  return(pmax(0,c(a,a*(1/mu - 1))))
}

