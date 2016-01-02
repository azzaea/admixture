# SUMMARY
# -------
# This file contains some function definitions that didn't fit
# anywhere else. Here is an overview of the functions defined in this
# file:
#
#   sigmoid(x)
#   logit(x)
#   softmax(x)
#   softmax.inverse(y)
#   softmax.rows(x)
#   distribute(x,k)
#   caterase(s)
#   get.ncols.file(file.name,sep)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns the sigmoid of x. The sigmoid function is also known as the
# logistic link function. It is the inverse of logit(x).
sigmoid <- function (x)
  1/(1 + exp(-x))

# ----------------------------------------------------------------------
# The logit function, which is the inverse of the sigmoid function. I
# add a small constant to the numerator and denominator to ensure that
# the return value is always a finite number.
logit <- function (x)
  log((x + 1e-8)/(1 - x + 1e-8))

# ----------------------------------------------------------------------
# Return the softmax of vector c(0,x). 
softmax <- function (x) {
  x <- c(0,x)
  c <- max(x)
  x <- exp(x - c)
  return(x/sum(x))
}

# ----------------------------------------------------------------------
# This is the inverse of the softmax function on c(0,x).
softmax.inverse <- function (y) {
  z <- -log(y[1] + 1e-8)
  y <- y[-1]
  return(z + log(y + 1e-8))
}

# ----------------------------------------------------------------------
# Given n x k matrix X, return an matrix n x (k+1) matrix Y such that
# Y[i,] is the softmax of vector c(0,X[i,]).
softmax.rows <- function (X)
  t(apply(X,1,softmax))

# ----------------------------------------------------------------------
# Given n x (k+1) matrix Y, return an matrix n x k matrix X such that
# Y[i,] is the softmax of vector c(0,X[i,]).
softmax.inverse.rows <- function (Y)
  t(apply(Y,1,softmax.inverse))

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
# Get the number of columns in a table file delimited by some
# character specified by "sep".
get.ncols.file <- function (file.name, sep) {
  fc <- file(file.name)
  n  <- length(unlist(strsplit(readLines(fc,1)," ")))
  close(fc)
  return(n)
}
