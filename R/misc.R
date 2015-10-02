# SUMMARY
# -------
# This file contains some function definitions that didn't fit
# anywhere else. Here is an overview of the functions defined in this
# file:
#
#   distribute(x,k)
#   caterase(s)
#   remove.filename.suffix(x)
#   get.ncols.file(file.name,sep)
#   multiply.bycol(X,a)
#   sigmoid(x)
#   softmax.rows(X)
#   logpexp(x)
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
# Removes the filename extension; e.g. "panel.new.txt" becomes
# "panel.new".
remove.filename.suffix <- function (x) {
  y <- strsplit(x,"[.]")[[1]]
  n <- length(y)
  return(paste(y[-n],collapse = "."))
}

# ----------------------------------------------------------------------
# Get the number of columns in a table file delimited by some
# character specified by "sep".
get.ncols.file <- function (file.name, sep) {
  fc <- file(file.name)
  n  <- length(unlist(strsplit(readLines(fc,1)," ")))
  close(fc)
  return(n)
}

# ----------------------------------------------------------------------
# Multiply each column X[,i] by a[i].
multiply.bycol <- function (X, a)
  X * matrix(a,nrow(X),ncol(X),byrow = TRUE)

# ----------------------------------------------------------------------
# Returns the sigmoid of x. The sigmoid function is also known as the
# logistic link function. It is the inverse of logit(x).
sigmoid <- function (x)
  1/(1 + exp(-x))

# ----------------------------------------------------------------------
# This returns softmax(c(0,X[i,])) for each row i of matrix i.
softmax.rows <- function (X)
  exp(cbind(0,X)) / rowSums(exp(cbind(0,X)))

# ----------------------------------------------------------------------
# Returns log(1 + exp(x)). For large x, logpexp(x) should be
# approximately x. The computation is performed in a numerically
# stable manner.
logpexp <- function (x) {

  # For large entries, log(1 + exp(x)) is effectively the same as x.
  y <- x

  # Find entries of x that are not large. For these entries, compute
  # log(1 + exp(x)).
  i    <- which(x < 16)
  y[i] <- log(1 + exp(x[i]))
  return(y)
}

# ----------------------------------------------------------------------
# Return the parameters of the Beta distribution that yield a random
# variable with mean mu and variance s.
beta.mean.var <- function (mu, s) {
  a <- mu^2*((1-mu)/s - 1/mu)
  return(pmax(0,c(a,a*(1/mu - 1))))
}

