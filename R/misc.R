# SUMMARY
# -------
# This file contains some function definitions that didn't fit
# anywhere else. Here is an overview of the functions defined in this
# file:
#
#   distribute(x,k)
#   caterase(s)
#   get.ncols.file(file.name,sep)
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
# Get the number of columns in a table file delimited by some
# character specified by "sep".
get.ncols.file <- function (file.name, sep) {
  fc <- file(file.name)
  n  <- length(unlist(strsplit(readLines(fc,1)," ")))
  close(fc)
  return(n)
}
