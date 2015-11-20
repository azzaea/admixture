# SUMMARY
# -------
# This file contains functions for reading data files used in the R
# scripts. Here is an overview of the functions defined in this file:
#
#   read.traw.file(file.name)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Read the genotype data from the .traw file. After the header line,
# each line of this space-delimited file corresponds to a single
# marker. The first 6 columns of the file, from left to right, are (1)
# chromosome, (2) marker id, (3) genetic distance in cM, (4) base-pair
# position, (5) counted allele, and (6) other allele(s). After that,
# each column gives the genotype "dosage", or the number of "counted"
# alleles in the genotype.
#
# The return value is a list with two elements: a data frame
# containing the marker map positions ("map"), and the p x n matrix of
# genotype allele counts, where p is the number of markers, and n is
# the number of samples in the panel.
read.traw.file <- function (file.name) {

  # Get the number of markers.
  p <- get.ncols.file(file.name,sep = " ") - 6

  # Read the marker data from the .traw file.
  out <- fread(file.name,sep = " ",header = TRUE,verbose = FALSE,
               showProgress = FALSE,stringsAsFactors = FALSE,
               colClasses = c("integer","character","numeric","integer",
                 "character","character",rep("numeric",p)))

  # Get the marker map positions.
  class(out) <- "data.frame"
  map        <- out[1:6]
  names(map) <- c("chr","id","dist","pos","A1","A2")
  
  # Get the n x p genotype matrix.
  geno          <- out[-(1:6)]
  rownames(out) <- map$id
  rm(out)
  geno <- t(geno)
  geno <- as.matrix(geno)
  
  # Return a list containing (1) the marker map positions, and (2) the
  # genotype matrix.
  return(list(map = map,geno = geno))
}
