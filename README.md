![Admixture estimates in Human Genome Diversity Panel with K = 7
  ancestral populations](hgdp.gif)

## admixture

*A simple EM implementation of the
[ADMIXTURE](http://dx.doi.org/10.1101/gr.094052.109) model in
R, plus extensions.*

ADMIXTURE is widely used to estimate population structure from
genotype data in part because the computation scales well to
whole-genome genotype data. The 1,000,000+ customers who have taken
the AncestryDNA test have all received their [ethnicity
estimate](http://dna.ancestry.com/resource/whitePaper/AncestryDNA-Ethnicity-White-Paper)
from ADMIXTURE.

- Pros and cons of ADMIXTURE software and this R package: facilitates
  development of extensions, including the one we have developed that
  encourages *sparse* admixture estimates (pro); the motivation for
  developing a quasi-Newton optimization algorithm is that the much
  simpler EM implementation often converges more slowly to the
  solution (con);

- Explain why I didn't develop this as an R package.

- We have a procedure for choosing the L0-penalty strength via
cross-validation, but it isn't demonstrated yet in the R scripts. For
some details, see unction **calc.geno.error.R**.

- Tested using R version 3.2.2.

The admixture source code repository is free software: you can
redistribute it under the terms of the **MIT license**. All the files
in this project are part of varbvs. This project is distributed in the
hope that it will be useful, but **without any warranty**; without
even the implied warranty of **merchantability or fitness for a
particular purpose**. See file [LICENSE](LICENSE) for the full text of
the license.

[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)
[AncestryDNA](http://dna.ancestry.com)
San Francisco, California

### Getting started

First, build the the shared object (.so) files using the following
commands:

    R CMD SHLIB mcmc.c
    R CMD SHLIB admixture.c

I've written two scripts that demonstrate usage of the algorithm.

Script **example.admixture.R** uses the EM algorithm to predict
admixture proportions when we have a reference set of labeled,
single-origin samples.

![Admixture estimates in simulated genotype data](example-sim-error.gif)

### The admixture.em function

Estimate population-specific allele frequencies and admixture
proportions in unlabeled samples from genotypes.

#### Usage

    admixture.em(X, K, z = NULL, e = 0.001, a = 0, F = NULL, Q = NULL,
    			 tolerance = 1e-4, max.iter = 1000, exact.q = FALSE,	
                 cg = FALSE,mc.cores = 1, verbose = TRUE, T = 1)

#### Arguments
				 
The non-optional inputs are as follows:

   X   n x p genotype matrix, where n is the number of
       samples and p is the number of markers;

   K   number of ancestral populations;

   z   vector giving the population of origin for each of the samples
       (an integer between 1 and K), or NA if the sample is unlabeled.
       If set to NULL, or not specified, all samples are unlabeled.

The return value is a list with two list elements: F, the p x k
matrix of population-specific allele frequency estimates; and Q, the
n x k matrix of estimated admixture proportions, in which each row
of Q sums to 1. For labeled samples, the admixture proportions are
Q[i,k] = 1 when z[i] = k, otherwise all the other entries are
exactly zero.

There are two variations to the M-step update for the Q matrix. When
the number of ancestral populations is small (k < 20), it is
feasible to compute the L0-penalized estimate exactly by
exhaustively calculating the posterior probability for each possible
choice of the nonzero admixture proportions. Setting exact.q = TRUE
will activate this option. However, for larger k, it is not feasible
to compute the exact solution because the number of ways of choosing
nonzero admixture proportions is too large. Instead, setting exact.q
= FALSE computes an approximate solution using a simulated annealing
algorithm. In this case, it is necessary to set input T. For an
explanation of input T, see function update.q.sparse.approx.

The cg parameter specifies the M-step update for the F matrix. When
cg = FALSE, the binomial success rates are updated using the
standard M-step solution that is derived by finding the roots of the
partial derivatives of the expected complete log-likelihood. When cg
= TRUE, the standard M-step update is adjusted using the conjugate
gradient algorithm (specifically, using the Hestenes-Stiefel update
formula). In some cases, I've found that the conjugate gradient
upgrade can lead to improvements in the convergence rate of the EM
iterates.

#### Value

*Details about output go here.*

### Overview of the R files in this repository

*Details go here.*

### Output from running example.admixture.R

The R console output should look something like this:

    > source("example.admixture.R")
    Simulating data with the following settings:
    markers               200
    ancestral populations 20
    training samples      100
	test samples          300
	Generating genotype data for single-origin training samples.
	Generating genotype data for test samples.
	Setting 1% of the genotypes to NA.
	Computing maximum-likelihood admixture proportion estimates.
	iter delta-F delta-Q -beta-
	431 9.8e-05 6.7e-05 00.484
	Computation took 0.6 min.
	Computing L0-penalized admixture proportion estimates.
	iter delta-F delta-Q -beta-
	 119 9.9e-05 3.0e-05 00.473
	Computation took 0.2 min.
	Overlap between estimated and ground-truth admixture proportions:
	   0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.85 0.9 0.95
	ML 0   0   0   1   1  11  26  74  52   61  37   37
	L0 0   0   1   0   6   1  13  44  18   36  42  139

	Number of contributing ancestral populations (>1%):
	    ML
	true  1  2  3  4  5  6  7  8  9 10 11 12
	   1 19 32 25 13  7  2  2  0  0  0  0  0
	   2  0  4 10 22 34 18  5  6  0  1  0  0
	   4  0  0  0  0  6 15 25 24 15  9  3  3
	   
	    L0
	true  1  2  3  4  5
	   1 96  4  0  0  0
	   2  0 81 19  0  0
	   4  1  2 35 47 15

### Output from running example.sim.R

The R console output should look something like this:

    > source("example.sim.R")
    Loading allele frequency data.
    Generating population allele frequencies.
    Generating data.
    Fitting admixture model to data.
    iter delta-F delta-Q -beta-
    510 4.8e-04 4.1e-04 00.489
    Computation took 1.8 min.
    Fitting L0-penalized admixture model to data.
    iter delta-F delta-Q -beta-
    153 3.0e-04 4.5e-04 -0.112
    Computation took 1.2 min.

    Overlap between estimated and ground-truth admixture proportions:
    0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.85 0.9 0.95
	ML 0   0   0   0   0   0   0  53 170  197  76    4
	L0 0   0   0   0   0   0   0   0   4   38 189  269

	Number of contributing ancestral populations (>1%):
	       1   2  3   4   5   6  7  8 9 10
	true 254 246  0   0   0   0  0  0 0  0
	ML     1   2 38 105 143 118 57 27 9  0
	L0   186 217 84  12   1   0  0  0 0  0
			
