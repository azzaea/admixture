![Admixture estimates in Human Genome Diversity Panel with K = 7
  ancestral populations](hgdp.gif)

# admixture

*A simple EM implementation of the
[ADMIXTURE](http://dx.doi.org/10.1101/gr.094052.109) model in
R, plus extensions.*

The ADMIXTURE
[software](https://www.genetics.ucla.edu/software/admixture/) is
widely used to estimate population structure from genotype data in
part because the computation scales well to whole-genome genotype
data. For example, the million+ people who have taken the AncestryDNA
test have all received their [ethnicity
estimate](http://dna.ancestry.com/resource/whitePaper/AncestryDNA-Ethnicity-White-Paper)
from ADMIXTURE.

I have developed a simple, alternative implementation of ADMIXTURE
that computes maximum-likelihood estimates of the admixture
proportions and population allele frequencies using the expectation
maximization (EM) algorithm. (See
[admixture.barebones.R](R/admixture.barebones.R) and
[admixture.barebones.demo.R](R/admixture.barebones.demo.R) for a
extremely simple, or "bare bones", implementation that actually works,
albeit slowly!)

The ADMIXTURE software is implemented using a quasi-Newton method, and
will typically converge much more quickly to a solution than the EM
algorithm. I've modified the model to allow for genotype errors, and
this seems to help convergence to some extent. In any case, the hope
is that this very simple implementation will facilitate development of
extensions to ADMIXTURE. One extension I have developed here is a
modification to the optimization (M-step) that encourages *sparse*
admixture estimates.

This code was tested using R version 3.2.2.

The admixture source code repository is free software: you can
redistribute it under the terms of the **MIT license**. All the files
in this project are part of **admixture**. This project is distributed
in the hope that it will be useful, but **without any warranty**;
without even the implied warranty of **merchantability or fitness for
a particular purpose**. See file [LICENSE](LICENSE) for the full text
of the license.

[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
[AncestryDNA](http://dna.ancestry.com)<br>
San Francisco, California

### Getting started

Before using any of the functions in R, build the the shared object
(.so) files with the following commands:

    R CMD SHLIB mcmc.c
    R CMD SHLIB admixture.c

I have written three scripts to demonstrate usage of the EM algorithm.

Script [example.admixture.R](R/example.admixture.R) uses the EM
algorithm (function **admixture.em**) to predict admixture proportions
when we have a reference set of labeled, single-origin individuals.

Script [example.sim.R](R/example.sim.R) evaluates accuracy of the
admixture estimates in simulated genotype data, with and without the
L0-penalty term that encourages sparsity in the admixture
proportions. In this example, all the samples are unlabeled.

Script [predict.admix.hgdp.R](R/predict.admix.hgdp.R) uses the EM
function to estimate admixture proportions in the Human Genome
Diversity Panel (HGDP). For instructions on obtaining the HGDP
genotype data, see the comments at the top of this R script.

![Admixture estimates in simulated genotype data](example-sim-error.gif)

### The admixture.em function

Estimate admixture proportions in unlabeled samples from their
genotypes.

#### Usage

    admixture.em(X, K, z = NULL, e = 0.001, a = 0, F = NULL, Q = NULL,
    			 tolerance = 1e-4, max.iter = 1000, exact.q = FALSE,	
                 cg = FALSE,mc.cores = 1, verbose = TRUE, T = 1)

#### Arguments
				 
Input **X** is an n x p genotype matrix, where n is the number of
samples and p is the number of biallic genetic markers. Genotypes are
represented as allele counts, so all entries must be 0, 1 or
2. Missing values (NA) are also allowed.

Input **K** is a model parameter specifying the number of ancestral
populations.

Input **z** is a vector giving the population of origin (an integer
between 1 and K) for each of the samples, or NA is the sample is
unlabeled. If z is set to NULL, or is not specified, all samples are
treated as unlabeled.

Input **e** specifies the probably of a genotype error. It must be a
positive number. It can be small (e.g., 1e-6), but note that small
values tend to increase convergence time of the EM algorithm.

Input **a** specifies the strength of the L0-penalty term that
encourages sparsity in the admixture estimates. By default, a = 0,
which means that the L0-penalty term has no effect, and the
maximum-likelihood estimate is returned. I have implemented a
procedure for choosing the L0-penalty strength using cross-validation,
but this procedure isn't demonstrated yet in this code. For some
details on this cross-validation, see function **calc.geno.error**.

Inputs **F** and **Q** are the initial estimates of the population
allele frequencies and admixture proportions, respectively. If these
inputs aren't specified, these model parameters are randomly
initialized. For more details on F and Q, see below.

Input argument **tolerance** specifies the convergence tolerance of
the EM iterates. Convergence is reached when the maximum absolute
difference between the parameters at two successive iterations is less
than the specified tolerance. Input **max.iter** specifies the maximum
number of EM iterations.

There are two variations to the M-step update for Q. When the number
of ancestral populations is small (e.g., K < 20), it is feasible to
compute the optimal solution exactly by exhaustively calculating the
posterior probability for each possible choice of nonzero admixture
proportions.  Setting **exact.q = TRUE** activates this
option. However, for larger K, it is not feasible to compute the exact
solution because the number of ways of choosing nonzero admixture
proportions is too large. Instead, setting **exact.q = FALSE**
computes an approximate solution using simulated annealing. In this
case, it is necessary to set input **T**. For an explanation of T, see
function **update.q.sparse.approx**.

The **cg** parameter controls the M-step update for the F matrix. When
cg = FALSE, the allele frequencies are updated using the standard
M-step solution. When cg = TRUE, the standard M-step update is
adjusted using conjugate gradient with the Hestenes-Stiefel update
formula. I've found that this sometimes improves convergence of the EM
iterates.

Finally, **mc.cores** is the input to mclapply specifying the number
of cores to use in the multicore versions of the E and M-steps.

#### Value

The return value is a list with two elements:

**F**, the p x K matrix of population-specific allele frequency
estimates;

**Q**, the n x K matrix of admixture proportions, where n is the
number of individuals (samples). For labeled samples, the admixture
proportions are Q[i,k] = 1 for z[i] = k.
