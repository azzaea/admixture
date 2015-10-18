Put an example here of admixture in HGDP data.

admixture
=========

A simple EM implementation of the
[ADMIXTURE](http://dx.doi.org/10.1101/gr.094052.109) model in
R, plus a couple of extensions.

*Give overview of admixture repository here.*

*Items*:

- Show figure from running this on the simulated data.

- Show figure from running this on the HGDP data.

- Pros and cons of ADMIXTURE software and this R package.

- We have a procedure for choosing the L0-penalty strength via
  cross-validation, but it isn't demonstrated yet in the R scripts.

- Tested using R version 3.2.2.

- Authorship.

![Alt text](example-sim-error.gif)

### Getting started

*Text goes here.*

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
			
