admixture
=========

*Give overview of admixture repository here.*

*Items*:
- Tested using R version 3.2.2.
- Authorship.

### Output from running example.admixture.R

The output to the R console should look something like this:

### Output from running example.sim.R

The output to the R console should look something like this:

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
			
