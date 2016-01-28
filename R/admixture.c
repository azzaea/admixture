// These include files have a bunch of definitions to interface C
// routines in R.
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
SEXP genoprob_given_q_Call (SEXP fp, SEXP qp, SEXP ep, SEXP pxp) {

  // Get the inputs.
  double* f  = REAL(fp);   // Allele frequencies.
  double* q  = REAL(qp);   // Admixture proportions.
  double* px = REAL(pxp);  // Genotype probabilities.
  double  e  = *REAL(ep);  // Genotype error probability.

  // These variables keep track of loop iterations.
  R_xlen_t i, j;
  int      x;

  // Get the number of ancestral populations.
  R_xlen_t n = length(qp);

  // Intermediate quantities used to calculate the genotype probabilities.
  double u[2];
  u[0] = e;
  u[1] = 1 - 2*e;

  // Repeat for each combination of ancestral populations.
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (x = 0; x <= 2; x++)

	// Compute the probability that the allele count of the unphased
        // genotype is equal to x given the assignment (i,j) to the
        // population-of-origin indicators.
	px[x] += q[i]*q[j] *
          (u[x == 0] * (1-f[i]) * (1-f[j]) +
           u[x == 1] * (1-f[i]) * f[j] +
           u[x == 1] * f[i]     * (1-f[j]) +
           u[x == 2] * f[i]     * f[j]);

  return R_NilValue;
}

// -----------------------------------------------------------------
// This function is used by admixture.labeled.Estep.fast to compute
// the expected allele counts "n0" and "n1", which are used to update
// the population-level allele frequencies. This function is called in
// R using the .Call interface.
//
// Note that all the entries of input matrices n0 and n1 should be
// initialized to 0, and the population labels z should start at 0,
// not 1.
SEXP admixture_labeled_Estep_Call (SEXP Xp, SEXP Fp, SEXP zp, SEXP ep,
				   SEXP loglp, SEXP n0p, SEXP n1p) {

  // Get the inputs.
  double* X    = REAL(Xp);    // Genotype matrix.
  double* F    = REAL(Fp);    // Allele frequency matrix.
  double* z    = REAL(zp);    // Population indicators.
  double  e    = *REAL(ep);   // Genotype error probability.
  double* logl = REAL(loglp); // Log-likelihood.
  double* n0   = REAL(n0p);   // Expected counts of "0" allele.
  double* n1   = REAL(n1p);   // Expected counts of "1" allele.

  // Get the number of samples (n) and the number of markers (p).
  SEXP     Xdim = getAttrib(Xp,R_DimSymbol);
  R_xlen_t n    = INTEGER(Xdim)[0];
  R_xlen_t p    = INTEGER(Xdim)[1];

  // Posterior probabilities of all possible hidden (phased) genotype
  // configurations at a single locus.
  double r00, r01, r10, r11;

  // Some more variables for storing intermediate values used in the
  // calculations below.
  double   x;      // A genotype.
  double   f0, f1; // An allele frequency.
  double   w;      // Normalizing constant for posterior probabilities.
  R_xlen_t k;      // Population of origin.

  // Initialize the calculation of the log-likelihood.
  *logl = 0;

  // Intermediate quantities used to calculate the genotype likelihoods.
  double u[2];
  u[0] = e;
  u[1] = 1 - 2*e;

  // Repeat for each sample, and for each marker.
  for (R_xlen_t i = 0; i < n; i++)
    for (R_xlen_t j = 0; j < p; j++) {

      // Get the ancestral population of origin for sample i (k), the
      // genotype for sample i at marker j (x), and get the allele
      // frequency of marker j in population k (f1).
      k  = (R_xlen_t) z[i];
      x  = X[i + n*j];
      f1 = F[j + p*k];

      // Compute the posterior probabilities (up to a normalizing
      // constant) for all possible hidden (phased) genotype
      // configurations at each locus. Note that genotype
      // configurations (0,1) and (1,0) always have the same posterior
      // probability.
      f0  = 1 - f1;
      r00 = f0 * f0;
      r01 = f0 * f1;
      r10 = f1 * f0;
      r11 = f1 * f1;
      if (ISNA(x)) {
	r00 *= 0.3333;
	r01 *= 0.3333;
	r10 *= 0.3333;
	r11 *= 0.3333;
      } else {
	r00 *= u[x == 0];
	r01 *= u[x == 1];
	r10 *= u[x == 1];
	r11 *= u[x == 2];
      }
      
      // Compute the normalizing constant for the posterior
      // probabilities.
      w = r00 + r01 + r10 + r11;

      // Add the normalizing constant to the log-likelihood.
      *logl += log(w);

      // Add the posterior probabilities to the expected allele counts.
      n0[j + p*k] += 2*(r00 + r01)/w;
      n1[j + p*k] += 2*(r11 + r10)/w;
    }

  return R_NilValue;
}

// -----------------------------------------------------------------
// This function is used by admixture.Estep.fast to compute the
// expected allele counts "n0" and "n1" and the expected population
// counts "m", which are used to update the population-level allele
// frequencies and the admixture proportions. This function is called
// in R using the .Call interface.
//
// Note that all the entries of input matrix m should be initialized
// to 0.
SEXP admixture_unlabeled_Estep_Call (SEXP Xp, SEXP Fp, SEXP Qp, SEXP ep,
				     SEXP mp, SEXP n0p, SEXP n1p, SEXP loglp,
				     SEXP r00p, SEXP r01p, SEXP r10p, 
				     SEXP r11p) {

  // Get the inputs.
  double* X    = REAL(Xp);    // Genotype matrix.
  double* F    = REAL(Fp);    // Allele frequency matrix.
  double* Q    = REAL(Qp);    // Admixture proportions matrix.
  double  e    = *REAL(ep);   // Genotype error probability.
  double* m    = REAL(mp);    // Expected population counts.
  double* n0   = REAL(n0p);   // Expected counts of "0" allele.
  double* n1   = REAL(n1p);   // Expected counts of "1" allele.
  double* logl = REAL(loglp); // Log-likelihood.
  double* r00  = REAL(r00p);  // Posterior probability of (0,0) genotype.
  double* r01  = REAL(r01p);  // Posterior probability of (0,1) genotype.
  double* r10  = REAL(r10p);  // Posterior probability of (1,0) genotype.
  double* r11  = REAL(r11p);  // Posterior probability of (1,1) genotype.

  // Get the number of samples (n), the number of markers (p), and the
  // number of ancestral populations (K).
  SEXP     Xdim = getAttrib(Xp,R_DimSymbol);
  SEXP     Qdim = getAttrib(Qp,R_DimSymbol);
  R_xlen_t n    = INTEGER(Xdim)[0];
  R_xlen_t p    = INTEGER(Xdim)[1];
  R_xlen_t K    = INTEGER(Qdim)[1];

  // Some more variables to store intermediate values used in the
  // calculations below.
  double x;           // A genotype.
  double f1, f2, f12; // Allele frequencies.
  double q1, q2, q12; // Admixture proportions.  
  double w;           // Normalizing constant for posterior probabilities.
  double r1, r2;

  // These variables keep track of the loop iterations or locations in
  // an array.
  R_xlen_t i, j, k, k1, k2, d1, d2, t;
  R_xlen_t t1, t2;

  // This stores intermediate quantities used to calculate the genotype
  // likelihoods.
  double u[2];

  // These two variables store the number of ancestral populations
  // that contribute to an individual's genome (nd), and the indices
  // of the contributing ancestral populations (demes).
  int      nd;
  R_xlen_t demes[K];

  // Initialize the calculation of the log-likelihood.
  *logl = 0;

  // Repeat for each sample.
  for (i = 0; i < n; i++) {

    // Get indices of the ancestral populations contributing to the
    // genome of individual i.
    for (k = 0, nd = 0; k < K; k++) 
      if (Q[i + n*k] > 0) {
	demes[nd] = k;
	nd++;
      }

    // Repeat for each marker.
    for (j = 0; j < p; j++) {

      // The inner and outer foor loops iterate over each pair of
      // contributing ancestral populations.
      w = 0;
      for (k1 = 0; k1 < nd; k1++) {

	// Get the genotype for sample i at marker j, the allele
	// frequency of marker j in population k1, and the admixture
	// proportion for population k1 in sample i.
	x  = X[i + n*j];
	d1 = demes[k1];
	f1 = F[j + p*d1];
	q1 = Q[i + n*d1];

	// Calculate intermediate quantities.
	if (ISNA(x)) {
	  u[0] = 0.3333;
	  u[1] = 0.3333;
	} else {
	  u[0] = e;
	  u[1] = 1 - 2*e;
	}

	// Compute the posterior probabilities (up to a normalizing
	// constant) for all possible hidden (phased) genotype
	// configurations at each locus, given the assigment (k1,k1)
	// to the population-of-origin indicators.
	t      = k1*(1 + nd);
	f12    = f1 * f1;
	q12    = q1 * q1;
	r00[t] = q12 * u[x == 0] * (1 - 2*f1 + f12);
	r01[t] = q12 * u[x == 1] * (f1 - f12);
	r10[t] = q12 * u[x == 1] * (f1 - f12);
	r11[t] = q12 * u[x == 2] * f12;

	// Add these posterior probabilities to the normalizing
	// constant.
	w += r00[t] + r01[t] + r10[t] + r11[t];

	for (k2 = k1 + 1, t = k1*(1+nd) + nd, t2 = k1*(1+nd) + 1; 
	     k2 < nd; k2++, t += nd, t2++) {

	  // Get the allele frequency of marker j in population k2,
	  // and the admixture proportion for population k2 in sample i.
	  d2 = demes[k2];
	  f2 = F[j + p*d2];
	  q2 = Q[i + n*d2];

	  // Compute the posterior probabilities (up to a normalizing
	  // constant) for all possible hidden (phased) genotype
	  // configurations at each locus, given the assigment (k1,k2)
	  // to the population-of-origin indicators. Note that
	  // genotype configurations (0,1) and (1,0) do *not* have the
	  // same posterior probability here (except when k1 and k2
	  // happen to be the same).
	  f12    = f1 * f2;
	  q12    = q1 * q2;
	  r00[t] = q12 * u[x == 0] * (1 - f1 - f2 + f12);
	  r01[t] = q12 * u[x == 1] * (f2 - f12);
	  r10[t] = q12 * u[x == 1] * (f1 - f12);
	  r11[t] = q12 * u[x == 2] * f12;

	  // Add these posterior probabilities to the normalizing
	  // constant.
	  r00[t2] = r00[t];
	  r01[t2] = r10[t];
	  r10[t2] = r01[t];
	  r11[t2] = r11[t];
	  w += 2*(r00[t] + r01[t] + r10[t] + r11[t]);
	}
      }

      // Add the normalizing constant to the log-likelihood.
      *logl += log(w);

      // Add the (normalized) posterior probabilities to the expected
      // population counts and expected allele counts. The inner and
      // outer for loops iterate over all pairs of contributing
      // ancestral populations.
      for (k1 = 0; k1 < nd; k1++) {
	d1 = demes[k1];
	t  = k1*(nd + 1);
	r1 = 2*(r00[t] + r01[t])/w;
	r2 = 2*(r10[t] + r11[t])/w;
	m[i  + n*d1] += r1 + r2;
	n0[j + p*d1] += r1;
	n1[j + p*d1] += r2;
	
      	for (k2 = 0; k2 < k1; k2++) {
	  t1            = k1 + nd*k2;
	  n0[j + p*d1] += 2*(r00[t1] + r01[t1])/w;
	  n1[j + p*d1] += 2*(r10[t1] + r11[t1])/w;
	}

      	for (k2 = k1 + 1; k2 < nd; k2++) {
	  t1 = k1 + nd*k2;
	  d2 = demes[k2];
	  r1 = 2*(r00[t1] + r01[t1])/w;
	  r2 = 2*(r10[t1] + r11[t1])/w;
	  m[i  + n*d1] += r1 + r2;
	  m[i  + n*d2] += r1 + r2;
	  n0[j + p*d1] += r1;
	  n1[j + p*d1] += r2;
	}
      }
    }
  }

  return R_NilValue;
}
