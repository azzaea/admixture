// These include files have a bunch of definitions to interface C
// routines in R.
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Returns the L0-penalized log-likelihood given the category counts
// ("counts") and a binary indicator for the nonzero coordinates (x).
// Input a is the strength of the L0-penalty, n is the number of
// categories, and e is a small number included to ensure that log(0)
// is never computed.
double calclogp (const double* counts, const double* x, int n, 
		 double a, double e);

// Returns the sum of the first n elements in vector x.
double sumvector (const double* x, int n);

// Copy the first n elements of vector x into vector y.
void copyvector (const double* x, double* y, int n);

// Choose uniformly at random among the n elements of x that are equal
// to a. Input u is number drawn uniformly at random from [0,1].
int sampleset (const double* x, int n, double a, double u);

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// This function is used by mh.sparse.multinom.fast to implement the
// Metropolis-Hastings algorithm. This function is called in R using
// the .Call interface.
//
// Here I assume that number of coordinates is not large (does not
// need to be stored as an R_xlen_t).
SEXP mh_sparse_multinom_Call (SEXP xp, SEXP xnewp, SEXP xmap_p, SEXP Tp, 
			      SEXP countsp, SEXP ap, SEXP ep, SEXP up,
			      SEXP samples_p, SEXP s_out_p) {

  // Get the inputs.
  double* x       = REAL(xp);          // State of Markov chain.
  double* xnew    = REAL(xnewp);       // Proposed state.
  double* xmap    = REAL(xmap_p);      // Best solution so far.
  double* T       = REAL(Tp);          // Inverse temperatures.
  double* counts  = REAL(countsp);     // Expected category counts.
  double  a       = *REAL(ap);         // Strength of L0-penalty.
  double  e       = *REAL(ep);         // Small number.
  double* u       = REAL(up);          // Uniform random deviates.
  double* samples = REAL(samples_p);   // All Markov chain states.
  int     s_out   = *INTEGER(s_out_p); // If 1, output all states.

  // Get the number of coordinates (n) and the length of the Markov
  // chain to simulate (numiter).
  int      n       = (int) xlength(xp);
  R_xlen_t numiter = xlength(Tp);

  // Some more variables for storing intermediate values used in the
  // calculations below.
  int    n1;            // Number of nonzero coordinates.
  int    i, j;          // Coordinate indices.
  double f, fnew, fmap; // Log-probability up to normalizing constant.
  double r;             // Ratio of proposal probabilities.
  double A;             // Acceptance probability.

  // Compute the log-probability at the current solution (xmap).
  fmap = calclogp(counts,xmap,n,a,e);

  // Repeat for each temperature in the cooling schedule.
  for (R_xlen_t iter = 0; iter < numiter; iter++, u += 4) {

    // Get the number of nonzero coordinates.
    n1 = (int) sumvector(x,n);

    // Propose a new state (xnew).
    copyvector(x,xnew,n);
    if (n1 == n) {

      // DEATH MOVE (n1 = n)
      // When none of the coordinates are zero, the only possible move
      // we can make is a death, in which we randomly set one of the
      // nonzero coordinates to zero.
      i       = sampleset(x,n,1,u[1]);
      xnew[i] = 0;
      r       = (double) n/3;
    } else if (n1 == 1) {

      // When only one coordinate is nonzero, the only two possible
      // moves we can make are a birth or a swap.
      if (u[0] < 0.5) {
        
        // BIRTH MOVE (n1 = 1)
        // Randomly activate one of the zero coordinates.
        i       = sampleset(x,n,0,u[1]);
        xnew[i] = 1;
        r       = (double) (n-1)/3;
      } else {
        
        // SWAP MOVE (n1 = 1)
        // Switch the state of the only nonzero coordinate (i) with
        // one of the zeroed coordinates (j).
        i       = sampleset(x,n,1,u[1]);
        j       = sampleset(x,n,0,u[2]);
        xnew[i] = 0;
        xnew[j] = 1;
        r       = 1;
      }
    } else {

      // If we have ended up here, then there is at least one zero
      // coordinate and at least one nonzero coordinate. Randomly
      // execute a birth, death or swap move with equal probability.
      if (u[0] < (double) 1/3) {

        // BIRTH MOVE (n1 > 1)
        // Randomly activate one of the zero coordinates.
        i       = sampleset(x,n,0,u[1]);
        xnew[i] = 1;
        if (n1 == n-1)
          r = (double) 3/n;
        else
          r = (double) (n-n1)/(n1+1);
      } else if (u[0] > (double) 2/3) {

        // DEATH MOVE (n1 < n)
        // Randomly set one of the nonzero coordinates to zero.
        i       = sampleset(x,n,1,u[1]);
        xnew[i] = 0;
        if (n1 == 2)
          r = (double) 3/(n-1);
        else
          r = (double) n1/(n-n1+1);
      } else {

        // SWAP MOVE (n1 > 1)
        // Switch the state of nonzero coordinate i with zeroed
        // coordinate j.
        i       = sampleset(x,n,1,u[1]);
        j       = sampleset(x,n,0,u[2]);
        xnew[i] = 0;
	xnew[j] = 1;
	r       = 1;
      }
    }

    // Compute the acceptance probability, and move to the new state
    // according to this probability. Note that the acceptance
    // probability may be greater than 1, but this should make no
    // difference.
    f    = calclogp(counts,x,n,a,e);
    fnew = calclogp(counts,xnew,n,a,e);
    A    = r * exp(T[iter] * (fnew - f));
    if (u[3] < A) {

      // Move to the proposed state.
      copyvector(xnew,x,n);

      // Update the solution to the maximization problem.
      if (fnew > fmap) {
	copyvector(xnew,xmap,n);
	fmap = fnew;
      }
    }

    // If requested, store the state of the Markov chain.
    if (s_out)
      copyvector(x,samples + n*iter,n);
  }

  return R_NilValue;
}

// -----------------------------------------------------------------
// Returns the L0-penalized log-likelihood given the category counts
// and a binary indicator for the nonzero coordinates.
double calclogp (const double* counts, const double* x, int n, 
		 double a, double e) {
  double y = 0;  // Log-likelihood term.
  double z = 0;  // Normalizing constant.

  // Get the number of nonzero coordinates.
  double n1 = sumvector(x,n);

  // Compute the sum of the counts corresponding to the nonzero
  // coordinates.
  for (int i = 0; i < n; i++)
    z += x[i] * counts[i];

  if (z == 0) {

    // Compute log p(counts | x) for special case when all the counts
    // corresponding to the nonzero coordinates are zero.
    for (int i = 0; i < n; i++)
      y += counts[i] * log(x[i]/n1 + e);
  } else

    // Compute log p(counts | x) when at least one of the counts
    // corresponding to the nonzero coordinates is not zero.
    for (int i = 0; i < n; i++)
      y += counts[i] * log(x[i]*counts[i]/z + e);
  
  // Return the log-likelihood plus L0-penalty term.
  return(y - a*n1);
}

// -----------------------------------------------------------------
// Copy the first n elements of vector x into vector y.
void copyvector (const double* x, double* y, int n) {
  for (int i = 0; i < n; i++)
    y[i] = x[i];
}

// -----------------------------------------------------------------
// Returns the sum of the first n elements in vector x.
double sumvector (const double* x, int n) {
  double y = 0;
  for (int i = 0; i < n; i++)
    y += x[i];
  return y;
}

// -----------------------------------------------------------------
// Choose uniformly at random among the n entries of vector x that are
// equal to a.
int sampleset (const double* x, int n, double a, double u) {
  double r = 0; 

  // Get the number of entries equal to a.
  int na = 0;
  for (int i = 0; i < n; i++)
    na += (x[i] == a);

  // Choose one of the entries in x equal to a.
  for (int i = 0; i < n; i++) {
    r += (x[i] == a);
    if (na*u < r)
      return i;
  }

  return n-1;
}
