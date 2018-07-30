#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* Destructively updates the indices in sZ */
/* openMP is used*/
SEXP updateIndicesX(SEXP blocks, SEXP sneighbors, SEXP snneigh,
				   SEXP sk, SEXP sZ, SEXP scheck, SEXP sden)
{
	if (TYPEOF(sneighbors) != INTSXP)
		error("'neighbors' must be of type 'integer'.");
	if (TYPEOF(snneigh) != INTSXP)
		error("'nneigh' must be of type 'integer'.");
	if (TYPEOF(sk) != INTSXP)
		error("'k' must be of type 'integer'.");
    if (TYPEOF(sZ) != INTSXP)
		error("'Z' must be of type 'integer'.");
	if (TYPEOF(scheck) != REALSXP)
		error("'check' must be of type 'double'.");
    if (TYPEOF(sden) != REALSXP)
		error("'den' must be of type 'double'.");
    
	
	int *neighbors = INTEGER(sneighbors);
    int nneigh = asInteger(snneigh);
    if (nneigh <= 0)
		error("The number of neighbors must be positive.");
	int k = asInteger(sk);
	if (k <= 0)
		error("The number of components must be positive.");
	int *Z = INTEGER(sZ);
    int ldZ = LENGTH(sZ) / k;
    if (ldZ <= 0)
		error("The leading dimension of 'Z' must be positive.");
	int ldN = LENGTH(sneighbors) / nneigh;
	if (ldN <= 0)
		error("The leading dimension of 'neighbors' must be positive.");
	int ldD = LENGTH(sden) / k;
	if (ldD <= 0)
		error("The leading dimension of 'den' must be positive.");
	
	if (ldZ - 1 != ldN || ldZ - 1 != ldD || ldN != ldD)
		error("The leading dimension of 'Z', 'neighbors' and 'den' do not match.");
	
    int ldC = LENGTH(scheck) / k;
    double *den = REAL(sden);
    double *check = REAL(scheck);
    int nblocks = LENGTH(blocks);
    int b;
	
	
    GetRNGstate();
	
	
    for (b = 0; b < nblocks; b++) {
		SEXP spoints = VECTOR_ELT(blocks, b);
		int n = LENGTH(spoints);
		int *points = INTEGER(spoints);
		double *U = (double *) R_alloc(n, sizeof(double));
		/*double *prob = (double *) R_alloc(k, sizeof(double));*/
		/*  int *Ni = (int *) R_alloc(k, sizeof(int));*/
   
		int i;
             
		for (i = 0; i < n; i++)
			U[i] = unif_rand();
		

#pragma omp parallel for firstprivate(k, ldD, ldN, ldZ, ldC, nneigh, neighbors, \
				      Z, den, check, points, U)
		for (i = 0; i < n; i++) {
			int j, m;
			double s = 0.0;
			int number = 0;
			double prob[10];
			int Ni[10];
			/* compute the posterior weights for the different classes */
			for (j = 0; j < k; j++) {
				Ni[j] = 0;
				for (m = 0; m < nneigh; m++) {
					int mm = neighbors[points[i] - 1 + m * ldN] - 1;
					Ni[j] += Z[mm + j * ldZ];
				}
			}
			
            /*compute the index to use the check table*/            
			for (j = 0; j < k; j++) {
				number += Ni[j] * pow(nneigh+1, j);
			}
			
			
			for (j = 0; j < k; j++) {
				/* use tabled value check[or + j*ldC] = exp(beta * Ni[j]) */
				prob[j] = den[points[i] - 1 + j * ldD] * check[number + j * ldC];
				s += prob[j];
			}
			
			/* fix up the weights and convert to probabilities */
			if (s != 0.0 && R_FINITE(s))
				for (j = 0; j < k; j++)
					prob[j] = prob[j] / s;
			else
				for (j = 0; j < k; j++)
					prob[j] = 1.0 / k;
			
			/* generate new Z[points[i] - 1,] entries */
			for (j = 0; j < k; j++)
				Z[points[i] - 1 + j * ldZ] = 0;
			for (m = 0, s = 0.0, j = 0; j < k - 1; j++) {
				s += prob[j];
				if (U[i] > s) m++;
			}
			Z[points[i] - 1 + m * ldZ] = 1;
		}
    }
	
    PutRNGstate();
	
    return sZ;
}

/* Non-destructive version */
SEXP updateIndices(SEXP blocks, SEXP sneighbors, SEXP snneigh,
				   SEXP sk, SEXP sZ, SEXP scheck, SEXP sden)
{
	SEXP val;
    PROTECT(sZ = duplicate(sZ));
    val = updateIndicesX(blocks, sneighbors, snneigh,
						 sk, sZ, scheck, sden);
    UNPROTECT(1);
    return val;
}
