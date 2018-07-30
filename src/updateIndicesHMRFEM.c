#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* Destructively updates the indices in sZ */
SEXP updateIndicesHMRFEM(SEXP blocks, SEXP sneighbors, SEXP snumnei, SEXP smaxneigh,
						 SEXP beta, SEXP sk, SEXP sZ, SEXP sden)
{
    if (TYPEOF(sneighbors) != INTSXP)
		error("'neighbors' must be of type 'integer'.");
    if (TYPEOF(smaxneigh) != INTSXP)
		error("'nneigh' must be of type 'integer'.");
    if (TYPEOF(snumnei) != INTSXP)
		error("'numnei' must be of type 'integer'.");
    if (TYPEOF(sk) != INTSXP)
		error("'k' must be of type 'integer'.");
    if (TYPEOF(sZ) != INTSXP)
		error("'Z' must be of type 'integer'.");
    //if (TYPEOF(scheck) != REALSXP)
	//	error("'check' must be of type 'double'.");
    if (TYPEOF(sden) != REALSXP)
		error("'den' must be of type 'double'.");
    
    int *neighbors = INTEGER(sneighbors);
    int nneigh = asInteger(smaxneigh);
	int *numnei = INTEGER(snumnei);
	double d_beta = asReal(beta);

    if (nneigh <= 0)
		error("The number of neighbors must be positive.");
    int k = asInteger(sk);
    if (k <= 0)
		error("The number of components must be positive.");
    int *Z = INTEGER(sZ);
    int ldZ = LENGTH(sZ) / k;
    if (ldZ <= 0)
		error("The leading dimension of 'Z' must be positive.");
    
    int ldN = LENGTH(sneighbors) / nneigh; //nneigh is max number of neighbors
    if (ldN <= 0)
		error("The leading dimension of 'neighbors' must be positive.");
    int ldD = LENGTH(sden) / k;
    if (ldD <= 0)
		error("The leading dimension of 'den' must be positive.");
    
    if (ldZ - 1 != ldN || ldZ - 1 != ldD || ldN != ldD)
		error("The leading dimension of 'Z', 'neighbors' and 'den' do not match.");
    
    //int ldC = LENGTH(scheck) / k;
    double *den = REAL(sden);
    //double *check = REAL(scheck);
    double *prob = (double *) R_alloc(k, sizeof(double));
    double *log_prob = (double *) R_alloc(k, sizeof(double));
    int *Ni = (int *) R_alloc(k, sizeof(int));
    int nblocks = LENGTH(blocks);
    int b;
    
    for (b = 0; b < nblocks; b++) {
		SEXP spoints = VECTOR_ELT(blocks, b);
		int n = LENGTH(spoints);
		int *points = INTEGER(spoints);
		int i;
		
		for (i = 0; i < n; i++) {
			int j, m;
			double s = 0.0;
			int number = 0;
			/* compute the posterior weights for the different classes */
			for (j = 0; j < k; j++) {
				Ni[j] = 0;
				//for (m = 0; m < nneigh; m++) {
				for (m = 0; m < numnei[points[i] - 1]; m++) {
					int mm = neighbors[points[i] - 1 + m * ldN] - 1;
					Ni[j] += Z[mm + j * ldZ];
				}
				//Rprintf("point %d j %d val %d\n", points[i], j, Ni[j]);
			}
			
			/*compute the index to use the check table*/
			/*    
			for (j = 0; j < k; j++) {
				number += Ni[j] * pow(nneigh+1, j);
			}
			*/
			
			for (j = 0; j < k; j++) {
				double exp_val = (double) exp(d_beta * Ni[j]);
				int ix = points[i] - 1 + j * ldD;
				prob[j] = exp(log(den[points[i] - 1 + j * ldD]) + log(exp_val));
				log_prob[j] = log(den[points[i] - 1 + j * ldD]) + log(exp_val);
				s += prob[j];
			}
			
			/* generate new Z[points[i] - 1,] entries */
			for (j = 0; j < k; j++)
				Z[points[i] - 1 + j * ldZ] = 0;
            /*find the postion of prob where prob takes the maximum of prob*/
			s = log_prob[0];
			m = 0;
			for (j = 1; j < k; j++) {
				if (log_prob[j] > s){
					m = j;
					s = log_prob[j];
				}
			}

			Z[points[i] - 1 + m * ldZ] = 1;
		}
    }
    
    return sZ;
}

