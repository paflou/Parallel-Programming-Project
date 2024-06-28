#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <omp.h>

#define NUM_THREADS 8

#define MAXVARS		(250)	/* max # of variables	     */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */

/* prototype of local optimization routine, code available in torczon.c */
extern void mds(double *startpoint, double *endpoint, int n, double *val, double eps, int maxfevals, int maxiter,
         double mu, double theta, double delta, int *ni, int *nf, double *xl, double *xr, int *term);


/* global variables */
unsigned long funevals = 0;

/* Rosenbrock classic parabolic valley ("banana") function */
double f(double *x, int n)
{
    double fv;
    int i;

	/*avoid race conditions*/
	#pragma omp atomic
    funevals++;
	fv = 0.0;
    for (i=0; i<n-1; i++)   /* rosenbrock */
        fv = fv + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
		usleep(1);	/* do not remove, introduces some artificial work */

    return fv;
}


double get_wtime(void)
{
    struct timeval t;

    gettimeofday(&t, NULL);

    return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}

typedef struct calculations {
	/* information about the best point found by multistart */
	double pt[MAXVARS];
	double fx;
	int trial;
	int nt;
	int nf;
} calc;

int main(int argc, char *argv[])
{
	/* problem parameters */
	int nvars = 4;		/* number of variables (problem dimension) */
	int ntrials = 64;	/* number of trials */
	double lower[MAXVARS], upper[MAXVARS];	/* lower and upper bounds */

	omp_set_num_threads(NUM_THREADS);

	/* mds parameters */
	double eps = EPSMIN;
	int maxfevals = 10000;
	int maxiter = 10000;
	double mu = 1.0;
	double theta = 0.25;
	double delta = 0.25;

	double startpt[MAXVARS], endpt[MAXVARS];	/* initial and final point of mds */
	double fx;	/* function value at the final point of mds */
	int nt, nf;	/* number of iterations and function evaluations used by mds */


	/* local variables */
	int trial, i;
	double t0, t1;

	/* initialization of lower and upper bounds of search space */
	for (i = 0; i < MAXVARS; i++) lower[i] = -2.0;	/* lower bound: -2.0 */
	for (i = 0; i < MAXVARS; i++) upper[i] = +2.0;	/* upper bound: +2.0 */

	calc best;
	best.fx = 1e10;

	calc local_best;
	local_best.fx = 1e10;

	t0 = get_wtime();


	unsigned short buffer[3];

	#pragma omp parallel for private(fx, nt, nf, startpt, endpt, buffer) firstprivate(local_best) schedule(static, 1)
	for (trial = 0; trial < ntrials; trial++) {
		buffer[0] = (short)trial;
		buffer[1] = (short)trial + 1;
		buffer[2] = (short)trial + 2;

		/* starting guess for rosenbrock test function, search space in [-2, 2) */
		for (i = 0; i < nvars; i++) {
			startpt[i] = lower[i] + (upper[i]-lower[i])*erand48(buffer);
		}

		int term = -1;
    	mds(startpt, endpt, nvars, &fx, eps, maxfevals, maxiter, mu, theta, delta,
        &nt, &nf, lower, upper, &term);
#if DEBUG
		printf("\n\n\nMDS %d USED %d ITERATIONS AND %d FUNCTION CALLS, AND RETURNED\n", trial, nt, nf);
		for (i = 0; i < nvars; i++)
			printf("x[%3d] = %15.7le \n", i, endpt[i]);

		printf("f(x) = %15.7le\n", fx);
#endif
	
		/* keep the local best solution to minimize global best access*/
		if (fx < local_best.fx) {
			local_best.trial = trial;
			local_best.nt = nt;
			local_best.nf = nf;
			local_best.fx = fx;
			for (i = 0; i < nvars; i++)
				local_best.pt[i] = endpt[i];
		}

		#pragma omp critical
		{
		if(local_best.fx < best.fx) {
			best = local_best;
			//printf("Process %d, best result yet is f(x) = %15.7le\n", omp_get_thread_num(), best.fx);
		}
		}
	}
	t1 = get_wtime();

	printf("\n\nFINAL RESULTS:\n");
	printf("Elapsed time = %.3lf s\n", t1-t0);
	printf("Total number of trials = %d\n", ntrials);
	printf("Total number of function evaluations = %ld\n", funevals);
	printf("Best result at trial %d used %d iterations, %d function calls and returned\n", best.trial, best.nt, best.nf);
	for (i = 0; i < nvars; i++) {
		printf("x[%3d] = %15.7le \n", i, best.pt[i]);
	}
	printf("f(x) = %15.7le\n", best.fx);

	return 0;
}
