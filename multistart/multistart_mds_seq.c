#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

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


int main(int argc, char *argv[])
{
	/* problem parameters */
	int nvars = 4;		/* number of variables (problem dimension) */
	int ntrials = 64;	/* number of trials */
	double lower[MAXVARS], upper[MAXVARS];	/* lower and upper bounds */

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

	/* information about the best point found by multistart */
	double best_pt[MAXVARS];
	double best_fx = 1e10;
	int best_trial = -1;
	int best_nt = -1;
	int best_nf = -1;

	/* local variables */
	int trial, i;
	double t0, t1;

	/* initialization of lower and upper bounds of search space */
	for (i = 0; i < MAXVARS; i++) lower[i] = -2.0;	/* lower bound: -2.0 */
	for (i = 0; i < MAXVARS; i++) upper[i] = +2.0;	/* upper bound: +2.0 */

	unsigned short buffer[3];
	t0 = get_wtime();
	for (trial = 0; trial < ntrials; trial++) {
		buffer[0] = trial;
		buffer[1] = trial + 1;
		buffer[2] = trial + 2;

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
		printf("current result is f(x) = %15.7le, best result is f(x) = %15.7le\n", fx, best_fx);

		/* keep the best solution */
		if (fx < best_fx) {
			best_trial = trial;
			best_nt = nt;
			best_nf = nf;
			best_fx = fx;
			for (i = 0; i < nvars; i++)
				best_pt[i] = endpt[i];
		}
	}
	t1 = get_wtime();

	printf("\n\nFINAL RESULTS:\n");
	printf("Elapsed time = %.3lf s\n", t1-t0);
	printf("Total number of trials = %d\n", ntrials);
	printf("Total number of function evaluations = %ld\n", funevals);
	printf("Best result at trial %d used %d iterations, %d function calls and returned\n", best_trial, best_nt, best_nf);
	for (i = 0; i < nvars; i++) {
		printf("x[%3d] = %15.7le \n", i, best_pt[i]);
	}
	printf("f(x) = %15.7le\n", best_fx);

	return 0;
}
