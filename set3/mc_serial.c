#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

inline double f(double x)
{
	return sin(cos(x));
}

// WolframAlpha: integral sin(cos(x)) from 0 to 1  =	0.738643
//							0.73864299803689018
//							0.7386429980368901838000902905852160417480209422447648518714116299

int main(int argc, char *argv[])
{
	double a = 0.0;
	double b = 1.0;
	unsigned long n = 24e8;
	long tseed = time(0);

	if (argc == 2) {
		tseed = atol(argv[1]);
	}
	else if (argc == 3) {
		n = atol(argv[1]);
		tseed = atol(argv[2]);
	}

	const double h = (b-a)/n;
	const double ref = 0.73864299803689018;
	double res = 0;
	double t0, t1;
	unsigned long i;

	srand48(tseed);

	t0 = omp_get_wtime();
	for (i = 0; i < n; i++) {
		double xi;
		xi = drand48();
		res += f(xi);
	}
	res *= h;
	t1 = omp_get_wtime();

	printf("Result=%.16f Error=%e Rel.Error=%e Time=%lf seconds\n", res, fabs(res-ref), fabs(res-ref)/ref, t1-t0);
	return 0;
}
