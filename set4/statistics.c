#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


void compute_max_density(double *rho_, int N)
{
	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (rho_[i*N + j] > max_rho)
			{
				max_rho = rho_[i*N + j];
				max_i = i;
				max_j = j;
			}
		}

	printf("=====================================\n");
	printf("Output of compute_max_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);
}


void compute_max_density_omp(double *rho_, int N)
{
	// TODO 2: parallelize this function with OpenMP. (The provided code is identical to that of compute_max_density()).

	// rho_: matrix of size NxN, allocated as one dimensioal array. rho[i*N+j] corresponds to rho[i][j]
	// This routine finds the value of max density (max_rho) and its location (max_i, max_j) - it assumes there are no duplicate values
	double max_rho;
	int max_i,  max_j;

	max_rho = rho_[0];
	max_i = 0;
	max_j = 0;

	#pragma omp parallel for
	for(int k=0;k<N*N;k++)
	{
		int i = k/N;
		int j = k%N;

		if (rho_[i*N + j] > max_rho){
		#pragma omp critical
		{
			max_rho = rho_[i*N + j];
			max_i = i;
			max_j = j;
		}
		}
	}

	printf("=====================================\n");
	printf("Output of compute_max_omp_density():\n");
	printf("Max rho: %.16f\n", max_rho);
	printf("Matrix location: %d %d\n", max_i, max_j);
}



int main(int argc, char *argv[])
{
	int N = 4096;

	double *rho = (double *)malloc(N*N*sizeof(double));	// NxN 'density' matrix


	// matrix initialization
	srand48(1);
	for (int i = 0 ; i < N; i++) 
	for (int j = 0 ; j < N; j++) 
		rho[i*N+j] = drand48();			// these will be density values when we study diffusion


	compute_max_density(rho, N);			// sequential version 

	compute_max_density_omp(rho, N);		// your parallel version: must produce the same results
	

	free(rho);

	return 0;
}

