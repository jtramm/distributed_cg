#include "cg_header.h"

int main(int argc, const char** argv)
{
	double start, stop;

	// n x n square grid
	int n = 75;

	// Dimension of operator matrix and vectors is n^2
	int N = n*n;

	printf("Solving Poisson Equation on %d x %d domain...\n", n, n);

	/////////////////////////////////////////////////////////////////
	// Run "naive" CG solve, where A is fully stored in dense format
	/////////////////////////////////////////////////////////////////
	
	printf("Running Dense CG solve...\n");
	
	// Allocate full A matrix and vectors
	double ** A = matrix( N );
	double *  x = (double*) calloc(N, sizeof(double));
	double *  b = (double*) calloc(N, sizeof(double));
	printf("Dense Memory  = %.2lf MB\n", (N*N+2*N)*sizeof(double)/1024.0/1024.0);

	// Compute elements of 'A' matrix (Poisson Operator)
	fill_A(A, N);

	// Compute elements of boundary condition vector 'b'
	fill_b(b, N);

	// Run Dense CG Solve
	start = get_time();
	cg_dense(A, x, b, N);
	stop = get_time();
	printf("Dense Runtime = %.2lf seconds\n", stop-start);

	// Save Solution Vector to File
	save_vector(x,N, "dense.out");

	// Free A matrix
	matrix_free(A);

	/////////////////////////////////////////////////////////////////
	// Run optimized CG solve, where A is assembled on the fly (OTF)
	/////////////////////////////////////////////////////////////////
	
	printf("Running Sparse CG solve...\n");
	
	// reset vectors
	memset(x, 0, N * sizeof(double));
	memset(b, 0, N * sizeof(double));
	printf("Sparse Memory  = %.2lf MB\n", (2*N)*sizeof(double)/1024.0/1024.0);

	// Compute elements of boundary condition vector 'b'
	fill_b(b, N);

	// Run Sparse CG Solve
	start = get_time();
	cg_sparse_poisson(x, b, N);
	stop = get_time();
	printf("Sparse Runtime = %.2lf seconds\n", stop-start);

	// Save Solution Vector to File
	save_vector(x,N, "sparse.out");

	// Free vectors
	free(x);
	free(b);
	
	return 0;
}
