#include "conj_header.h"

int main(int argc, const char** argv){
	
	FILE *outfile;
		
	// n x n square grid
	int n = 100;

	// Dimension of operator matrix and vectors is n^2
	int N = n*n;
	
	// Allocate full A matrix and vectors
	double ** A = matrix( N );
	double *  x = (double*) calloc(N, sizeof(double));
	double *  b = (double*) calloc(N, sizeof(double));

	// Compute elements of 'A' matrix (Poisson Operator)
	fill_A(A, N);

	// Compute elements of boundary condition vector 'b'
	fill_b(b, N);

	conjgrad(A, x, b, N);

	save_vector(x,N);

	// Run "naive" CG solve, where A is fully stored in dense format
	
}
