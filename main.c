#include "cg_header.h"

int main(int argc, char* argv[])
{
	int nprocs = 1;
	int mype = 0;

	#ifdef MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	// n is the physical domain size;
	int n;

	// Read n from first command line argument
	if( argc != 2 )
	{
		printf("Please provide physical domain dimension as first argument to program, e.g.:\n\t$> ./cg 75\n");
		return 1;
	}
	else
		n = atoi(argv[1]);
	
	// Make sure n is reasonable
	assert(n > 0 && isfinite(n));

	printf("Solving Poisson Equation on %d x %d domain...\n", n, n);

	// Run dense CG solve
	run_dense(n);

	// Run sparse CG solve
	run_sparse(n);

	// Run parallel CG solve
	run_parallel_sparse(n, mype, nprocs);

	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}

/////////////////////////////////////////////////////////////////
// Run "naive" CG solve, where A is fully stored in dense format
/////////////////////////////////////////////////////////////////
void run_dense(long n)
{
	printf("Running Dense CG solve...\n");

	// Dimension of operator matrix and vectors is n^2
	int N = n*n;

	// Allocate full A matrix and vectors
	double ** A = matrix( N );
	double *  x = (double*) calloc(N, sizeof(double));
	double *  b = (double*) calloc(N, sizeof(double));
	printf("Dense Memory  = %.2lf MB\n", (N*N+5*N)*sizeof(double)/1024.0/1024.0);

	// Compute elements of 'A' matrix (Poisson Operator)
	fill_A(A, N);

	// Compute elements of boundary condition vector 'b'
	fill_b(b, N);

	// Run Dense CG Solve
	double start = get_time();
	cg_dense(A, x, b, N);
	double stop = get_time();
	printf("Dense Runtime = %.2lf seconds\n", stop-start);

	// Save Solution Vector to File
	save_vector(x,N, "dense.out");

	// Free A matrix
	matrix_free(A);

	// Free vectors
	free(x);
	free(b);
}

/////////////////////////////////////////////////////////////////
// Run optimized CG solve, where A is assembled on the fly (OTF)
/////////////////////////////////////////////////////////////////
void run_sparse(long n)
{
	printf("Running Sparse CG solve...\n");

	// Dimension of operator matrix and vectors is n^2
	int N = n*n;

	// reset vectors
	double *  x = (double*) calloc(N, sizeof(double));
	double *  b = (double*) calloc(N, sizeof(double));
	printf("Sparse Memory  = %.2lf MB\n", (5*N)*sizeof(double)/1024.0/1024.0);

	// Compute elements of boundary condition vector 'b'
	fill_b(b, N);

	// Run Sparse CG Solve
	double start = get_time();
	cg_sparse_poisson(x, b, N);
	double stop = get_time();
	printf("Sparse Runtime = %.2lf seconds\n", stop-start);

	// Save Solution Vector to File
	save_vector(x,N, "sparse.out");

	// Free vectors
	free(x);
	free(b);
}

/////////////////////////////////////////////////////////////////
// Run Domain Decomposed Parallel CG solve w/MPI
// where A is assembled on the fly (OTF)
/////////////////////////////////////////////////////////////////
void run_parallel_sparse(long n, int mype, int nprocs )
{
	printf("Parallel Sparse Solver Not Yet Implemented...\n");

}
