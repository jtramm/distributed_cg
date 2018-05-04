#include "cg_header.h"

// Parallel Distributed MPI Version of
// Conjugate Gradient Solver Function for Ax = b
// where A is the 2D poisson matrix operator that
// is assembled on the fly (not stored)
// x = 1D solution vector
// b = 1D vector
// N = dimension
void parallel_cg_sparse_poisson(double * x, double * b, long N, int mype, int nprocs)
{
}

// Parallel Distributed MPI Version of
// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void parallel_matvec_OTF( double * v, double * w, long N, int mype, int nprocs )
{
}

// Parallel Distributed MPI Version of
// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double parallel_dotp( double * a, double * b, long N, int mype, int nprocs)
{
	return 0;
}

// Parallel Distributed MPI Version of
// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void parallel_axpy( double alpha, double * w, double beta, double * v, long N, int mype, int nprocs)
{
}

// Parallel Distributed MPI Version of
// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void parallel_fill_b(double * b, long N, int mype, int nprocs)
{
}
