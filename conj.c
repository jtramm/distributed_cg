#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Allocates 2-D Contiguous Matrix */
double ** matrix( long N )
{
	double *data = (double *) calloc( N*N, sizeof(double) );
	double **M  = (double **) malloc( N  * sizeof(double*));

	for( int i = 0; i < N; i++ )
		M[i] = &data[i*N];

	return M;
}

/* Free's 2-D Contiguous Matrix */
void matrix_free( double ** M)
{
	free(M[0]);
	free(M);
}

// Matrix Vector Product for b = A * x
void matvec( double * b, double ** A, double * x, long n )
{
	// Wipe solution vector clean
	memset( x, 0, n*sizeof(double));
	
	for (long i = 0; i < n; i++)
		for (long j = 0; j < n; j++)
			b[i] += (A[i][j] * x[j]);
}

// dot product of c = a * b
double dotp( double * a, double * b, long n)
{
	double c = 0.0;
	for( long i = 0; i < n; i++ )
		c += a[i]*b[i];

	return c;
}

/* MATLAB Function
function [x] = conjgrad(A, b, x)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end
*/

// solves x = alpha * x + beta * v (overwrites what's in x)
// alpha, beta = scalars
// x, v = 1D vectors
void axpy( double alpha, double * x, double beta, double * v, long n)
{
	for( long i = 0; i < n; i++ )
		x[i] = alpha * x[i] + beta * v[i];
}

void conjgrad(double ** A, double * x, double * b, long n)
{
	// r = -A*x + b
	double * r = (double *) malloc( n*sizeof(double));
	matvec(r, A, x, n);
	axpy(-1.0, r, 1.0, b, n);

    //p = r;
	double * p = (double *) malloc( n*sizeof(double));
	memcpy(p, r, n*sizeof(double));

    //rsold = r' * r;
	double rsold = dotp(r, r, n); 

	// Ap
	double * Ap = (double *) malloc( n*sizeof(double));

	for( long iter = 0; iter < n; iter++ )
	{
        //Ap = A * p;
		matvec(Ap, A, p, n);
	
        //alpha = rsold / (p' * Ap);
		double alpha = rsold / dotp(p, Ap, n);

        //x = x + alpha * p;
		axpy(1.0, x, alpha, p, n);

        //r = r - alpha * Ap;
		axpy(1.0, r, -alpha, Ap, n);

        double rsnew = dotp(r,r,n);

        if( sqrt(rsnew) < 1.0e-10 )
            break;

        //p = (rsnew / rsold) * p + r;
		axpy(rsnew/rsold, p, 1.0, r, n);

        //rsold = rsnew;
		rsold = rsnew;
	}

	free(r);
	free(p);
	free(Ap);
}

void fill_A(double ** A, long N)
{
	long n = sqrt(N);

	for( long i = 0; i < N; i++ )
	{
		for( long j = 0; j < N; j++ )
		{
			// Compute x-coordinate of index
			long x = j % n;

			// main diagonal
			if(i == j)
				A[i][j] = 4.0;
			// left diagonal
			else if (i == j+1 )
			{
				A[i][j] = -1.0;
				if(x == n-1)
					A[i][j] = 0.0;
			}
			// right diagonal
			else if (i == j-1 )
			{
				A[i][j] = -1.0;
				if(x == 0)
					A[i][j] = 0.0;
			}
			// far left diagonal
			else if(j+n == i)
				A[i][j] = -1.0;
			// far right diagonal
			else if(j-n == i)
				A[i][j] = -1.0;
			// Otherwise, A is 0
			else
				A[i][j] = 0.0;
		}
	}
}

double get_b(long i, long j, long n)
{
	double up = 0.0;
	double down = 2.0;
	double left = 0.5;
	double right = 1.0;

	long idx = i*n + j;

	if( j == 0 )
		return left;
	else if( j == n-1 )
		return right;
	else if( i == 0 )
		return up;
	else if( i == n-1 )
		return down;
	
	return 0.0;
}

void fill_b(double * b, long N)
{
	long n = sqrt(N);
	for( long i = 0; i < n; i++ )
		for( long j = 0; j < n; j++ )
		{
			long idx = i*n + j;
			b[idx] = get_b(i,j,n);
		}
}

void print_matrix(double ** A, long N )
{
	for( long i = 0; i < N; i++ )
	{
		for( long j = 0; j < N; j++ )
		{
			printf("%4.1lf ", A[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, const char** argv){
	
	FILE *outfile;
		
	// n x n square grid
	int n = 4;

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

	// Solve for 'x' solution vector
	print_matrix(A, N);
	
	/*
	// print
	char buf[128];
	sprintf(buf, "A.dat", s);
	outfile = fopen(buf, "wb");
	fwrite(&b[0],sizeof(double),n*n,outfile);
	fclose(outfile);
	*/
}
