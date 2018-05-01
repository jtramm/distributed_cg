#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

void print_vector(double * x, long N )
{
	long n = sqrt(N);
	long idx = 0;
	for( long i = 0; i < n; i++ )
	{
		for( long j = 0; j < n; j++ )
		{
			printf("%4.1lf ", x[idx]);
			idx++;
		}
		printf("\n");
	}
}

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
	memset( b, 0, n*sizeof(double));
	
	for (long i = 0; i < n; i++)
		for (long j = 0; j < n; j++)
			b[i] += (A[i][j] * x[j]);
}

// "On the fly" Matrix Vector Product for b = A * x
// where 'A' is the known poisson operator matrix
// that is applied without explicit storage
void matvec_OTF( double * b, double * x, long N )
{
	long n = sqrt(N);

	// Wipe solution vector clean
	memset( b, 0, N*sizeof(double));
	
	for( long i = 0; i < N; i++ )
	{
		long x_i = i%n;

		// Far left diagonal
		double far_left_diag = 0.0;
		if( i >= n )
			far_left_diag = -1.0 * x[i-n];

		// left diagonal
		double left_diag = 0.0;
		if( x_i != 0 )
			left_diag = -1.0 * x[i-1];

		// Main diagonal
		double main_diag = 4.0 * x[i];

		// Right diagonal
		double right_diag = 0.0;
		if( x_i != n-1 )
			right_diag = -1.0 * x[i+1];

		// Far right diagonal
		double far_right_diag = 0.0;
		if( i < N - n )
			far_right_diag = -1.0 * x[i+n];

		b[i] = far_left_diag + left_diag + main_diag + right_diag + far_right_diag;
	}
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

void conjgrad(double ** A, double * x, double * b, long N)
{
	// r = -A*x + b
	double * r = (double *) malloc( N*sizeof(double));
	matvec_OTF(r, x, N);
	axpy(-1.0, r, 1.0, b, N);
	//printf("r = p at beginning:\n");
	//print_vector(r, N);

    //p = r;
	double * p = (double *) malloc( N*sizeof(double));
	memcpy(p, r, N*sizeof(double));

    //rsold = r' * r;
	double rsold = dotp(r, r, N); 

	// Ap
	double * Ap = (double *) malloc( N*sizeof(double));

	long iter = 0;
	for( iter = 0; iter < N; iter++ )
	{
        //Ap = A * p;
		matvec_OTF(Ap, p, N);
		//printf("Ap after iter %ld\n", iter);
		//print_vector(Ap, N);
		//break;
	
        //alpha = rsold / (p' * Ap);
		double alpha = rsold / dotp(p, Ap, N);
		//printf("alpha = %lf after iter %ld\n", alpha, iter);

        //x = x + alpha * p;
		axpy(1.0, x, alpha, p, N);

        //r = r - alpha * Ap;
		axpy(1.0, r, -alpha, Ap, N);
		//printf("r after iter %ld\n", iter);
		//print_vector(r, N);

        double rsnew = dotp(r,r,N);

        if( sqrt(rsnew) < 1.0e-10 )
            break;

        //p = (rsnew / rsold) * p + r;
		axpy(rsnew/rsold, p, 1.0, r, N);

        //rsold = rsnew;
		rsold = rsnew;
		//printf("x after iter %ld\n", iter);
		//print_vector(x, N);
	}
	printf("converged in %ld iterations.\n", iter);

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

	printf("A matrix:\n");
	print_matrix(A, N);

	printf("b vector:\n");
	print_vector(b, N);
	
	printf("Running CG...\n");
	conjgrad(A, x, b, N);

	printf("Solution matrix:\n");
	print_vector(x, N);

	
	/*
	// print
	char buf[128];
	sprintf(buf, "A.dat", s);
	outfile = fopen(buf, "wb");
	fwrite(&b[0],sizeof(double),n*n,outfile);
	fclose(outfile);
	*/
}
