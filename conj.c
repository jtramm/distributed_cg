#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
	
	for (long i = 0; i < n; i++){
		for (long j = 0; j < n; j++){
			b[i] += (A[i][j] * x[j]);
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

int main(int argc, const char** argv){
	
	FILE *outfile;
		
	// n x n square grid
	int n = 65;
	int nsteps = 100;
	double alpha = 0.1;
	double dt = 1;
	double dx = 1;
	double C = (alpha*dt) / (dx*dx);
	
	// stuff
	double *A = (double*)calloc(n*n*n*n, sizeof(double));
	double *x = (double*)calloc(n*n, sizeof(double));
	double *b = (double*)calloc(n*n, sizeof(double));
	
	
	// init grid
    int x0 = n/2;
    int y0 = n/2;
    double a = 1.0f;
    double sig_x = 5.0f;
    double sig_y = 5.0f;
	for (int j = 0; j < n; j++){
		for (int i = 0; i < n; i++){
			b[j * n + i] = a * exp(-1.0f * ((((i-x0)*(i-x0)) / (2*sig_x*sig_x)) + (((j-y0)*(j-y0)) / (2*sig_y*sig_y))));                                                        
		}
	}
	
	// set up A
	for (int i = 0; i < n*n; i++){
		A[i * n*n + i] = 1 + 4*C;
		if (i+1 < n) A[i * n*n + i+1] = -C;
		if (i+n < n*n) A[i * n*n + i+n] = -C;
		if (i > 0) A[i * n*n + i-1] = -C;
		if (i >= n) A[i * n*n + i-n] = -C;
	}
	
	// advanace time
	for (int s = 0; s < nsteps; s++){
		printf("step = %d\n", s);
		
		// reset boundaries, i think this might cause issues with discontinuties
		for (int i = 0; i < n; i++){
			b[0 * n + i] = 0;
			b[(n-1) * n + i] = 0;
			b[i * n + 0] = 0;
			b[i * n + (n-1)] = 0;
		}
		
		// add source term
		for (int j = 0; j < n; j++){
			for (int i = 0; i < n; i++){
				b[j * n + i] += source(i,j);
			}
		}
		
		// print
		char buf[128];
		sprintf(buf, "frame%d.dat", s);
		outfile = fopen(buf, "wb");
		fwrite(&b[0],sizeof(double),n*n,outfile);
		fclose(outfile);
		
		// solve and rotate
		solCG(n*n, A, b, x);
		memcpy(b, x, n*n*sizeof(double));
		
	}

}
