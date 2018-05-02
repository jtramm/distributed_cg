#include "conj_header.h"

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

void save_vector(double * x, long N )
{
	FILE * fp = fopen("out.txt", "w");
	long n = sqrt(N);
	long idx = 0;
	for( long i = 0; i < n; i++ )
	{
		for( long j = 0; j < n; j++ )
		{
			fprintf(fp, "%.9le", x[idx]);
			idx++;
			if( j != n-1 )
				fprintf(fp, " ");
		}
		if( i != n - 1 )
			fprintf(fp, "\n");
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
