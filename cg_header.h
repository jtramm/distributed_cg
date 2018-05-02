#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// serial.c
void cg_dense(double ** A, double * x, double * b, long N);
void cg_sparse_poisson(double * x, double * b, long N);
void matvec( double * v, double ** M, double * w, long N );
void matvec_OTF( double * v, double * w, long N );
double dotp( double * a, double * b, long N);
void axpy( double alpha, double * w, double beta, double * v, long N);
void fill_A(double ** A, long N);
double find_b(long i, long j, long n);
void fill_b(double * b, long N);

// utils.c
void print_matrix(double ** A, long N );
void print_vector(double * x, long N );
void save_vector(double * x, long N, char * fname );
double ** matrix( long N );
void matrix_free( double ** M);
double get_time(void);
