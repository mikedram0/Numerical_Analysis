#include <complex.h>

typedef double complex COMPL;

double falsepoint( double a , double b , double (*func)(double)  , double epsilon);

double illinois( double a , double b , double (*func)(double)  , double epsilon);

double bisec( double a , double b , double (*func)(double)  , double epsilon);

//double temn( double a, double b, double (*func)(double), double epsilon);
COMPL temn( COMPL a, COMPL b, COMPL (*func)(COMPL), double epsilon);
//COMPL mod_temn( COMPL a, COMPL b, COMPL (*func)(COMPL,COMPL*,int), double epsilon, COMPL* roots, int n);
int matmul(int MA, int NA, double A[], int MB, int NB, double  B[], double result[]);

void copyto( int M, int N, double begin[], double end[]  );

void GS(int N,double A[], double B[], double x[], double epsilon);

void jacobi(int N,double A[], double B[], double x[], double epsilon);

void setmat(int M, int N, double A[],double val);

void setmatc(int M, int N, COMPL A[],COMPL val);

void changediag(int M, COMPL B[], COMPL num );

COMPL det(int N, COMPL A[], COMPL lambda);

COMPL mod_temn_matrix( COMPL a, COMPL b, COMPL (*func)(COMPL,COMPL*,int,int,COMPL*), double epsilon, COMPL* roots, int n, int N, COMPL A[]);

COMPL char_pol (COMPL x, COMPL roots[], int n, int N, COMPL A[]);

void eigenval(int N, COMPL A[], COMPL roots[] , COMPL a, COMPL b );

void eigenvec(int N, COMPL A[], COMPL roots[], COMPL result[], double epsilon);

int maketriangle( int N, double A[]);

void triangle_pivot( int N, double A[], int Nb, double B[]);

int backtrack( int N, double A[], int Nb, double B[], double epsilon);

int gauss( int N, double A[],int Nb, double B[], double epsilon);

int COMPLgauss( int N, COMPL A[],int Nb, COMPL B[], double epsilon);

void print( int M,int N, double A[]);

void printc( COMPL x);

void COMPLprint( int M,int N, COMPL A[]);

COMPL muller( COMPL x0, COMPL x1, COMPL x2, COMPL (*func)( COMPL),double epsilon);

double stablepoint(double a, double  (*func)(double), double epsilon);

double NR(double a, double (*func)( double), double (*dfunc)(double), double  epsilon);

double mean(double* x, int n);

void lsq(double* x, double* y, int n, double *a, double *b);

void fit( double *x, double *y, int n, double *a, double *b, char type);

void polynomial( double *x, double *y, double *res,int M, int N );

void linspace(double a, double b, int N, double *x);

void range(double a, double b, double h, double *x);

int nonlinear(double *xi,  double (*F[])(double*, int), int N, double epsilon);

double polynomial_interpolation( double *x, double *y, int N, double xi, char method);

void polynomialgauss (double* x, double* y, int N, double* xi,int Ni, double* result);

void poldiv(double *x, double *y, int M,int N,int Ni, double *xi, double *result);

double inter2d( double *x, int Nx, double *y, int Ny, double* f, double xi, double yi );