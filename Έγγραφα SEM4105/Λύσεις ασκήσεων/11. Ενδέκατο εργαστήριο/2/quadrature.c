#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double const pi = 3.14159265358979323846;

// max abs value of a(row,:), b(row)
double maxval(size_t n, double const a[], double const b[], size_t row)
{
  size_t mj = 0;
    
  for (size_t j = 1; j < n; ++j) {
    if (fabs(a[row+j*n]) > fabs(a[row+mj*n])) {
      mj=j;
    }
  }

  double g = fmax(fabs(b[row]), fabs(a[row+mj*n]));
    
  return g;
}

void swap(double * a, double *b)
{
  double t = *a;
  *a = *b;
  *b = t;
}

void pivot(size_t k, size_t n, double a[], double b[])
{
  size_t mi = k;
  double mvmi = maxval(n,a,b,mi);
  for (size_t i = k+1; i < n; ++i) {
    double  mvi = maxval(n,a,b,i);
    
    if (fabs(a[i+k*n] / mvi) > fabs(a[mi+k*n] / mvmi) ) {
      mi = i;
      mvmi = mvi;
    }
  }
  
  for (size_t j=k; j < n; ++j) {
    swap(&a[mi+j*n], &a[k+j*n]);
  }
  swap(&b[mi], &b[k]);
  
}



void triang(size_t n, double a[], double  b[])
{
  for (size_t k = 0; k < n-1; ++k) {

    pivot(k,n,a,b);
	
    for (size_t i = k+1; i < n; ++i) {
      double const g = -a[i+n*k] / a[k+n*k];
      for (size_t j = k; j < n; ++j) {
	a[i+j*n] += g * a[k +j*n];
      }
      b[i] += g * b[k];
    }
  }
}


void backsub(size_t n, double const a[], double const b[], double x[])
{
  for (size_t i = n; i-- > 0; ) {
    double s = b[i];
    for (size_t j = i+1; j < n; ++j) {
      s -= a[i+j*n] * x[j];
    }
    x[i] = s / a[i+i*n];
  }	
}


void gauss(size_t n, double a[], double b[], double x[])
{
  triang(n, a, b);
  backsub(n,a,b,x);
}


/* Gauss elimination with pivoting: end */



void sys(size_t n, double const x[], double a, double b,
	 double A[], double B[])
{
  for (size_t i=0; i < n; ++i) {
	for (size_t j=0; j < n; ++j) {
	    A[i+j*n] = pow(x[j],i); 
	}
	B[i] = (pow(b,i+1)-pow(a,i+1))/(i+1);
    }
}



double f(double x)
{    
  return pow(x,3) * sin(pi * x);
}


int main()
{
    double const x[] = {-0.9,-0.7,-0.4,0.1,0.4,0.8,0.9};
    
    size_t const n = sizeof(x)/sizeof(x[0]);

    double *a = malloc(n*n*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *w = malloc(n*sizeof(double));
    
    double alimit = -1.0;
    double blimit = 1.0;
    
    sys(n, x, alimit, blimit, a, b);

    gauss(n,a,b,w);

    double s = 0.0;
    for (size_t i=0; i < n; ++i) {
	s += w[i] * f(x[i]);
    }

    double const correct = (2.0 -12.0 /(pi*pi)) / pi;

    printf("%g %g\n", s, correct);

    return 0;
}
