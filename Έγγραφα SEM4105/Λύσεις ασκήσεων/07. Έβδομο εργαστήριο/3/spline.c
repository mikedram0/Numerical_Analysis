#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>



void normalize(size_t n, double a[], double b[], size_t k)
{
  for (size_t j = k; j < n; ++j) {
    double maxv = fabs(b[j]);
    for (size_t p = k; p < n; ++p) {
      double v = fabs(a[j+n*p]);
      if (v > maxv) {
	maxv = v;
      }
    }
	
    for (size_t p = k; p < n; ++p) {
      a[j+n*p] /= maxv;
    }
    b[j] /= maxv;
  }
}



void swap(size_t n, double a[], double b[], size_t k, size_t p)
{
  double t = b[p];
  b[p] = b[k];
  b[k] = t;
  
  for (size_t  j = k; j < n; ++j) {
    t = a[p+n*j];
    a[p+n*j] = a[k+n*j];
    a[k+n*j] = t;
  }
}



void pivot(size_t n, double a[], double b[], size_t k)
{
  normalize(n,a,b,k);

  size_t maxp = k;
  for (size_t p = k+1; p < n; ++p) {
    if (fabs(a[p+n*k]) > fabs(a[maxp+n*k])) {
      maxp = p;
    }
  }

  if (maxp != k) {
    swap(n,a,b,k,maxp);
  }
}



void triang(size_t n, double a[], double b[])
{    
  for (size_t k = 0; k < n-1; ++k) {
    pivot(n,a,b,k);
    
    for (size_t i = k+1; i < n; ++i) {
      double ell = -a[i+n*k]/a[k+n*k];
      for (size_t j = k; j < n; ++j) {
	a[i+n*j] += ell * a[k+n*j];
      }
      b[i] += ell * b[k];
    }
  }
}


void backsub(size_t n, double a[], double b[], double x[])
{
  for (size_t i = n; i-- > 0; ) {
    x[i] = b[i];
    for (size_t j = i+1; j < n; ++j) {
      x[i] -= a[i+n*j] * x[j];
    }
    x[i] /= a[i+n*i];
  }
}


void gauss(size_t n, double A[], double b[], double x[])
{
  triang(n,A,b);
  backsub(n,A,b,x);
}


/*

  Για i=0,...,n-1 :   
  dx_i = x_{i+1}-x_i

  Για i=0,...,n-1 :

  a_i dx_i^2 + b_i dx_i + c_i = (f_{i+1} - f_i) / dx_i


  Για i=0,...,n-2 :

  3 a_i dx_i^2 + 2b_i dx_i + c_i - c_{i+1} = 0


  Για i=0,...,n-2 :

  3a_i dx_i + b_i - b_{i+1} = 0



  b_0 = 0

  3 a_{n-1} dx_{n-1} + b_{n-1} = 0

*/
void fill(size_t n, double const x[], double const y[], double A[], double b[])
{
  size_t const m=3*n;
  for (size_t k=0; k < m*m; ++k) {
    A[k] = 0.0;
  }

  for (size_t k=0; k < m; ++k) {
    b[k] = 0.0;
  }
  
  double * dx = malloc(n*sizeof(double));

  for (size_t i=0; i < n; ++i) {
    dx[i] = x[i+1]-x[i];
  }
    
  for (size_t i=0; i < n; ++i) {
    size_t I = i;
    size_t J = 3*i;

    A[I+m*J]     = dx[i]*dx[i];
    A[I+m*(J+1)] = dx[i];
    A[I+m*(J+2)] = 1.0;

    b[I] = (y[i+1]-y[i])/dx[i];
  }

  for (size_t i=0; i < n-1; ++i) {
    size_t const I = n+i;
    size_t const J = 3*i;

    A[I+m*J]     = 3.0*dx[i]*dx[i];
    A[I+m*(J+1)] = 2.0*dx[i];
    A[I+m*(J+2)] = 1.0;
    A[I+m*(J+5)] = -1.0;
  }
    

  for (size_t i=0; i < n-1; ++i) {
    size_t const I = 2*n-1+i;
    size_t const J = 3*i;
	
    A[I+m*J]     = 3.0*dx[i];
    A[I+m*(J+1)] = 1.0;
    A[I+m*(J+4)] = -1.0;

  }

  A[(m-2)+m*1] = 1.0;
    
  A[(m-1)+m*(m-1-2)] = 3.0*dx[n-1];
  A[(m-1)+m*(m-1-1)] = 1.0;

   
  free(dx);
}



void spline_eval(size_t n, double const x[], double const y[],
		 double sol[])
{
  size_t m = 3*n;

  double *A = malloc(m*m*sizeof(double));
  double *b = malloc(m*sizeof(double));
    
  fill(n,x,y,A,b);
  
  gauss(m,A,b,sol);
  
  free(b);
  free(A);
}



double spline(size_t n, double const x[],
	      double const y[], double const sol[], double v)
{
  assert((v >= x[0]) && (v <= x[n]));

  size_t i=0;
  while (v > x[i]) {
    ++i;
  }
  if (i!=0) {
    --i;
  }
  
  
  double ca = sol[3*i];
  double cb = sol[3*i+1];
  double cc = sol[3*i+2];
  
  double vx = v-x[i];
  return ((ca*vx+cb)*vx+cc)*vx+y[i];
}


int main()
{
  size_t n;
  FILE * in = fopen("points.dat", "r");

  fscanf(in, "%ld", &n);

  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
    
  for (size_t i = 0; i < n; ++i) {
    fscanf(in, "%lg%lg", x+i, y+i);
  }
  --n;

  fclose(in);
  
  double *sol = malloc(3*n * sizeof(double));
	
  spline_eval(n,x,y,sol);

  size_t const m = 15; 

  double z = x[0];
  double const h = (x[n]-x[0]) / (m-1);

    
  for (size_t i=0; i < m; ++i) {
    if (i==m-1) z = x[n];
	
    double w = spline(n,x,y,sol, z);

    printf("%.15lg %.15lg\n", z, w);
	
    z+=h;
  }

  free(sol);
  free(y);
  free(x);

  return 0;
    
}






