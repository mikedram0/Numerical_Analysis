#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


bool check(size_t n, double const y[], double tol)
{
  bool found = true;
  for (size_t i=0; i < n; ++i) {
    found = found && (fabs(y[i]) < tol);
  }

  return found;
}


void gauss_jacobi(size_t n, double const a[], double const b[], double x[])
{
  for (size_t i = 0; i < n; ++i) {
    x[i] = 0.0;
  }

  double * y = malloc(n * sizeof(double));
      
  while (1) {
    for (size_t i=0; i < n; ++i) {
      y[i] = b[i];
      for (size_t j= 0; j < n; ++j) {
	y[i] -= a[i+n*j] * x[j];
      }
      y[i] /= a[i+i*n];
    }
	
    for (size_t i=0; i < n; ++i) {
      x[i] += y[i];
    }

    if (check(n,y,1e-7)) {
      break;
    }
  }

  free(y);
}



void gauss_seidel(size_t n, double const a[], double const b[], double x[])
{
  for (size_t i = 0; i < n; ++i) {
    x[i] = 0.0;
  }
  
  double * y = malloc(n * sizeof(double));
  
  while (1) {
    for (size_t i=0; i < n; ++i) {
      y[i] = b[i];
      for (size_t j= 0; j < n; ++j) {
	y[i] -= a[i+n*j] * x[j];
      }
      y[i] /= a[i+i*n];

      x[i] += y[i];
    }

    if (check(n,y,1e-7)) {
      break;
    }
  }

  free(y);
}


int main()
{
  size_t const n = 4;

  double * a = malloc(n*n * sizeof(double));

  a[0+n*0] = 12.1;
  a[0+n*1] = 3.9;
  a[0+n*2] = 0.3;
  a[0+n*3] = -4.1;

  a[1+n*0] = 4.3;
  a[1+n*1] = -11.3;
  a[1+n*2] = 0.8;
  a[1+n*3] = 1.5;

  a[2+n*0] = 1.0;
  a[2+n*1] = -2.8;
  a[2+n*2] = 14.3;
  a[2+n*3] = -8.1;
    
  a[3+n*0] = 2.4;
  a[3+n*1] = 6.1;
  a[3+n*2] = -1.1;
  a[3+n*3] = 12.5;


  double * b = malloc(n * sizeof(double));
    
  b[0] = 1.2;
  b[1] = 2.3;
  b[2] = 3.4;
  b[3] = 4.5;

  double * x = malloc(n * sizeof(double));


  gauss_jacobi(n,a,b,x);
    
  printf("solution with jacobi:\n");
  for (size_t i = 0; i < n; ++i) {
    printf("%g ", x[i]);
  }
  printf("\n");

    
  gauss_seidel(n,a,b,x);
    
  printf("solution with seidel:\n");
  for (size_t i = 0; i < n; ++i) {
    printf("%g ", x[i]);
  }
  printf("\n");

  free(x);
  free(b);
  free(a);
  
  return 0;
}
