#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

bool iszero(double x)
{
  return fabs(x) < 1e-10;
}

// maxval in row i, from k column to end
double rowmaxval(size_t n, double const a[], double const b[],
		 size_t k, size_t row)
{
  double maxv = fabs(b[row]);

  for (size_t j = k; j < n; ++j) {
    double const v = fabs(a[row+n*j]);
    if (v > maxv) {
      maxv = v;
    }
  }

  return maxv;
}



void normalize(size_t n, double a[], double b[], size_t k)
{
  for (size_t i = k; i < n; ++i) {
    double const maxv = rowmaxval(n, a, b, k, i);

    if (iszero(maxv)) {
      continue;
    }
	
    for (size_t j = k; j < n; ++j) {
      a[i+n*j] /= maxv;
    }
    b[i] /= maxv;
  }
}


void swap(double * a, double * b)
{
  double t = *a;
  *a = *b;
  *b = t;
}


void swaparr(size_t n, double a[], double b[], size_t k, size_t i)
{
  for (size_t j = k; j < n; ++j) {
    swap(&a[i + n*j], &a[k + n *j]);
  }
  swap(&b[i], &b[k]);
}



void triang(size_t n, double a[], double  b[])
{
  for (size_t k = 0; k < n-1; ++k) {	
    for (size_t i = k+1; i < n; ++i) {
      double const ell = -a[i+n*k] / a[k+n*k];
      a[i+k*n]  = 0.0;
      for (size_t j = k+1; j < n; ++j) {
	a[i+j*n] += ell * a[k +j*n];
      }
      b[i] += ell * b[k];
    }
  }
}


// ax=b
double solve_linear(double a, double b)
{
  if (iszero(a)) {
    if (iszero(b)) {
      fprintf(stderr, "%s", "Το σύστημα έχει άπειρες λύσεις\n");
	} else {
            fprintf(stderr, "%s", "Το σύστημα δεν έχει λύση\n");
    }
  }
  
  return b/a;
}



void backsub(size_t n, double const a[], double const b[], double x[])
{
  for (size_t i = n; i-- > 0; ) {
    x[i] = b[i];
    for (size_t j = i+1; j < n; ++j) {
      x[i] -= a[i+j*n] * x[j];
    }
    x[i] = solve_linear(a[i+n*i], x[i]);
  }	
}


void gauss(size_t n, double a[], double b[], double x[])
{
  triang(n, a, b);
  backsub(n,a,b,x);
}


int main()
{
  FILE * in = fopen("gauss1.txt", "r");
  size_t n;
  fscanf(in, "%zd", &n);
  
  double * a = malloc(n * n * sizeof(double) );
  double * b = malloc(    n * sizeof(double) );
  double * x = malloc(    n * sizeof(double) );
  
  for (size_t i=0; i<n; ++i) {
    for (size_t j=0; j<n; ++j) {
      fscanf(in, "%lf", &a[i+j*n]);
    }
  }

  for (size_t i = 0; i<n; ++i) {
    fscanf(in, "%lf", &b[i]);
  }

  gauss(n,a,b,x);

  for (size_t i = 0; i<n; ++i) {
    printf("%g\n", x[i]);
  }

  free(x);
  free(b);
  free(a);

  return 0;
}
