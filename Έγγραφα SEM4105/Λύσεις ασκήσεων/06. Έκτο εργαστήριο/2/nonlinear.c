#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* ---------------- GAUSS   begin ------------------ */

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

void normalize(double a[], double b[], size_t k, size_t n)
{
  for (size_t j = k; j < n; ++j) {
    double maxv = fabs(b[j]);
    for (size_t p = k; p < n; ++p) {
      size_t const v = fabs(a[j+n*p]);
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

void swaparr(double a[], double b[], size_t n, size_t k, size_t p)
{
  for (size_t j = k; j < n; ++j) {
    swap(&a[p + n*j], &a[k + n *j]);
  }
  swap(&b[p], &b[k]);
}


void pivot(size_t k, size_t n, double a[], double b[])
{
  normalize(a,b,k,n);
  
  size_t maxp = k;
  for (size_t p = k+1; p < n; ++p) {
    if (fabs(a[p+n*k]) > fabs(a[maxp+n*k])) {
      maxp = p;
    }
  }

  if (maxp != k) {
    swaparr(a,b,n,k,maxp);
  }

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


/* ---------------- GAUSS  end  ------------------ */


void values(double const x[], double f[])
{
    f[0] = x[0]             + x[1]             + x[2] - 3.0;
    f[1] = x[0]*x[0] * x[1] + x[1]*x[1] * x[2] + x[2]*x[2] * x[0] - 4.0;
    f[2] = x[0]*x[0]        + x[1]*x[1]        + x[2]*x[2]  - 5.0;
}

void derivs(int n, double const x[], double a[])
{
    a[0+n*0] = 1.0;
    a[0+n*1] = 1.0;
    a[0+n*2] = 1.0;

    a[1+n*0] = 2.0 * x[0] * x[1] + x[2]*x[2];
    a[1+n*1] = 2.0 * x[1] * x[2] + x[0]*x[0];
    a[1+n*2] = 2.0 * x[2] * x[0] + x[1]*x[1];

    a[2+n*0] = 2.0 * x[0];
    a[2+n*1] = 2.0 * x[1];
    a[2+n*2] = 2.0 * x[2];
}


int check(int n, double const z[], double const eps[])
{
  int k = 1;
  for (int i = 0; i < n; ++i) {
    if (fabs(z[i]) > eps[i]) {
      k = 0;
      break;
    }
  }

  return k;
}


void output(int n, double const x[])
{
  for (int i=0; i < n; ++i) {
    printf("%f ", x[i]);
  }
  printf("\n");
}



int main()
{
  size_t const n = 3;

  double * eps = malloc(n*sizeof(double));
  double * x   = malloc(n*sizeof(double));

  double * b   = malloc(n*sizeof(double));
  double * y   = malloc(n*sizeof(double));
  double * a   = malloc(n*n*sizeof(double));

  x[0] = 1.0;
  x[1] = 1.5;
  x[2] = 2.1;
  
  for (size_t i = 0; i < n; ++i) {
    eps[i] = 1e-7;
  }
  
  while (1) {
    values(x,b);
    if (check(n,b,eps)) {
      break;
    }
    derivs(n,x,a);
    gauss(n,a,b,y);
    
    for (size_t i=0; i < n; ++i) {
      x[i] -= y[i];
    }
  }

  output(n,x);

  free(eps);
  free(x);
  free(a);
  free(b);
  free(y);
  
  return 0;
}
