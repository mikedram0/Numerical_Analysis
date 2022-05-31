#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// max abs value of a(row,:), b(row)
double maxval(int n, double const a[], double const b[], int row)
{
  int mj = 0;
    
  for (int j = 1; j < n; ++j) {
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

void pivot(int k, int n, double a[], double b[])
{
  int mi = k;
  double mvmi = maxval(n,a,b,mi);
  for (int i = k+1; i < n; ++i) {
    double  mvi = maxval(n,a,b,i);
    
    if (fabs(a[i+k*n] / mvi) > fabs(a[mi+k*n] / mvmi) ) {
      mi = i;
      mvmi = mvi;
    }
  }
  
  for (int j=k; j < n; ++j) {
    swap(&a[mi+j*n], &a[k+j*n]);
  }
  swap(&b[mi], &b[k]);
  
}



void triang(int n, double a[], double  b[])
{
  for (int k = 0; k < n-1; ++k) {

    pivot(k,n,a,b);
	
    for (int i = k+1; i < n; ++i) {
      double const g = -a[i+n*k] / a[k+n*k];
      for (int j = k; j < n; ++j) {
	a[i+j*n] += g * a[k +j*n];
      }
      b[i] += g * b[k];
    }
  }
}


void backsub(int n, double const a[], double const b[], double x[])
{
  for (int i = n-1; i >= 0; --i) {
    double s = b[i];
    for (int j = i+1; j < n; ++j) {
      s -= a[i+j*n] * x[j];
    }
    x[i] = s / a[i+i*n];
  }	
}


void gauss(int n, double a[], double b[], double x[])
{
  triang(n, a, b);
  backsub(n,a,b,x);
}


double pol(int n, double const a[], double xbar)
{
    double r = 0.0;
    for (int i = 0; i < n; ++i) {
	r += a[i] * pow(xbar, i);
    }
    return r;
}

void minmax(double const x[], int n, double * minp, double * maxp)
{
  *minp = x[0];
  *maxp = x[0];

  for (int i = 1; i < n; ++i) {
    if (*minp > x[i]) {
      *minp = x[i];
    }
    if (*maxp < x[i]) {
      *maxp = x[i];
    }
  }
}


//     Read data file
void read_data(char const fname[], double *x[], double *y[], int * np)
{
  FILE * in = fopen(fname, "r");

  fscanf(in, "%d", np);
  *x = malloc(*np * sizeof(double));
  *y = malloc(*np * sizeof(double));
   
  for (int i = 0; i < *np; ++i) {
    fscanf(in, "%lf%lf", &(*x)[i], &(*y)[i]);
  }

  fclose(in);
}




int main()
{
  int n;
  double *x, *y;
  read_data("points.dat", &x,&y,&n);

  double * a = malloc(n * n * sizeof(double) );

  for (int i = 0; i<n; ++i) {
    for (int j = 0; j<n; ++j) {
      a[i+j*n] = pow(x[i], j);
    }
  }

  double * coef = malloc(    n * sizeof(double) );    

  gauss(n, a, y, coef);

    int const m = 100;

    double minv, maxv;
    minmax(x,n, &minv, &maxv);
    double const step = (maxv-minv)/(m-1);

  for (int i = 0; i < m; ++i) {
    double xout = minv + step * i;
    double yout = pol(n,coef,xout);

    printf("%g %g %g\n", xout, yout, sin(xout));
  }

  free(coef);
  free(x);
  free(y);
  free(a);
    
  return 0;
}


