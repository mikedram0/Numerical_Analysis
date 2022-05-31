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



/*     Read data file */
void read_data(char const * fname, double **x, double **y, int * nptr)
{
  FILE *in = fopen(fname, "r");

  fscanf(in,"%d\n", nptr);

  *x = malloc(*nptr * sizeof(double));
  *y = malloc(*nptr * sizeof(double));
  
  for (int i = 0; i < *nptr; ++i) {
    double a, b;
    fscanf(in,"%lf%lf", &a, &b);
    (*x)[i] = a;
    (*y)[i] = b;
  }
  fclose(in);
}



/*  f^{(order)}_i (x) */
double f(int i, int order, double x)
{
  double v;
    
  switch (order) {
  case 0:    /* => function value */
    v = pow(x,i);
    break;
	
  case 1:    /* => first derivative */
    if (i == 0) {
      v = 0.0;
    } else {
      v= i * pow(x,i-1);
    }
    break;

  case 2: /* => second derivative */
    if (i == 0 || i == 1) {
      v = 0.0;
    } else {
      v= i*(i-1)*pow(x,i-2);
    }
    break;
  default:
    fprintf(stderr, "%s", "not implemented\n");
    exit(-1);
      
  }

  return v;
}


void sys(int n, double const x[], double xbar, int order, double a[],
	 double b[])
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i+j*n] = f(i, 0, x[j]);
    }
    b[i] = f(i,order,xbar);
  }

}

double correct(int order, double x)
{
  switch (order) {
  case 0:
    return pow(sin(x),2);
  case 1:
    return 2*sin(x)*cos(x);
  case 2:
    return 2*cos(2*x);
  default:
    fprintf(stderr, "%s", "not implemented\n");
    exit(-1);
  }

}

int main(void)
{
  int n;
  double *x, *y;
  read_data("points2.dat", &x, &y, &n);

  double *a = malloc(n*n*sizeof(double));
  double *b = malloc(n*sizeof(double));
  double *sol = malloc(n*sizeof(double));

  double xbar = 2.5;
    
  for (int order = 1; order <= 2; ++order) {
    sys(n, x, xbar, order, a, b);
	
    gauss(n,a,b,sol);
	
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
      s += sol[i] * y[i];
    }
    
    free(x);
    free(y);
    free(sol);
    free(b);
    free(a);
    
    printf("%d derivative at %g is %g\n", order, xbar, s);
    printf("difference from correct value %g\n", s - correct(order,xbar));
  }

  return 0;
}
