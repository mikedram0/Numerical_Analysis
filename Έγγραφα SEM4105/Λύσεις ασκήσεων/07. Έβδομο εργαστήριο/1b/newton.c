#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void create(char const fname[])
{
  int const n = 15;
  double const a=2.0, b=4.0, step=(b-a)/(n-1);

  FILE * out = fopen(fname, "w");    
  fprintf(out, "%d\n", n);
    
  for (int i = 0; i < n; ++i) {
    double x = a + step * i;
    double y = sin(x);
    fprintf(out, "%f %f\n", x, y);
  }

  fclose(out);
}



//     Read data file
void read_data(char const fname[], double *x[], double *y[], int * np)
{
  FILE * in = fopen(fname, "r");

  fscanf(in, "%d", np);
  *x = malloc(*np * sizeof(double));
  *y = malloc(*np * sizeof(double));
   
  for (int i = 0; i < *np; ++i) {
    fscanf(in, "%lf%lf", &((*x)[i]), &((*y)[i]));
  }

  fclose(in);
}


double q(int i, double const x[], double z)
{
  double r = 1.0;
  for (int j = 0; j < i; ++j) {
    r *= z-x[j];
  }

  return r;
}


void coefficients(double const x[], double const y[], double a[], int n)
{
  for (int i = 0; i < n; ++i) {
    double s = y[i];
    for (int j = 0; j < i; ++j) {
      s -= a[j] * q(j,x,x[i]);
    }
    a[i] = s / q(i, x, x[i]);
  }
}


double poly(int n, double const a[], double const x[], double z)
{
  double r = 0.0;
  for (int i=0; i < n; ++i) {
    r += a[i] * q(i,x,z);
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

int main(void)
{
  create("points.dat");
  
  int n;
  double *x, *y;
  read_data("points.dat", &x,&y,&n);

  int const m = 100;

  double minv, maxv;
  minmax(x,n, &minv, &maxv);
  double const step = (maxv-minv)/(m-1);

  double * a = malloc(n*sizeof(double));
  coefficients(x,y,a,n);
    
  for (int i = 0; i < m; ++i) {
    double xout = minv + step * i;
    double yout = poly(n,a,x, xout);

    printf("%g %g\n", xout, yout);
  }
  free(x);
  free(y);
  free(a);

  return 0;
}
