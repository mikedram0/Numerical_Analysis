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
    fscanf(in, "%lf%lf", &(*x)[i], &(*y)[i]);
  }

  fclose(in);
}


double ell(int i, double z, double const x[], int n)
{
  double p = 1.0;
  for(int j = 0; j < n; ++j) {
    if (j==i) {
      continue;
    }

    p *= (z-x[j]) / (x[i]-x[j]);
  }

  return p;
}


double poly(int n, double const x[], double const y[], double z)
{
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    s += ell(i,z,x,n) * y[i];
  }

  return s;
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
  double *x;
  double *y;
  read_data("points.dat", &x,&y, &n);

  int const m = 100;

  double minv, maxv;
  minmax(x,n, &minv, &maxv);
  double const step = (maxv-minv)/(m-1);

  for (int i = 0; i < m; ++i) {
    double xout = minv + step * i;
    double yout = poly(n,x,y, xout);

    printf("%g %g\n", xout, yout);
  }

  free(x);
  free(y);
}
