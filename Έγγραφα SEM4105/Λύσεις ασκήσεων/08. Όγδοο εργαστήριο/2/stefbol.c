#include <math.h>
#include <stdio.h>
#include <stdlib.h>



double mylog(double x)
{
  return log(x);
}



//     Read data file
void read_data(char const * fname, double **x, double **y, int * n)
{
  FILE *in = fopen(fname, "r");

  fscanf(in,"%d\n", n);

  *x = malloc(*n * sizeof(double));
  *y = malloc(*n * sizeof(double));
  
  for (int i=0; i < *n; ++i) {
    double a, b;
    fscanf(in,"%lf%lf", &a, &b);
    (*x)[i] = a;
    (*y)[i] = b;
  }
  fclose(in);
}



//     Y = A1 * X + a0
void lsqlin(int n, double const x[], double const y[],
	    double * a1, double * a0, double * r2)
{
  double sx = 0.0;
  double sy = 0.0;
  double sx2 = 0.0;
  double sy2 = 0.0;
  double sxy = 0.0;

  for (int i=0; i < n; ++i) {
    sx += x[i];
    sy += y[i];
    sx2 += x[i]*x[i];
    sy2 += y[i]*y[i];
    sxy += x[i]*y[i];
  }

  double const d = n * sx2 - sx * sx;
  
  *a1 = (n * sxy - sx * sy) / d; 
  *a0 = (sy - *a1 * sx) / n;
  *r2 = pow(n * sxy - sx * sy,2) / d / (n * sy2 - sy*sy);
}


//     Y = a0 X^A1 => LN(Y) = LN(a0) + A1 * LN(X)
void lsqpow(int n, double const x[], double const y[],
	    double * a1, double * a0, double * r2)
{
  for (int i = 0; i < n; ++i) {
    if ((x[i] <= 0.0) || (y[i] <= 0.0)) {
      fprintf(stderr, "Data not appropriate\n");
      return; 
    }
  }
  
  double *zx = malloc(n * sizeof(double));
  double *zy = malloc(n * sizeof(double));

  double za0;

  for (int i=0; i < n; ++i) {
    zx[i] = mylog(x[i]);
    zy[i] = mylog(y[i]);
  }

  lsqlin(n,zx, zy, a1, &za0, r2);
  
  *a0 = exp(za0);

  free(zx);
  free(zy);
}


int
main()
{
  int n;
  double *x, *y;
  read_data("stefbol.dat", &x, &y, &n);
  
  double a1, a0, r2;
  lsqpow(n,x, y, &a1, &a0, &r2);
    
  double const embado = 0.05e-4;
  printf("Σταθερά Stefan-Boltzmann: %g J/K^4/m^2/s\n",  a0 / embado);
  printf("Εκθέτης: %g\n", a1);
  printf("r^2 = %g\n", r2);

  free(x);
  free(y);
  
}


