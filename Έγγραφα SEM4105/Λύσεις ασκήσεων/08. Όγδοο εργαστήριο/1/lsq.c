#include <math.h>
#include <stdio.h>
#include <stdlib.h>



double mylog(double x)
{
  return log(x);
}

double myexp(double x)
{
  return exp(x);
}



/*     Read data file */
void read_data(char const * fname, double **x, double **y, int * nptr)
{
  FILE *in = fopen(fname, "r");

  fscanf(in,"%d\n", nptr);

  *x = malloc(*nptr * sizeof(double));
  *y = malloc(*nptr * sizeof(double));
  
    for (int i=0; i < *nptr; ++i) {
	double a, b;
	fscanf(in,"%lf%lf", &a, &b);
	(*x)[i] = a;
	(*y)[i] = b;
    }
    fclose(in);
}



/*     Y = A1 * X + a0 */
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


/*     Y = a0 X^A1 => LN(Y) = LN(a0) + A1 * LN(X) */
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


/*     Y = a0 + a1 * LOG(X)   */
void lsqlog(int n, double const x[], double const y[],
	    double * a1, double * a0, double * r2)
{
  for (int i = 0; i < n; ++i) {
    if (x[i] <= 0.0) {
      fprintf(stderr, "Data not appropriate\n");
      return; 
    }
  }
  
  double *zx = malloc(n * sizeof(double));

  for (int i=0; i < n; ++i) {
    zx[i] = mylog(x[i]);
  }

  lsqlin(n,zx, y, a1, a0, r2);
  
  free(zx);
}



/*     Y = a0 + a1 * EXP(X) */
voidlsqexp(int n, double const x[], double const y[],
	   double * a1, double * a0, double * r2)
{
  double *zx = malloc(n * sizeof(double));

  for (int i=0; i < n; ++i) {
    zx[i] = myexp(x[i]);
  }

  lsqlin(n,zx, y, a1, a0, r2);
  
  free(zx);
}

int main(void)
{
  int n;
  double *x, *y;
  read_data("points.dat", &x, &y, &n);

  printf("Μέθοδος:\n");
  printf("1) y = a * x + b\n");
  printf("2) y = a * x^b\n");
  printf("3) y = a + b * exp(x)\n");
  printf("4) y = a + b * ln(x)\n");
  printf("Δώσε μέθοδο.\n");

  int method;
  scanf("%d", &method);
  
  double a1, a0, r2;
  
  switch (method) {
  case 1:
    lsqlin(n, x, y, &a1, &a0, &r2);
    printf("y = %g * x + %g\n", a1, a0);
    break;
  case 2:
    lsqpow(n,x, y, &a1, &a0, &r2);
     printf("y = %g * x^%g\n", a0, a1);
    break;
  case 3:
    lsqexp(n,x, y, &a1, &a0, &r2);
    printf("y = %g + %g * exp(x)\n", a0, a1);
    break;
  case 4:
    lsqlog(n,x, y, &a1, &a0, &r2);
    printf("y = %g + %g * ln(x)\n", a0, a1);
    break;
  }

  printf("r^2 = %g\n", r2);

  free(x);
  free(y);
  
}


