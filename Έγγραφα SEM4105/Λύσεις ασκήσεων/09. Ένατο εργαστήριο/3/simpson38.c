#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>


double fint(double x)
{
    return (sin(x) - cos(x)) / 2.0 * exp(x);
}



double simpson38(double const x[], double const y[])
{
    double const h = x[1]-x[0];
    
    return 3.0*h/8.0 * (y[0]+3.0*(y[1] +y[2]) + y[3]);
}


double simpson13(size_t n, double const y[], double a, double b)
{
  --n;
  
  if (n%2 != 0) {
    fprintf(stderr, "error\n");
  }
    
  double const h = (b-a)/n;
  
  double s = y[0] + y[n];

  double s1 = 0.0;
  for (size_t i = 1; i < n; i+=2) {
	s1 += y[i];
    }

    s += 4.0 * s1;

    double s2 = 0.0;
    for (size_t i = 2; i < n; i+=2) {
	s2 += y[i];
    }

    s += 2.0 * s2;
    
    return s*h/3.0;
}



int main()
{
  FILE * in = fopen("points.dat", "r");
   
  size_t n;
  fscanf(in, "%zd", &n);


  double * x = malloc(n * sizeof(double));
  double * y = malloc(n * sizeof(double));

    
  for (size_t i = 0; i < n; ++i) {
    fscanf(in,"%lf %lf", &x[i], &y[i]);
  }

  fclose(in);
  
  double s1 = simpson38(x, y);

  double s2 = simpson13(n-3, &y[3], x[3], x[n-1]);

  double integr = s1+s2;
  
  double correct = fint(x[n-1]) - fint(x[0]);

  printf("%.15g %.15g %.15g\n", integr, correct, integr-correct);

  free(x);
  free(y);
  
}
