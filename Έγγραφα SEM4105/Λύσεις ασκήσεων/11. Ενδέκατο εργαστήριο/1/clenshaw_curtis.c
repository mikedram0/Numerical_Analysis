#include <stdio.h>
#include <math.h>
#include <stdlib.h>



double const  pi = 3.14159265358979323846;


double c(int i, int n)
{
    return (i%n == 0) ? 1.0 : 2.0;
}

double b(int j, int n)
{
    return (j%(n/2) == 0) ? 1.0 : 2.0;
}




void coefs(int n, double w[])
{
  int i, j;

  --n;

  for (i=0; i<=n; ++i) {
    w[i] = 0.0;
    for (j=0; j <= n/2; ++j) {
      w[i] += b(j,n) / (1-4.0*j*j) * cos(2.0*i*j*pi/n);
    }
    
    w[i] *= c(i,n)/n;
  }
}


/*
  \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x
*/
double f(double x)
{
  return 2.0/ (1.0+4.0*x*x);
}



double calc(int n)
{
  int i;
  double s;
  double * w = malloc((n+1) * sizeof(double));
  
  coefs(n+1,w);

  s = 0.0;
  for (i=0; i < n+1; ++i) {
    double x = cos(i*pi/n);
    s += w[i] * f(x);
  }
    
  return s;
}



int main()
{
    double const correct = 2.0*atan(2.0);

    int n = 3;
    while (1) {
	double s = calc(n);
	double diff = s-correct;
	printf("%g %g %g\n", s, correct, diff);
	if (fabs(diff) < 1e-12) {
	  break;	    
	}
	++n;	
    }
    printf("%d\n", n);

    return 0;
}


