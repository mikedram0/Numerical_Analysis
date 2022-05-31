#include <math.h>
#include <stdio.h>
#include <complex.h>


double  g(double x, double y)
{
  return 4.0 * x*x - pow(y,3) + 28.0;
}


double h(double x, double y)
{
  return 3.0 * pow(x,3) + 4.0 * y*y - 145.0;
}


double complex f(double complex z)
{
  double x = creal(z);
  double y = cimag(z);

  double re = g(x,y);
  double im = h(x,y);

  return CMPLX(re,im);
}


double complex secant(double complex x1, double complex x2,
		      double toler,
		      double complex (*f)(double complex x))
{
  double complex f1 = f(x1);
  double complex f2 = f(x2);
  double complex x;
    
  do {
    x = (x2 * f1 - x1 * f2) / (f1 - f2);

    x1 = x2;
    f1 = f2;

    x2 = x;
    f2 = f(x);

  } while (cabs(f2) > toler);

  return x;
}

int main(void)
{
  double complex z1 = CMPLX(0.0,0.0);
  double complex z2 = CMPLX(1.0,1.0);

  double complex z = secant(z1, z2, 1e-8, f);

  printf("%g %g\n", creal(z) , cimag(z));

  return 0;
}
