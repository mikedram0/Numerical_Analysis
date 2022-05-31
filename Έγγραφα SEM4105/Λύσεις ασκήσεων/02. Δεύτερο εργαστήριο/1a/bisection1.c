#include <stdio.h>
#include <math.h>


/*
  Bisection method to locate roots of 

  f(x)=x^3+4x^2-10 in [1,2].
*/


double f(double x)
{
  return pow(x,3) + 4.0*pow(x,2) - 10.0;
}

int samesign(double a, double b)
{
    // return a*b > 0.0;
  return copysign(a,b) == a;
}

int main(void)
{
  const double toler = 1e-8;

  //  initial limits
  double a = 1.0;
  double b = 2.0;

  double fa = f(a);
  double fb = f(b);

  if (samesign(fa,fb)) {
    return -1;
  }

  while (1) {
    double x = (a+b)/2.0;
    double fnval = f(x);

    printf("%s%.12lf\n%s%le\n", "The current approximation is ", x,
	   "The current value of F(x) is ", fnval);

    //  Check if root is found
    if (fabs(fnval) < toler) {
      break;
    }

    if (samesign(fnval, fa) > 0.0) {
      a = x;
      fa = fnval;
    } else {
      b = x;
      fb = fnval;
    }	
  }

  return 0;
}
