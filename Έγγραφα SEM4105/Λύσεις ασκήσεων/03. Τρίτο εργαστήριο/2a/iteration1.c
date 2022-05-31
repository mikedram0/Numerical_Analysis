#include <stdio.h>
#include <math.h>

double f(double x)
{
  return x*x - 6.0 * x + 5.0;
}



/*  
    if f(x) above then 
    f(x) = 0 => x = g(x) 
    with g as follows:
*/
double g(double x)
{
  return (x*x + 5.0) / 6.0;
}




/*
  Fixed point method to locate a root of 
   
  f(x)=x^2-6x+5
*/


int
main(void)
{
  double const toler = 1E-8;   //  a small constant number
  double const vfar = 1e5;  // a very large number

  double x = 1.3;  // initial guess for x
  while (1) {
    double fnval = f(x); // function value 
	
    if (fabs(x) > vfar) {
      printf("X is very large; probably g(x) is not appropriate\n");
      break;
    }

    printf("The current approximation is %.12g\n", x);
    printf("The current value of F(x) is %g\n\n", fnval);

    if (fabs(fnval) < toler) {
      break;
    }

    x = g(x);
  }

  return 0;
}
