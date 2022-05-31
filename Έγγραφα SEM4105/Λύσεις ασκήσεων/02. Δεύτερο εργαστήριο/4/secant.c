#include <math.h>
#include <stdio.h>

double secant(double x1, double x2, double toler, double (*func)(double x))
{
    double f1 = func(x1);
    double f2 = func(x2);
    double x;
    
    do {
	x = x2 - f2 * (x2 - x1) / (f2 - f1);

	x1 = x2;
	f1 = f2;

	x2 = x;
	f2 = func(x);

    } while (fabs(f2) > toler);

    return x;
}


double f(double x)
{
    return 3.0 * log(x) + 5.0;
}


    // ///   SECANT   ///////
int main(void)
{
  double const toler = 1e-8; // A small constant number

//     Initial approximations
  double a = 0.1;
  double b = 0.2;
  
  double x = secant(a,b,toler, f);
  
  printf("A root is approximately %lf\n", x);
  printf("The function value is %le\n", f(x));

  return 0;
}
