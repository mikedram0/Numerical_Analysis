#include <stdio.h>
#include <math.h>


double f(double x, double y)
{
  return cos(x)-x*sin(x);
}


double forward_euler(double x0, double y0, double x1,
		     double (*f)(double x, double y))
{
  return y0 + (x1-x0) * f(x0, y0);    
}
		     

double correct(double x)
{
    return 2.0 + x*cos(x);
}


int main()
{
  double const a = 0.0;
  double const yinit = 2.0;

  double const b = 3.0;
   

  double h = 0.01;

  double x0 = a;
  double y0 = yinit;

  printf("%g %g %g\n", x0, y0, correct(x0));

  while (x0 < b) {
    if (x0 +h > b) {
      h = b-x0;
    }
    
    y0 = forward_euler(x0, y0, x0+h, f);
    x0 += h;
    printf("%g %g %g\n", x0, y0, correct(x0));
  }
    
  return 0;
}
