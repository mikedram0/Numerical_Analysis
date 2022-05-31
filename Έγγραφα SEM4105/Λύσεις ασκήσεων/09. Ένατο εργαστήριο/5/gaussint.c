#include <stdio.h>
#include <math.h>


double f(double x)
{
  return pow(x,3) * exp(-x);
}


/*     D FINT(X) / D X = F */
double fint(double x)
{
  return -(6.0 + 6.0 * x + 3.0 * x*x + x*x*x) * exp(-x);
}


double gauss2(double a, double b, double (*f)(double x))
{
  double const w[2] = {1.0, 1.0};
  double const x[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
  
  double sum = 0.0;
  for (size_t i = 0; i < 2; ++i) {
    double y = ( (b-a) * x[i] + (b+a) ) / 2.0;
    sum += w[i] * f(y);
  }

  return sum * (b-a) / 2.0;
}

double gauss3(double a, double b, double (*f)(double x))
{
  double const w[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
  double const x[3] = {-sqrt(0.6), 0.0, sqrt(0.6)};
  
  double sum = 0.0;
  for (size_t i = 0; i < 3; ++i) {
    double y = ( (b-a) * x[i] + (b+a) ) / 2.0;
    sum += w[i] * f(y);
  }

  return sum * (b-a) / 2.0;
}


/*   gauss integration  */
int main(void)
{
  double const a = 2.1;
  double const b = 5.2;

  printf("Gauss 2 points: %f\n", gauss2(a, b, f));
  printf("Gauss 3 points: %f\n", gauss3(a, b, f));
  printf("correct: %f\n",  fint(b) - fint(a));

  return 0;
}


