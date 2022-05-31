#include <stdio.h>
#include <math.h>


double f(double x)
{
  return -2.0 + 6.2 * x - 4.0 * x*x + 0.7 * pow(x,3);
}



double g(double x)
{
  return pow(x,10) - 0.95;
}


    

int
bisection(double a, double b, //     Initial limits
	  double (*f)(double x), // function
	  double toler, // tolerance
	  double * xfin, // final approximation
	  int * it // Number of iterations needed
	  )
{
  double fa = f(a);
  double fb = f(b);

  if (fa * fb > 0.0) {
    return -1;
  }

  for (int iter = 1; ; ++iter) {
    double x = (a+b)/2.0;
    double fnval = f(x);

    //  Check if root is found
    if (fabs(fnval) < toler) {
      *it = iter;
      *xfin = x;
      break;
    }

    if (fnval * fa > 0.0) {
      a = x;
      fa = fnval;
    } else {
      b = x;
      fb = fnval;
    }	
  }

  return 0;
}




//     Regula falsi (false position) method

int
regfal(double a, double b, //     Initial limits
       double (*f)(double x), // function
       double toler, // tolerance
       double * xfin, // final approximation
       int * it // Number of iterations needed
       )
{
  double fa = f(a), fb = f(b);

  if (fa * fb >= 0.0) {
    fprintf(stderr, "limits are wrong\n");
    return -1;
  }
    
  for (int iter = 1; ; ++iter) {
    double x = a - fa * (b-a) / (fb-fa);
    double fx = f(x);

    //     check if root is found
    if (fabs(fx) < toler) {
      *it = iter;
      *xfin = x;
      break;
    }

    //     new limits
    if (fa * fx > 0.0) {
      a = x;
      fa = fx;
    } else {
      b = x;
      fb = fx;
    }
  }

  return 0;
}




int main(void)
{
  double x;
  const double toler = 1e-6;
  int iter;


  int err = regfal(0.4, 0.6, f, toler, &x, &iter);
  if (err != 0) {
    fprintf(stderr, "error in regfal\n");
    return -1;
  }
  printf("Question b: X is %lf\n", x);


  printf("Question c: \n");

  err = regfal(0.0, 1.4, g, toler, &x, &iter);
  if (err != 0) {
    fprintf(stderr, "error in regfal\n");
    return -1;
  }
  printf("x is %lf by false position method, in %d iterations\n", x, iter);

  err = bisection(0.0, 1.4, g, toler, &x, &iter);
  if (err != 0) {
    fprintf(stderr, "error in bisection\n");
    return -1;
  }
  printf("x is %lf by bisection method, in %d iterations\n", x, iter);

  return 0;
}


