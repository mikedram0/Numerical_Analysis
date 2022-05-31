#include <stdio.h>
#include <math.h>
#include <stddef.h>


double g(double x)
{
    return sin(x);
}


double simpson(double (*f)(double x), double a, double b, size_t n)
{
    if (n%2 != 0) {
      fprintf(stderr, "error\n");
    }
    
    double const h = (b-a)/n;

    double s = f(a)+f(b);

    double s1 = 0.0;
    for (size_t i = 1; i < n; i+=2) {
	s1 += f(a+i*h);
    }
    s += 4.0 * s1;

    double s2 = 0.0;
    for (size_t i = 2; i < n; i+=2) {
	s2 += f(a+i*h);
    }
    s += 2.0 * s2;
    
    return s*h/3.0;
}


int main()
{
    double const pi = 3.14159265358979323846;
    
    size_t n = 2;
    double const correct = 2.0;
    for (size_t k = 1; k < 10; ++k) {
        double integr = simpson(g, 0.0, pi, n);

	printf("%zd %g %g\n",n, integr, correct - integr);
	n*=2;
    }

    return 0;
}


