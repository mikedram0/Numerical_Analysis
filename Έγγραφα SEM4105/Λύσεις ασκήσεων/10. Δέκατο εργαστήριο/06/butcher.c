#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double runge_kutta(size_t s,
		   double const A[], double const b[], double const c[],
		   double x0, double y0, double x1,
		   double (*f)(double x, double y))
{
    double const h = x1 - x0;

    double * k = malloc(s*sizeof(double));

    for (size_t i = 0; i < s; ++i) {
	double sum = y0;
	for (size_t j = 0; j < i; ++j) {
	    sum += A[i*s+j] * k[j];
	}
	k[i] = h * f(x0+c[i] * h, sum);
    }
    
    double y1 = y0;
    for (size_t i = 0; i < s; ++i) {
	y1 += b[i] * k[i];
    }

    free(k);
    
    return y1;
}


double correct(double x)
{
    return x + sqrt(1+2.0*x*x);
}

double f(double x, double y)
{
    return (y+x)/(y-x);
}



int main()
{
    double a = 0.0;
    double b = 2.0;
    double h = 0.1;
    double yinit = 1.0;

    double x = a;
    double y = yinit;

    size_t s = 4;

    double rk_A[] = {  0.0, 0.0, 0.0, 0.0,
	    0.5, 0.0, 0.0, 0.0,
	    0.0, 0.5, 0.0, 0.0,
	    0.0, 0.0, 1.0, 0.0 	    };
    
    double rk_b[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

    double rk_c[] = {0.0, 0.5, 0.5, 1.0};

    while (x < b) {
	if (x+h > b) {
	    h = b-x;
	}

        y = runge_kutta(s, rk_A, rk_b, rk_c, x, y, x+h, f);
	
	x += h;
	
	printf("%.5lf %.15lf %.15lf\n", x,y, correct(x));
    }

    return 0;
}


