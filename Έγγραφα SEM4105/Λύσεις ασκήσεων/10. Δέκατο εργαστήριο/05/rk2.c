#include <stdio.h>
#include <math.h>


double correct(double x)
{
    return 1.0 - exp(-x) + x*x -x;    
}

double f(double x, double y)
{
    return x*x + x-y;
}


double heun(double x0, double y0, double x1,
	    double (*f)(double x, double y))
{
    double h = x1-x0;
    
    double k1 = h * f(x0,y0);
    double k2 = h * f(x1, y0 + k1);

    return y0 + (k1 + k2)/2.0;
}

double ralston(double x0, double y0, double x1,
	       double (*f)(double x, double y))
{
    double h = x1-x0;
    
    double k1 = h * f(x0,y0);
    double k2 = h * f(x0 + 2.0/3.0 * h, y0 + 2.0 / 3.0 * k1);

    return y0 + (k1 + 3.0 * k2) / 4.0;
}

int main()
{
    double  a = 0.0;
    double  b = 0.6;
    double  h = 0.01;
    double  yinit = 0.0;

    double x = a;
    double yh = yinit;
    double yr = yinit;

    printf("x\t\tHeun\t\tRalston\t\tcorrect\n");
    while (x < b) {
	if (x+h > b) {
	    h = b-x;
	}

        yh = heun(x, yh, x+h, f);
	yr = ralston(x, yr, x+h, f);
	
	x += h;

	printf("%.5lf %.15lf %.15lf %.15lf\n", x, yh, yr, correct(x));
    }    
}
