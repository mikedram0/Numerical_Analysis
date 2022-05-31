#include <stdio.h>
#include <math.h>
#include <complex.h>


// Find the roots of f(x)=0 with f(x)=sin(x) - x^2
// using the Muller's method.
// with complex arithmetic

double complex f(double complex x);

    
int main(void)
{
    double const epsilon = 1e-6;
    double const delta = 1e-6;

    double complex x0 = 0.1; 
    double complex x1 = 0.2;
    double complex x2 = 0.5;

    double complex f0 = f(x0);
    double complex f1 = f(x1);
    double complex f2 = f(x2);

    double complex x,t;
    
    while (1) {
	double complex w1 = (f2-f1)/(x2-x1);
	double complex w0 = (f2-f0)/(x2-x0);

	double complex a = (w1-w0)/(x1-x0);
	double complex b = w0+a*(x2-x0);
	double complex c = f2;

	double complex p = csqrt(b*b-4.0*a*c);
	
	double complex d1 = b + p;
	double complex d2 = b - p;

	double complex d = (cabs(d1) > cabs(d2) ? d1 : d2);

	x = x2 - 2.0*c / d;

	t = f(x);
	    
	//        termination condition
	if ((cabs(x - x2) < delta) || (cabs(t) < epsilon)) {
	    break;
	}

	x0 = x1;
	x1 = x2;
	x2 = x;

	f0 = f1;
	f1 = f2;
	f2 = t;
    }

    printf("%s%g\t%g\n", "H riza einai ", creal(x), cimag(x));
    t = f(x);
    printf("%s%g\t%g\n", "H timi ths synarthshs einai ", creal(t), cimag(t));

    return 0;
}


//     define f(x)
double complex f(double complex x)
{
    return csin(x) - x*x;
}
