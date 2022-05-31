#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
     Epilisi tis f(x)=0 me f(x)=sin(x)-x**2
     me th methodo Newton-Raphson.
*/


double f(double x);
double df(double x);


int
main()
{    
    double const toler = 1e-8; //     A small constant number

    double x = 1.0;     //    Initial approximation 

    while (fabs(f(x)) > toler) {
	x = x - f(x) / df(x);
    } 

    printf("%s%g\n","The root is approximately ", x);
    printf("%s%g\n", "The function value is ", f(x));

    return 0;
    
}


//    The function
double
f(double x)
{
    return sin(x)-x*x;
}



//     The first derivative 
double
df(double x)
{
    return cos(x) - 2.0 * x;
}
