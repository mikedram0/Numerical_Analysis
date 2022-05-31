#include <iostream>
#include <cmath>

/*
     Epilisi tis f(x)=0 me f(x)=sin(x)-x**2
     me th methodo Newton-Raphson.
*/


double f(double x);
double df(double x);


int
main()
{    
    double constexpr toler{1e-8}; //     A small constant number

    double x{1.0};     //    Initial approximation 

    while (std::abs(f(x)) > toler) {
	x = x - f(x) / df(x);
    } 

    std::cout <<"The root is approximately " << x << '\n';
    std::cout <<"The function value is " << f(x) << '\n';

}


//    The function
double
f(double x)
{
    return std::sin(x)-x*x;
}



//     The first derivative 
double
df(double x)
{
    return std::cos(x) - 2.0 * x;
}
