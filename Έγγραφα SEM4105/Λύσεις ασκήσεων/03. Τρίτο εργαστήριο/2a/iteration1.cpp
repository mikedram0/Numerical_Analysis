#include <iostream>
#include <cmath>

double f(double x)
{
    return x*x - 6.0 * x + 5.0;
}



/*  
    if f(x) above then 
    f(x) = 0 => x = g(x) 
    with g as follows:
*/
double g(double x)
{
    return (x*x + 5.0) / 6.0;
}




/*
  Fixed point method to locate a root of 
   
   f(x)=x^2-6x+5
*/


int main()
{
    double constexpr toler{1E-8};   //  a small constant number
    double constexpr vfar{1e5};  // a very large number

    double x;
    std::cin >> x;

    while (true) {
	double fnval{f(x)}; // function value 
	
	if (std::abs(x) > vfar) {
	    std::cout<< "X is very large; probably g(x) is not appropriate\n";
	    break;
	}

	std::cout.precision(12);
	std::cout << "The current approximation is " << x <<'\n';
	std::cout << "The current value of F(x) is " << fnval << "\n\n";

	if (std::abs(fnval) < toler) {
	    break;
	}

	x = g(x);
    }
}
