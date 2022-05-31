#include <cmath>
#include <iostream>

/*
  Bisection method to locate roots of 

  f(x)=x^3+4x^2-10 in [1,2].
*/


double f(double x)
{
    return std::pow(x,3) + 4.0*std::pow(x,2) - 10.0;
}

bool samesign(double a, double b)
{
    // return a*b > 0.0;
    return std::signbit(a) == std::signbit(b);
}

int main()
{
    double constexpr toler{1e-8};

    //  initial limits
    double a{1.0};
    double b{2.0};

    double fa{f(a)};
    double fb{f(b)};

    if (samesign(fa, fb)) {
	return -1;
    }

    std::cout.precision(12);
    
    while (true) {
	double x {(a+b)/2.0};
	double fnval{f(x)};

	std::cout << "The current approximation is " << x
		  << "\nThe current value of F(x) is " << fnval
		  <<'\n';

	//  Check if root is found
	if (std::abs(fnval) < toler) {
	    break;
	}

	if (samesign(fnval, fa)) {
	    a = x;
	    fa = fnval;
	} else {
	    b = x;
	    fb = fnval;
	}	
    }
	
}
