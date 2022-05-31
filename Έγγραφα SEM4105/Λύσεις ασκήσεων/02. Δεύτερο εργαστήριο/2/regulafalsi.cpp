#include <iostream>
#include <cmath>
#include <functional>

double f(double x)
{
    return -2.0 + 6.2 * x - 4.0 * x*x + 0.7 * std::pow(x,3);
}

double g(double x)
{
    return std::pow(x,10) - 0.95;
}


bool samesign(double a, double b)
{
    // return a*b > 0.0;
    return std::signbit(a) == std::signbit(b);
}
    

int bisection(double a, double b, //     Initial limits
	      std::function<double (double)> f, // function
	      double toler, // tolerance
	      double & x, // final approximation
	      int & iter // Number of iterations needed
    )
{
    double fa{f(a)};
    double fb{f(b)};

    if (samesign(fa,fb)) {
	return -1;
    }

    for (iter = 1; ; ++iter) {
	x = (a+b)/2.0;
	double fnval{f(x)};

	//  Check if root is found
	if (std::abs(fnval) < toler) {
	    break;
	}

	if (samesign(fnval,fa)) {
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
	  std::function<double (double)> f, // function
	  double toler, // tolerance
	  double & x, // final approximation
	  int & iter // Number of iterations needed
    )
{
    double fa{f(a)}, fb{f(b)};

    if (samesign(fa,fb)) {
	std::cerr << "limits are wrong\n";
	return -1;
    }
    
    for (iter = 1; ; ++iter) {
  	x = a - fa * (b-a) / (fb-fa);
	double fx = f(x);

	//     check if root is found
	if (std::abs(fx) < toler) {
	    break;
	}

	//     new limits
	if (samesign(fa,fx)) {
	    a = x;
	    fa = fx;
	} else {
	    b = x;
	    fb = fx;
	}
    }

    return 0;
}




int main()
{
    double x;
    double constexpr toler{1e-6};
    int iter;


    int err = regfal(0.4, 0.6, f, toler, x, iter);
    if (err != 0) {
	std::cerr << "error in regfal\n";
	return -1;
    }
    std::cout << "Question b: X is " << x << '\n';


    std::cout << "Question c: \n";

    err = regfal(0.0, 1.4, g, toler, x, iter);
    if (err != 0) {
	std::cerr << "error in regfal\n";
	return -1;
    }
    std::cout << "x is " << x << " by false position method,"
	      << " in " << iter << " iterations\n";

    err = bisection(0.0, 1.4, g, toler, x, iter);
    if (err != 0) {
	std::cerr << "error in bisection\n";
	return -1;
    }

    std::cout << "x is " << x << " by bisection method,"
	      << " in " << iter << " iterations\n";

}


