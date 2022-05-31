#include <cmath>
#include <iostream>
#include <functional>

double secant(double x1, double x2, double toler,
	      std::function<double (double)> func)
{
    double f1{func(x1)};
    double f2{func(x2)};
    double x;
    
    do {
	x = (x2 * f1 - x1 * f2) / (f1 - f2);

	x1 = x2;
	f1 = f2;

	x2 = x;
	f2 = func(x);

    } while (std::abs(f2) > toler);

    return x;
}


double f(double x)
{
    return 3.0 * std::log(x) + 5.0;
}


    // ///   SECANT   ///////
int main()
{
    double constexpr toler{1e-8}; // A small constant number

//     Initial approximations
    double a{0.1};
    double b{0.2};
    
    double x = secant(a,b,toler, f);

    std::cout << u8"Μια ρίζα είναι περίπου " << x << '\n';
    std::cout << u8"με τιμή συνάρτησης "  << f(x) << '\n';
}
