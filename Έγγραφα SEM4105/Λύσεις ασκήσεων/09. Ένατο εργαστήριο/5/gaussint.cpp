#include <iostream>
#include <cmath>
#include <functional>
#include <array>


double f(double x)
{
    return std::pow(x,3) * std::exp(-x);
}


/*     D FINT(X) / D X = F */
double fint(double x)
{
    return -(6.0 + 6.0 * x + 3.0 * x*x + x*x*x) * std::exp(-x);
}


double gauss2(double a, double b, std::function<double(double)> f)
{
    std::array<double,2> w{1.0, 1.0};
    std::array<double,2> x{-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
  
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
	double y = ( (b-a) * x[i] + (b+a) ) / 2.0;
	sum += w[i] * f(y);
    }

    return sum * (b-a) / 2.0;
}

double gauss3(double a, double b, std::function<double(double)> f)
{
std::array<double,3> w{5.0/9.0, 8.0/9.0, 5.0/9.0};
std::array<double,3> x{-std::sqrt(0.6), 0.0, std::sqrt(0.6)};
  
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
	double y = ( (b-a) * x[i] + (b+a) ) / 2.0;
	sum += w[i] * f(y);
    }

    return sum * (b-a) / 2.0;
}


/*   gauss integration  */
int main()
{
    double const a = 2.1;
    double const b = 5.2;
    
    std::cout << "Gauss 2 points: " <<  gauss2(a, b, f) << '\n';
    std::cout << "Gauss 3 points: " <<  gauss3(a, b, f) << '\n';
    std::cout << "correct: " << fint(b) - fint(a) << '\n';
}


