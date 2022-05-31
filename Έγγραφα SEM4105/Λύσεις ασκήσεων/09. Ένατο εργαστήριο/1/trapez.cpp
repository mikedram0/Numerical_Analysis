#include <iostream>
#include <cmath>
#include <cstddef>
#include <functional>


double g(double x)
{
    return std::sin(x);
}


double trapez(std::function<double(double)> f,
	      double a, double b, std::size_t n)
{
    double const h{(b-a)/n};

    double s{(f(a)+f(b))/2.0};

    for (std::size_t i{1}; i < n; ++i) {
	s += f(a+i*h);
    }

    return s*h;
}


int main()
{
    double constexpr pi{3.14159265358979323846};

    std::size_t n{2};

    double constexpr correct = 2.0;
    for (std::size_t k{1}; k < 10; ++k) {
        double integr = trapez(g, 0.0, pi, n);

	std::cout << n << ' ' << integr << ' '
		  << correct - integr << '\n';
	n*=2;
    }
    
}


