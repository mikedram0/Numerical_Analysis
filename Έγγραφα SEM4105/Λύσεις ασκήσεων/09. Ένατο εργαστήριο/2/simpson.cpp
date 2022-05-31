#include <iostream>
#include <cmath>
#include <functional>
#include <cstddef>


double g(double x)
{
    return std::sin(x);
}


double simpson(std::function<double (double)> f,
	       double a, double b, std::size_t n)
{
    if (n%2 != 0) {
	std::cerr << "error\n";
    }
    
    double const h{(b-a)/n};

    double s{f(a)+f(b)};

    double s1{0.0};
    for (std::size_t i{1}; i < n; i+=2) {
	s1 += f(a+i*h);
    }

    s += 4.0 * s1;

    double s2{0.0};
    for (std::size_t i{2}; i < n; i+=2) {
	s2 += f(a+i*h);
    }

    s += 2.0 * s2;
    
    return s*h/3.0;
}


int main()
{
    double constexpr pi{3.14159265358979323846};
    
    std::size_t n{2};
    double constexpr correct = 2.0;
    for (std::size_t k{1}; k < 10; ++k) {
        double integr = simpson(g, 0.0, pi, n);

	std::cout << n << ' ' << integr << ' '
		  << correct - integr << '\n';
	n*=2;
    }
    
}


