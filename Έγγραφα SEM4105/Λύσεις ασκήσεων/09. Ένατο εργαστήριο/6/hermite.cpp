#include <iostream>
#include <cmath>
#include <array>


double factorial(int n)
{
    return (n==0 ? 1.0 : factorial(n-1)*n);
}



double hermite(int n, double x)
{
    double r;
    
    switch (n) {
    case 0:
	r = 1.0;
	break;
    case 1:
	r = 2.0 * x;
	break;
    case 2:
	r = 4.0 * x * x - 2.0;
	break;
    case 3:
	r = x * (8.0 * x * x - 12.0);
	break;
    case 4:
	r = x * x * (16.0 * x * x - 48.0) + 12.0;
	break;
    default:
	std::cerr << "not implemented\n";
	break;
    }

    return r;
}
      

double weight(int n, double x)
{
    double constexpr pi{3.14159265358979323846};
    
    return std::pow(2.0,n-1) * factorial(n) * std::sqrt(pi)
	/ std::pow(n * hermite(n-1,x), 2);
  
}


double f(double x)
{
    return x*x;
}



int main()
{
    std::szie_t constexpr n{4};
    double constexpr pi{3.14159265358979323846};
    
    std::array<double,n> x;
    x[0] = std::sqrt((3.0-std::sqrt(6.0)) / 2.0);
    x[1] = std::sqrt((3.0+std::sqrt(6.0)) / 2.0);
    x[2] = -x[0];
    x[3] = -x[1];

    double s{0.0};
    for (std::array<double,n>::size_type i{0}; i < n; ++i) {
	double const w = weight(n, x[i]);
	double const fval = f(x[i]);
	s += w * fval;
    }

    std::cout.precision(15);

    std::cout << "the integral is " << s << '\n'
	      << "the correct value is "<< std::sqrt(pi) / 2.0 << '\n';
    std::cout << "exactly the same (why?)\n";

}
