#include <fstream>
#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>



using vec = std::vector<double>;



double fint(double x)
{
    return (std::sin(x) - std::cos(x)) / 2.0 * std::exp(x);
}



double simpson38(vec const & x, vec const & y)
{
    double const h = x[1]-x[0];
    
    return 3.0*h/8.0 * (y[0]+3.0*(y[1] +y[2]) + y[3]);
}



double simpson13(vec const & y, double a, double b)
{
    auto n = y.size()-1;
    
    if (n%2 != 0) {
	std::cerr << "error\n";
    }
    
    double const h{(b-a)/n};

    double s{y[0] + y[n]};

    double s1{0.0};
    for (std::size_t i{1}; i < n; i+=2) {
	s1 += y[i];
    }

    s += 4.0 * s1;

    double s2{0.0};
    for (std::size_t i{2}; i < n; i+=2) {
	s2 += y[i];
    }

    s += 2.0 * s2;
    
    return s*h/3.0;
}





int main()
{
    std::ifstream in{"points.dat"};

    std::size_t n;
    in >> n;

    vec x(n);
    vec y(n);

    for (vec::size_type i{0}; i < n; ++i) {
	in >> x[i] >> y[i];
    }

    double s1 = simpson38(x, y);

    for (std::size_t i{0}; i <3; ++i){
	y.erase(y.begin());
    }

    double s2 = simpson13(y,x[3],x[n-1]);

    double integr = s1+s2;
    
    double correct = fint(x[n-1]) - fint(x[0]);

    std::cout.precision(15);
    
    std::cout  << integr << ' ' << correct
	       << ' ' << integr-correct<<'\n';
}
