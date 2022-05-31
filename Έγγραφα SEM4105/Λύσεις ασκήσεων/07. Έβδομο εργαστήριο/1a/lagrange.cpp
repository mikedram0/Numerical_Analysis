#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstddef>


using vec = std::vector<double>;

void create(const std::string & fname)
{
    constexpr std::size_t n{15};
    constexpr double a{2.0}, b{4.0}, step{(b-a)/(n-1)};
    
    std::ofstream out{fname};

    out << n << '\n';
    for (std::size_t i{0}; i < n; ++i) {
	double x{a + step * i};
	double y{std::sin(x)};
	out << x << ' ' << y << '\n';
    }
}



//     Read data file
void read_data(const std::string & fname, vec & x, vec & y)
{
    std::ifstream in{fname};

    std::size_t n;
    in >> n;
    
    x.reserve(n);
    y.reserve(n);

    for (std::size_t i{0}; i < n; ++i) {
	double a, b;
	in >> a >> b;
	x.push_back(a);
	y.push_back(b);
    }
}


double ell(std::size_t i, double z, const vec & x)
{
    double p{1.0};
    
    for(std::size_t j{0}; j < x.size(); ++j) {
	if (j==i) {
	    continue;
	}

	p *= (z-x[j]) / (x[i]-x[j]);
    }

    return p;
}


double poly(const vec & x, const vec & y, double z)
{
    double s{0.0};
    for (std::size_t i{0}; i < x.size(); ++i) {
	s += ell(i,z,x) * y[i];
    }

    return s;
}

int main()
{
    create("points.dat");

    vec x, y;
    read_data("points.dat", x,y);

    constexpr std::size_t m{100};

    const auto p = std::minmax_element(x.cbegin(), x.cend());

    const double & min = *p.first;
    const double & max = *p.second;

    const double step = (max-min)/(m-1);

    for (std::size_t i{0}; i < m; ++i) {
	double xout = min + step * i;
	double yout = poly(x,y, xout);

	std::cout << xout << ' ' << yout << '\n';
    }
}
