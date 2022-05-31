#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>


using vec = std::vector<double>;

void create(std::string fname)
{
    std::size_t constexpr n{15};
    double constexpr a{2.0}, b{4.0}, step{(b-a)/(n-1)};
    
    std::ofstream out{fname};

    out << n << '\n';
    for (std::size_t i{0}; i < n; ++i) {
	double x{a + step * i};
	double y{std::sin(x)};
	out << x << ' ' << y << '\n';
    }
}



//     Read data file
void read_data(std::string fname, vec & x, vec & y)
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


double q(std::size_t i, vec const & x, double z)
{
    double r{1.0};
    for (std::size_t j{0}; j < i; ++j) {
	r *= z-x[j];
    }

    return r;
}


void coefficients(vec const & x, vec const & y, vec & a)
{
    for (std::size_t i{0}; i < y.size(); ++i) {
	double s{y[i]};
	for (std::size_t j{0}; j < i; ++j) {
	    s -= a[j] * q(j,x,x[i]);
	}
	a[i] = s / q(i, x, x[i]);
    }
}


double poly(vec const & a, vec const & x, double z)
{
    double r{0.0};
    for (std::size_t i{0}; i < a.size(); ++i) {
	r += a[i] * q(i,x,z);
    }
    return r;
}


int main()
{
    create("points.dat");
    
    vec x, y;
    read_data("points.dat", x,y);

    std::size_t constexpr m{100};

    auto const p = std::minmax_element(x.cbegin(), x.cend());

    double const & min = *(p.first);
    double const & max = *(p.second);
    double const step = (max-min)/(m-1);

    vec a(x.size());
    coefficients(x,y,a);
    
    for (std::size_t i{0}; i < m; ++i) {
	double xout = min + step * i;
	double yout = poly(a,x, xout);

	std::cout << xout << ' ' << yout << '\n';
    }
}
