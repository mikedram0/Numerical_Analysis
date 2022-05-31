#include <vector>
#include <iostream>
#include <numeric>
#include <fstream>
#include <functional>
#include <cmath>
#include <algorithm>
#include <string>


using vec = std::vector<double>;

double mylog(double x)
{
    return std::log(x);
}

bool allpositive(vec const & x)
{
    auto gt0 = [] (auto & x) -> bool { return x > 0; };
	
    return std::all_of(x.cbegin(), x.cend(), gt0);
}


//     Read data file
void read_data(std::string const & fname, vec & x, vec & y)
{
    std::ifstream in{fname};

    int n;
    in >> n;
    
    x.reserve(n);
    y.reserve(n);

    for (int i{0}; i < n; ++i) {
	double a, b;
	in >> a >> b;
	x.push_back(a);
	y.push_back(b);
    }
}


//     Y = A1 * X + a0
void lsqlin(vec const & x, vec const & y,
	    double & a1, double & a0, double & r2)
{
    auto const n = x.size();

    auto const sx = std::accumulate(x.cbegin(), x.cend(), 0.0);
    auto const sy = std::accumulate(y.cbegin(), y.cend(), 0.0);
    auto const sx2 = std::inner_product(x.cbegin(), x.cend(), x.cbegin(), 0.0);
    auto const sy2 = std::inner_product(y.cbegin(), y.cend(), y.cbegin(), 0.0);
    auto const sxy = std::inner_product(x.cbegin(), x.cend(), y.cbegin(), 0.0);

    auto const d = n * sx2 - sx * sx;

    a1 = (n * sxy - sx * sy) / d; 
    a0 = (sy - a1 * sx) / n;

    r2 = std::pow(n * sxy - sx * sy,2) / d / (n * sy2 - sy*sy);
}




//     Y = a0*X^A1 => LN(Y) = LN(a0) + A1 * LN(X)
void lsqpow(vec const & x, vec const & y,
	    double & a1, double & a0, double & r2)
{
    if (!allpositive(x) || !allpositive(y)) {
	std::cerr << "Data not appropriate\n";
	return; 
    }

    auto const n = x.size();
    vec zx(n);
    vec zy(n);
    double za0;

    std::transform(x.cbegin(), x.cend(), zx.begin(), mylog);
    std::transform(y.cbegin(), y.cend(), zy.begin(), mylog);    

    lsqlin(zx, zy, a1, za0, r2);
  
    a0 = std::exp(za0);
}


int main()
{
    vec x, y;
    read_data("stefbol.dat", x,y);

    double a1, a0, r2;
    lsqpow(x, y, a1, a0, r2);
    
    double constexpr embado{0.05e-4};

    std::cout << u8"Σταθερά Stefan-Boltzmann: "
	      << a0 / embado << " J/K^4/m^2/s\n";
    std::cout << u8"Εκθέτης: " << a1 <<'\n';
    std::cout << "r^2 = " << r2 << '\n';
}

