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

double myexp(double x)
{
    return std::exp(x);
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


//     Y = a0 X^A1 => LN(Y) = LN(a0) + A1 * LN(X)
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


//     Y = a0 + a1 * LOG(X)
void lsqlog(vec const & x, vec const & y,
	    double & a1, double & a0, double & r2)
{
    if (!allpositive(x)) {
	std::cerr << "Data not appropriate\n";
	return; 
    }

    auto const n = x.size();
    vec zx(n);

    std::transform(x.cbegin(), x.cend(), zx.begin(), mylog);

    lsqlin(zx, y, a1, a0, r2);
  
}


//     Y = a0 + a1 * EXP(X)
void lsqexp(vec const & x, vec const & y,
	    double & a1, double & a0, double & r2)
{
    auto const n = x.size();
    vec zx(n);
        
    std::transform(x.cbegin(), x.cend(), zx.begin(), myexp);
    
    lsqlin(zx, y, a1, a0, r2);
}


int main()
{
    vec x, y;
    read_data("points.dat", x, y);

    std::cout <<  u8"Μέθοδος:\n"
	      << "1) y = a * x + b\n"
	      << "2) y = a * x^b\n"
	      << "3) y = a + b * exp(x)\n"
	      << "4) y = a + b * ln(x)\n"
	      << u8"Δώσε μέθοδο.\n";

    int method;
    std::cin >> method;

    double a1, a0, r2;

    switch (method) {
    case 1:
	lsqlin(x, y, a1, a0, r2);
	std::cout << "y = " << a1 << " x + " << a0 << '\n';
	break;
    case 2:
	lsqpow(x, y, a1, a0, r2);
	std::cout << "y = " << a0 << " * x^" << a1 << '\n';
	break;
    case 3:
	lsqexp(x, y, a1, a0, r2);
	std::cout << "y = " << a0 << " + " << a1 <<  " * exp(x)" <<'\n';
	break;
    case 4:
	lsqlog(x, y, a1, a0, r2);
	std::cout << "y = " << a0 << " + " << a1 << " * ln(x)" << '\n';
	break;
    }

    std::cout << "r^2 = " << r2 << '\n';
}


