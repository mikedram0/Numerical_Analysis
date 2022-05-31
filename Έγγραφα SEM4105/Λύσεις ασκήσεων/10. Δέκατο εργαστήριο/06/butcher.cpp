#include <vector>
#include <cmath>
#include <iostream>
#include <functional>
#include <numeric>


using vec = std::vector<double>;


double runge_kutta(vec const& A, vec const & b, vec const & c,
		   double x0, double y0, double x1,
		   std::function<double (double,double)> f)
{
    auto const s = c.size();
   
    double const h{x1 - x0};

    vec k(s);

    for (std::size_t i = 0; i < s; ++i) {
	double sum{y0};
	for (std::size_t j = 0; j < i; ++j) {
	    sum += A[i*s+j] * k[j];
	}
	k[i] = h * f(x0+c[i] * h, sum);
    }

    return std::inner_product(b.cbegin(), b.cend(), k.cbegin(), y0);
}


double correct(double x)
{
    return x + std::sqrt(1+2.0*x*x);
}

double f(double x, double y)
{
    return (y+x)/(y-x);
}



int main()
{
    double constexpr a{0.0};
    double constexpr b{2.0};
    double h{0.1};
    double constexpr yinit{1.0};

    double x{a};
    double y{yinit};


    vec rk_A{
	0.0, 0.0, 0.0, 0.0,
	0.5, 0.0, 0.0, 0.0,
	0.0, 0.5, 0.0, 0.0,
	0.0, 0.0, 1.0, 0.0
    };
    
    vec rk_b{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

    vec rk_c{0.0, 0.5, 0.5, 1.0};

    std::cout.precision(15);
    
    while (x < b) {
	if (x+h > b) {
	    h = b-x;
	}

        y = runge_kutta(rk_A, rk_b, rk_c, x, y, x+h, f);
	
	x += h;
	
	std::cout << x << ' ' << y << ' ' << correct(x) << '\n';
    }    
}


