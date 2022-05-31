#include <iostream>
#include <cmath>
#include <functional>


double correct(double x)
{
    return 1.0 - std::exp(-x) + x*x -x;    
}

double f(double x, double y)
{
    return x*x + x-y;
}


double heun(double x0, double y0, double x1,
	    std::function<double (double,double)> f) 
{
    auto h = x1-x0;
    
    auto k1 = h * f(x0,y0);
    auto k2 = h * f(x1, y0 + k1);

    return y0 + (k1 + k2)/2.0;
}

double ralston(double x0, double y0, double x1,
	       std::function<double (double,double)> f) 
{
    auto h = x1-x0;
    
    auto k1 = h * f(x0,y0);
    auto k2 = h * f(x0 + 2.0/3.0 * h, y0 + 2.0 / 3.0 * k1);

    return y0 + (k1 + 3.0 * k2) / 4.0;
}

int main()
{
    double constexpr a{0.0};
    double constexpr b{0.6};
    double h{0.01};
    double constexpr yinit{0.0};

    double x{a};
    double yh{yinit};
    double yr{yinit};

    while (x < b) {
	if (x+h > b) {
	    h = b-x;
	}

        yh = heun(x, yh, x+h, f);
	yr = ralston(x, yr, x+h, f);
	
	x += h;
	
	std::cout << x << ' ' << yh << ' ' << yr
		  << ' ' << correct(x) << '\n';
    }    
}
