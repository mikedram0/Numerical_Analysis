#include <iostream>
#include <cmath>
#include <functional>


double f(double x, double y)
{
    return std::cos(x)-x*std::sin(x);
}


double forward_euler(double x0, double y0, double x1,
		     std::function<double (double, double)> f)
{
    return y0 + (x1-x0) * f(x0, y0);    
}
		     
double correct(double x)
{
    return 2.0 + x*std::cos(x);
}

int main()
{
    double constexpr a{0.0};
    double constexpr yinit{2.0};

    double constexpr b{3.0};
   

    double h = 0.01;

    double x0{a};
    double y0{yinit};
    
    std::cout << x0 << ' ' << y0 << ' ' << correct(x0) << '\n';
    while (x0 < b) {
	if (x0 +h > b) {
	    h = b-x0;
	}
	y0 = forward_euler(x0, y0, x0+h, f);
	x0 += h;

	std::cout << x0 << ' ' << y0 << ' ' << correct(x0) << '\n';
    }    
}
