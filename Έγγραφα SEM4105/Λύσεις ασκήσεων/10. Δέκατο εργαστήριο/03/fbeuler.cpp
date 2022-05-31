#include <cmath>
#include <iostream>
#include <functional>


double f(double x, double y)
{
    return std::cos(x) - std::sin(y) + x*x;
}



double forward_euler(double x0, double y0, double x1,
		     std::function<double (double, double)> f)
{
    return y0 + (x1-x0) * f(x0,y0);
}


// y1 = y0 + (x1-x0) f(x1, y1)
double g(double y1, double x0, double y0, double x1,
	 std::function<double (double, double)> f)
{
    return y1 - y0 - (x1-x0) * f(x1,y1);
}


 
double backward_euler(double x0, double y0, double x1,
		      std::function<double (double, double)> f)
{
    // secant for g(y, x0,y0,x1,f) = 0

    // need two points: one is y1 from forward_euler
    double y1 = forward_euler(x0, y0, x1, f);
    //the other is y2, close to y1
    double y2 = y1 * 1.1;
    
    double y; // this will be the approximation to the unknown (y1)
    while (true) {
	double g1 = g(y1, x0, y0, x1, f);
	double g2 = g(y2, x0, y0, x1, f);
	
	y = y2 -  g2 * (y2-y1) / (g2 - g1);
	
	if (std::abs(g(y, x0, y0, x1, f)) < 1e-8) {
	    break;
	}

	y1 = y2;
	y2 = y;
    }

    return y;    
}



int main()
{
    double constexpr a{-1.0};
    double constexpr yinit{3.0};

    double constexpr b{1.0};
   
    double  h = 0.01;

    double x0{a};
    double y0f{yinit};
    double y0b{yinit};

    std::size_t const n = std::lround((b-a)/h);

    std::cout.precision(8);
    std::cout << "#X\tforward\t\tbackward\n";
    std::cout << x0 << '\t' << y0f << '\t' << y0b << '\n';
    
    for (std::size_t i{0}; i <=n; ++i) {
	// change h if x0 is beyond b:
	if (x0+h>b) {
	    h = b-x0;
	}
    
	y0f = forward_euler (x0, y0f, x0+h, f);
	y0b = backward_euler(x0, y0b, x0+h, f);
	x0 += h;

	std::cout << x0 << '\t' << y0f << '\t' << y0b << '\n';
    }
}
