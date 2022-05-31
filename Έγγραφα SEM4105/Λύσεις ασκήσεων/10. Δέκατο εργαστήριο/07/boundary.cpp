#include <iostream>
#include <functional>
#include <cstdlib>

double secant(double x1, double x2, double toler,
	      std::function<double (double)> func)
{
    double f1{func(x1)};
    double f2{func(x2)};
    double x;
    
    do {
	x = (x2 * f1 - x1 * f2) / (f1 - f2);

	x1 = x2;
	f1 = f2;

	x2 = x;
	f2 = func(x);	
    } while (std::abs(f2) > toler);

    return x;
}



double f(double x, double y)
{
    return (x-2*y)/(x+2*y);
}


double forward_euler(double x0, double y0, double x1,
		     std::function<double (double, double)> f)
{
    return y0 + (x1-x0) * f(x0,y0);
}



double backward_euler(double x0, double y0, double x1,
		      std::function<double (double, double)> f)
{
    // need two points: one is y1 from forward_euler
    double y1 = forward_euler(x0, y0, x1, f);
    //the other is y2, close to y1
    double y2 = y1 * 1.1;

    double constexpr tolerance{1e-9};

    auto const & g = [&x0,&y0,&x1,&f] (double y1) -> double
			 { return y1 - y0 - (x1-x0) * f(x1,y1); };

    double y = secant(y1, y2, tolerance, g);    

    return y;    
}

double g(double yinit)
{
    double xinit = 0.0;
    double xfin  = 2.0;

    double x0 = xinit;
    double y0 = yinit;

    double h = 1e-4;
    
    while (x0 < xfin) {
	double x1 = x0 + h;
	if (x1 > xfin) x1 = xfin;
	
	double y1 = backward_euler(x0,y0,x1,f);

	y0 = y1;
	x0 = x1;
    }

    return y0 - yinit;
}
    
int main()
{
    double yinit = secant(0.8, 1.2, 1e-9, g);

    std::cout << yinit << ' ' << g(yinit) << '\n';
}
    
