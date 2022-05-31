#include <cmath>
#include <iostream>
#include <complex>
#include <functional>


using complex = std::complex<double>;



double  g(double x, double y)
{
    return 4.0 * x*x - std::pow(y,3) + 28.0;
}


double h(double x, double y)
{
    return 3.0 * std::pow(x,3) + 4.0 * y*y - 145.0;
}


complex f(complex z)
{
    double x = std::real(z);
    double y = std::imag(z);

    double re = g(x,y);
    double im = h(x,y);

    return complex{re,im};
}


complex secant(complex x1, complex x2, double toler,
	      std::function<complex (complex)> func)
{
    complex f1{func(x1)};
    complex f2{func(x2)};
    complex x;
    
    do {
	x = (x2 * f1 - x1 * f2) / (f1 - f2);

	x1 = x2;
	f1 = f2;

	x2 = x;
	f2 = func(x);

    } while (std::abs(f2) > toler);

    return x;
}

int main()
{
    complex z1{0.0,0.0};
    complex z2{1.0,1.0};

    complex z = secant(z1, z2, 1e-8, f);

    std::cout << z.real() << ' ' << z.imag() <<'\n';
}
