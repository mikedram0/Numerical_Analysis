#include <cmath>
#include <iostream>
#include <complex>
#include <functional>



using complex = std::complex<double>;

// Find the roots of f(x)=0 with f(x)=sin(x) - x^2
// using the Muller's method.
// with complex arithmetic

complex f(complex x);

complex muller(complex x0, complex x1, complex x2,
	       std::function<complex (complex)> f)
{
    double constexpr epsilon{1e-6};
    double constexpr delta{1e-6};

    complex f0{f(x0)};
    complex f1{f(x1)};
    complex f2{f(x2)};

    complex x;
    
    while (true) {
	complex const w1{(f2-f1)/(x2-x1)};
	complex const w0{(f2-f0)/(x2-x0)};

	complex const a{(w1-w0)/(x1-x0)};

	complex const b{w0+a*(x2-x0)};
	complex const c{f2};

	complex const p{b*b-4.0*a*c};
	
	complex const d1{b + std::sqrt(p)};
	complex const d2{b - std::sqrt(p)};

	complex const d{std::abs(d1) > std::abs(d2) ? d1 : d2};

	x = x2 - 2.0*c / d;

	complex t{f(x)};
	    
	//        termination condition
	if ((std::abs(x - x2) < delta) || (std::abs(t) < epsilon)) {
	    break;
	}

	x0 = x1;
	x1 = x2;
	x2 = x;

	f0 = f1;
	f1 = f2;
	f2 = t;
    }

    return x;
}


int main()
{
    complex x0{0.1, 0.0};
    complex x1{0.2, 0.0};
    complex x2{0.5, 0.0};

    complex x = muller(x0,x1,x2, f);

    std::cout << u8"Η ρίζα είναι " << x << '\n';
    std::cout << u8"Η τιμή της συνάρτησης είναι " << f(x) <<'\n';
}


//     define f(x)
complex f(complex x)
{
    return std::sin(x) - x*x;
}
