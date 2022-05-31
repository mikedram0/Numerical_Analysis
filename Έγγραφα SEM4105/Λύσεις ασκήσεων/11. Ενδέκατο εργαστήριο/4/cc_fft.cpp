#include <cmath>
#include <cassert>
#include <complex>
#include <vector>
#include <fstream>
#include <cstddef>
#include <iterator>
#include <iostream>


using complex = std::complex<double>;
using vec = std::vector<complex>;




void fillv(vec & v) 
{
    auto n = v.size();

    for (vec::size_type k{0}; k < n/2; ++k) {
	v[k] = 2.0 / (1-4.0*k*k) -1.0 / (n*n-1.0 + n%2);
    }

    v[n/2] = (n-3.0)/( n - n%2 -1)
	- 1 + 1.0 / (n*n-1.0 + n%2) * ((2 - n%2)*n-1);
    
    for (vec::size_type k{1}; k <= (n-1)/2; ++k) {
	v[n-k] = v[k];
    }
}



double constexpr pi{3.14159265358979323846};

constexpr bool isPowerOf2(std::size_t n)
{
    return n && !(n & (n - 1));
}



// with auxiliary vectors
void fft(vec const & f, vec & C)
{
    auto const N = f.size();

    if (N == 1) {
	C[0] = f[0];

	return;
    }

    assert(isPowerOf2(N));

    vec fhalf(N/2), ce(N/2), co(N/2);

    for (vec::size_type i{0}; i < N/2; ++i) {
	fhalf[i] = f[2*i];
    }
    fft(fhalf, ce);

    for (vec::size_type i{0}; i < N/2; ++i) {
	fhalf[i] = f[2*i+1];
    }
    fft(fhalf, co);

    auto const pol = std::polar(1.0, -2.0*pi/N);
    complex coef{1.0};
    for (vec::size_type i{0}; i < N/2; ++i) {

	co[i] *= coef;
	
	C[i]     = 0.5*(ce[i] + co[i]);

	C[i+N/2] = C[i] - co[i];
	
	coef *= pol;
    }
    
}



// no auxiliary vectors
template<typename constiter, typename iter>
void fft(constiter fbeg, std::size_t N, iter cbeg,
	 std::size_t inc = 1) 
{
    auto const n = N / inc;

    if (n == 1) {
	*cbeg = *fbeg;
	return;
    }

    assert(isPowerOf2(n));

    auto const Feven = fbeg;
    auto const Fodd = std::next(fbeg, inc);

    auto Ceven = cbeg;
    auto Codd  = std::next(cbeg,n/2);

    fft(Feven, N, Ceven, 2*inc);
    fft(Fodd,  N, Codd,  2*inc);

    auto const pol = std::polar(1.0, -2.0*pi/n);

    complex coef{1.0};
    for (std::size_t i{0}; i < n/2; ++i) {
	*Codd *= coef;

	*Ceven = 0.5*(*Ceven + (*Codd));

	*Codd = *Ceven -  (*Codd);
	
	++Ceven;
	++Codd;

	coef *= pol;
    }
}



void coefs(vec & w)
{
    auto const n = w.size() - 1;

    vec f(n);

    fillv(f);

    
//    fft(f,w);
    fft(f.cbegin(), n, w.begin());

    w.back() = w.front();
}





/*
  \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x
*/
double f(double x)
{
    return 2.0/ (1.0+4.0*x*x);
}



double calc(std::size_t n)
{
    vec w(n+1);
    coefs(w);
    
    vec::value_type s{0.0};
    for (vec::size_type i{0}; i < w.size(); ++i) {
	double x = std::cos(i*pi/n);
	s += w[i] * f(x);
    }
    
    return s.real();
}
    

void check()
{
    vec::size_type constexpr n{1024};
    
    vec f(n);
    fillv(f);
    
    vec C(f.size());

//    fft(f,C);
    fft(f.cbegin(), n, C.begin());
    
    vec w(n+1);
    coefs(w);

    std::ofstream out{"wC.dat"};
    out.precision(15);
    for (std::size_t i=0; i < w.size(); ++i) {
	out << w[i] << ' ' << C[i] << ' ' << w[i]-C[i] << '\n';
    }
}

int main()
{
    check();
    
    double const correct = 2.0*std::atan(2.0);
    std::cout.precision(15);
    std::size_t n = 2;
    while (true) {
	double s = calc(n);
	double diff = s-correct;
	std::cout << s << ' ' << correct << ' ' << diff << '\n';
	if (std::abs(diff) < 1e-12) {
	    break;	    
	}
	n*=2;
    }
    std::cout << n << '\n';
}

