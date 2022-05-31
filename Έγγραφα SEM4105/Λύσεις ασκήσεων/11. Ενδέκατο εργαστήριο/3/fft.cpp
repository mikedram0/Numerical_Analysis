#include <cmath>
#include <cassert>
#include <complex>
#include <vector>
#include <fstream>
#include <cstddef>
#include <iterator>

#include <iostream>



double constexpr pi{3.14159265358979323846};

using complex = std::complex<double>;
using vec = std::vector<complex>;

constexpr bool isPowerOf2(std::size_t n)
{
    return n && !(n & (n - 1));
}



// with auxiliary vectors
void fft(vec const & f, vec & C)
{
    auto const N = f.size();
    assert(isPowerOf2(N));

    if (N == 1) {
	C[0] = f[0];

	return;
    }

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
    for (vec::size_type m{0}; m < N/2; ++m) {

	co[m] *= coef;
	
	C[m]     = 0.5*(ce[m] + co[m]);

	// C[m+N/2] = 0.5*(ce[m] - co[m]) = 0.5*(ce[m] + co[m]) - co[m];
	C[m+N/2] = C[m] - co[m];
	
	coef *= pol;
    }
    
}


// no auxiliary vectors
template<typename constiter, typename iter>
void fft(constiter fbeg, std::size_t N, iter cbeg, std::size_t inc = 1) 
{
    auto const n = N / inc;
    assert(isPowerOf2(n));

    if (n == 1) {
	*cbeg = *fbeg;
	return;
    }

    auto const Feven = fbeg;
    auto const Fodd = std::next(fbeg, inc);

    auto Ceven = cbeg; // C[m]
    auto Codd  = std::next(cbeg,n/2); // C[m+n/2]

    fft(Feven, N, Ceven, 2*inc);
    fft(Fodd,  N, Codd,  2*inc);

    auto const pol = std::polar(1.0, -2.0*pi/n);

    complex coef{1.0};
    for (std::size_t m{0}; m < n/2; ++m) {
	*Codd *= coef;

	*Ceven = 0.5*(*Ceven + *Codd);

	// C[m+N/2] = 0.5*(ce[m] - co[m]) = 0.5*(ce[m] + co[m]) - co[m];

	*Codd = *Ceven -  *Codd;

	++Ceven;
	++Codd;

	coef *= pol;
    }
}



int main()
{
    vec::size_type constexpr N{1024};

    vec f(N);    
    for (vec::size_type i{0} ; i < N; ++i) {
	f[i] = -0.5 + static_cast<double>(i)/N;
    }

    vec C(N);

//    fft(f,C);    
    fft(f.cbegin(),N,C.begin());
    
    std::ofstream out{"fourier.dat"};
    out.precision(15);
    for (vec::size_type i{0} ; i < N; ++i) {
	complex correct{0.0, 1.0/(2.0*i*pi)};
	if (i==0) correct = 0.0;
	out << C[i] << ' ' << correct << '\n';

    }
}
