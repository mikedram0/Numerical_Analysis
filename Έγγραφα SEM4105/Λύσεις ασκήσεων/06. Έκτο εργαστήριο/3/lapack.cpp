#include <complex>
#include <vector>
#include <iostream>


using complex = std::complex<double>;

extern "C"
void 
dgeev_(const char & jobvl, const char & jobvr, const int & n, 
       double a[], const int & lda, double wr[], double wi[], double vl[], 
       const int & ldvl, double vr[], const int & ldvr, double work[], 
       const int & lwork, int & info);


    
void eival(std::vector<double> & a, int n, std::vector<complex> & w)
{
    std::vector<double> wr(n);
    std::vector<double> wi(n);
    int constexpr ldvl{1};
    std::vector<double>  vl(ldvl*n);
    int constexpr ldvr{1};
    std::vector<double>  vr(ldvr*n);

    double temp;
    int lwork = -1;
    int info = 0;
    
    dgeev_('n', 'n', n, a.data(), n, wr.data(), wi.data(),
	   vl.data(), ldvl, vr.data(), ldvr, &temp, lwork, info);
    
    lwork = static_cast<int>(temp);
    std::vector<double> work(lwork);

    info = 0;
    
    dgeev_('n', 'n', n, a.data(), n, wr.data(), wi.data(),
	   vl.data(), ldvl, vr.data(), ldvr, work.data(), lwork, info);

    if (info == 0) {
	for(int i=0; i < n; ++i) {
	    w[i] = complex{wr[i], wi[i]};
	}
    }
}


int main()
{
    int constexpr n{3};
    std::vector<double> a(n*n);

    a[0+3*0] = 6.3;
    a[0+3*1] = 2.1;
    a[0+3*2] = 4.15;

    a[1+3*0] = 3.1;
    a[1+3*1] = 5.14;
    a[1+3*2] = 1.03;

    a[2+3*0] = -11.0;
    a[2+3*1] = 12.3;
    a[2+3*2] = -8.8;

    std::vector<complex> w(n);
    
    eival(a,n,w);

    for (auto const & z : w) {
	std::cout << z.real() << ' ' << z.imag() <<'\n';
    }
}
