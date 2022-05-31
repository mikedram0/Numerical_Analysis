#include <vector>
#include <iostream>
#include <cmath>


double constexpr pi{3.14159265358979323846};


constexpr double c(int i, int n)
{
    return (i%n == 0) ? 1.0 : 2.0;
}

constexpr double b(int j, int n)
{
    return (j%(n/2) == 0) ? 1.0 : 2.0;
}




void coefs(std::vector<double> & w)
{
    auto const n = w.size() - 1;
    
    for (int i=0; i<=n; ++i) {
	w[i] = 0.0;
	for (int j=0; j <= n/2; ++j) {
	    w[i] += b(j,n) / (1-4.0*j*j) * std::cos(2.0*i*j*pi/n);
	}
	
	w[i] *= c(i,n)/n;
    }
}


/*
  \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x
 */
double f(double x)
{
    return 2.0/ (1.0+4.0*x*x);
}



double calc(int n)
{
    std::vector<double> w(n+1);
    coefs(w);
    
    double s{0.0};
    for (int i=0; i < w.size(); ++i) {
	double x = std::cos(i*pi/n);
	s += w[i] * f(x);
    }
    
    return s;
}
    
int main()
{
    double const correct = 2.0*std::atan(2.0);
    std::cout.precision(15);
    int n = 3;
    while (true) {
	double s = calc(n);
	double diff = s-correct;
	std::cout << s << ' ' << correct << ' ' << diff << '\n';
	if (std::abs(diff) < 1e-12) {
	    break;	    
	}
	++n;	
    }
    std::cout << n << '\n';
}


