#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>

double constexpr pi{3.14159265358979323846};


using matrix=std::vector<double>;
using vec=std::vector<double>;


void normalize(matrix & a, vec & b, vec::size_type k)
{
    auto const n = b.size();

    for (auto j = k; j < n; ++j) {
	auto maxv = std::abs(b[j]);
	for (auto p = k; p < n; ++p) {
	    auto const v = std::abs(a[j+n*p]);
	    if (v > maxv) {
		maxv = v;
	    }
	}
	
	for (auto p = k; p < n; ++p) {
	    a[j+n*p] /= maxv;
	}
	b[j] /= maxv;
    }
}


void swap(matrix & a, vec & b, vec::size_type k, vec::size_type p)
{
    auto const n = b.size();
    
    for (auto j = k; j < n; ++j) {
	std::swap(a[p + n*j], a[k + n *j]);
    }
    std::swap(b[p], b[k]);
}


void pivot(matrix & a, vec & b, vec::size_type k)
{
    auto const n = b.size();

    normalize(a,b,k);

    // 2
    auto maxp = k;
    for (auto p = k+1; p < n; ++p) {
	if (std::abs(a[p+n*k]) > std::abs(a[maxp+n*k])) {
	    maxp = p;
	}
    }

    // 3
    if (maxp != k) {
	swap(a,b,k,maxp);
    }
}

void triang(matrix & a, vec & b)
{
    auto const n = b.size();
    
    for (vec::size_type k{0}; k < n-1; ++k) {
	pivot(a,b,k);
	
	for (vec::size_type i{k+1}; i < n; ++i) {
	    auto  ell = -a[i+n*k]/a[k+n*k];
	    for (vec::size_type j{k}; j < n; ++j) {
		a[i+n*j] += ell * a[k+n*j];
	    }
	    b[i] += ell * b[k];
	}
    }
}


void backsub(matrix const & a, vec const & b, vec & x)
{
    auto const n = b.size();

    for (auto i = n; i-- > 0; ) {
	x[i] = b[i];
	for (auto j = i+1; j < n; ++j) {
	    x[i] -= a[i+n*j] * x[j];
	}
	x[i] /= a[i+n*i];
    }
}


void gauss(matrix & a, vec & b, vec &x)
{
    triang(a, b);
    backsub(a,b,x);
}


/* Gauss elimination with pivoting: end */



void sys(vec const & x, double a, double b,  matrix & A, vec & B)
{
    auto n = x.size();
    
    for (vec::size_type i{0}; i < n; ++i) {
	for (vec::size_type j{0}; j < n; ++j) {
	    A[i+j*n] = std::pow(x[j],i); 
	}
	B[i] = (std::pow(b,i+1)-std::pow(a,i+1))/(i+1);
    }

}

double f(double x)
{
    return std::pow(x,3) * std::sin(pi * x);
}

int main()
{
    vec const x{-0.9,-0.7,-0.4,0.1,0.4,0.8,0.9};
    
    auto n = x.size();
    
    matrix a(n*n);
    vec b(n);
    vec w(n);

    double alimit = -1.0;
    double blimit = 1.0;
    
    sys(x, alimit, blimit, a, b);

    gauss(a,b,w);

    double s{0.0};
    for (vec::size_type i{0}; i < n; ++i) {
	s += w[i] * f(x[i]);
    }

    double constexpr correct = (2.0 -12.0 /(pi*pi)) / pi;

    std::cout << s << ' ' << correct <<'\n';
}
