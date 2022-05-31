#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>


using vec=std::vector<double>;


void triang(vec & a, vec & b)
{
    auto const n = b.size();
    
    for (std::size_t k{0}; k < n-1; ++k) {
	for (std::size_t i{k+1}; i < n; ++i) {
	    double ell{-a[i+n*k]/a[k+n*k]};
	    for (std::size_t j{k}; j < n; ++j) {
		a[i+n*j] += ell * a[k+n*j];
	    }
	    b[i] += ell * b[k];
	}
    }
}


void backsub(vec const & a, vec const & b, vec & x)
{
    auto const n = b.size();

    for (std::size_t i{n}; i-- > 0; ) {
	x[i] = b[i];
	for (std::size_t j{i+1}; j < n; ++j) {
	    x[i] -= a[i+n*j] * x[j];
	}
	x[i] /= a[i+n*i];
    }
}


int  main()
{
    std::ifstream in{"gauss1.txt"};
    std::size_t n;
    in >> n;

    vec a(n*n);
    vec b(n);
    vec x(n);

    for (std::size_t i{0}; i<n; ++i) {
	for (std::size_t j{0}; j<n; ++j) {
	    in >> a[i+j*n];
	}
    }

    for (auto & z : b) {
	in >> z;
    }


    triang(a,b);
    backsub(a,b,x);

    for (auto const & z : x) {
	std::cout << z << '\n';
    }
    
}
