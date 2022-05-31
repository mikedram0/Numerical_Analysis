#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <cstddef>


using matrix=std::vector<double>;
using vec=std::vector<double>;

bool iszero(matrix::value_type x)
{
    return std::abs(x) < 1e-10;
}


// maxval in row i, from k column to end
matrix::value_type rowmaxval(std::size_t n, matrix const & a, 
			     vec::size_type k, vec::size_type i)
{
    auto maxv = std::abs(a[i+n*k]);
			 
    for (auto j = k+1; j < n; ++j) {
	auto const v = std::abs(a[i+n*j]);
	if (v > maxv) {
	    maxv = v;
	}
    }
    
    return maxv;
}

void normalize(std::size_t n, matrix & a, vec::size_type k,
	       matrix::value_type & p)
{
    for (auto i = k; i < n; ++i) {
	auto const maxv = rowmaxval(n, a, k, i);
	
	if (iszero(maxv)) {
	    continue;
	}
	
	for (auto j = k; j < n; ++j) {
	    a[i+n*j] /= maxv;
	}
	p *= maxv;
    }
}


void swap(std::size_t n, matrix & a, vec::size_type k, vec::size_type i)
{
    for (auto j = k; j < n; ++j) {
	std::swap(a[i + n*j], a[k + n *j]);
    }
}


void pivot(std::size_t n, matrix & a, vec::size_type k, std::size_t & change,
	   matrix::value_type & p)
{
    normalize(n,a,k,p);

    auto maxp = k;
    for (auto i = k+1; i < n; ++i) {
	if (std::abs(a[i+n*k]) > std::abs(a[maxp+n*k])) {
	    maxp = i;
	}
    }

    if (maxp != k) {
	++change;
	swap(n,a,k,maxp);
    }
}

void triang(std::size_t n, matrix & a, std::size_t & change,
	    	       matrix::value_type & p)
{
    change = 0;
    p = 1.0;
    
    for (vec::size_type k{0}; k < n-1; ++k) {
	pivot(n,a,k, change, p);
	
	for (auto i = k+1; i < n; ++i) {
	    auto const ell = -a[i+n*k]/a[k+n*k];
	    a[i+n*k]=0;
	    for (auto j = k+1; j < n; ++j) {
		a[i+n*j] += ell * a[k+n*j];
	    }
	}
    }
}



double det(std::size_t n, matrix & a)
{
    std::size_t change;
    matrix::value_type p;
    triang(n,a,change, p);

    double d{change %2 == 0 ? p : -p};
    for (vec::size_type i{0}; i < n; ++i) {
	d *= a[i+i*n];
    }
    
    return d;
}


int
main()
{
    std::size_t constexpr n{4};
    vec a(n*n);

    a[0+n*0] = 2.1;
    a[0+n*1] = 3.9;
    a[0+n*2] = 0.3;
    a[0+n*3] = -4.1;
    a[1+n*0] = 4.3;
    a[1+n*1] = -1.3;
    a[1+n*2] = 0.8;
    a[1+n*3] = 1.5;
    a[2+n*0] = 1.0;
    a[2+n*1] = -2.8;
    a[2+n*2] = 4.3;
    a[2+n*3] = -8.1;
    a[3+n*0] = 2.4;
    a[3+n*1] = 6.1;
    a[3+n*2] = -1.1;
    a[3+n*3] = 12.5;

    std::cout << det(n,a) << '\n';
}

