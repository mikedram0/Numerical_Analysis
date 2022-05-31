#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <functional>

using matrix=std::vector<double>;



void normalize(matrix & a, matrix::size_type n, matrix::size_type k)
{
    for (auto i = k; i < n; ++i) {
	auto maxv = std::abs(a[i+n*k]);
	for (auto j = k+1; j < n; ++j) {
	    auto const v = std::abs(a[i+n*j]);
	    if (v > maxv) {
		maxv = v;
	    }
	}
	
	for (auto j = k; j < n; ++j) {
	    a[i+n*j] /= maxv;
	}
    }
}


void swap(matrix & a, matrix::size_type n, matrix::size_type k, matrix::size_type p)
{
    for (auto j = k; j < n; ++j) {
	std::swap(a[p + n*j], a[k + n *j]);
    }
}


void pivot(matrix & a, matrix::size_type n, matrix::size_type k,
	   matrix::size_type & change)
{
    normalize(a,n,k);

    auto maxp = k;
    for (auto p = k+1; p < n; ++p) {
	if (std::abs(a[p+n*k]) > std::abs(a[maxp+n*k])) {
	    maxp = p;
	}
    }

    if (maxp != k) {
	++change;
	swap(a,n,k,maxp);
    }
}

void triang(matrix & a, matrix::size_type n, matrix::size_type & change)
{
    for (matrix::size_type k{0}; k < n-1; ++k) {
	pivot(a,n,k, change);
	
	for (matrix::size_type i{k+1}; i < n; ++i) {
	    auto ell = -a[i+n*k]/a[k+n*k];
	    for (matrix::size_type j{k}; j < n; ++j) {
		a[i+n*j] += ell * a[k+n*j];
	    }
	}
    }
}




double det(matrix & a, matrix::size_type n)
{
    matrix::size_type change{0};
    triang(a,n,change);

    double d{1.0};
    for (int i{0}; i < n; ++i) {
	d *= a[i+i*n];
    }
    
    return (change%2 == 0? d : -d);
}



double f(matrix::size_type n, double x, matrix a)
{
    for(int i{0}; i < n; ++i) {
	a[i+n*i] -= x;
    }

    return det(a,n);
}



double secant(double x1, double x2, double toler,
	      std::function<double (matrix::size_type, double, matrix)> func,
	      matrix::size_type n, matrix const & a)
{
    double f1{func(n,x1,a)};
    double f2{func(n,x2,a)};
    double x;
    
    do {
	x = (x2 * f1 - x1 * f2) / (f1 - f2);

	x1 = x2;
	f1 = f2;

	x2 = x;
	f2 = func(n,x,a);

    } while (std::abs(f2) > toler);

    return x;
}


int main()
{
    int constexpr n{4};
    matrix a(n*n);

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

    double x1{-1.0};
    double x2{1.0};

    double x = secant(x1, x2, 1e-8, f, n, a);

    std::cout << u8"Μια ιδιοτιμή είναι η "<< x
	      << u8" με τιμή ορίζουσας " << f(n,x,a) <<'\n';
}

