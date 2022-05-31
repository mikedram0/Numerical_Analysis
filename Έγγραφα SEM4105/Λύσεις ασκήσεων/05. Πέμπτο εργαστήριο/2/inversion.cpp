#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstddef>


using matrix=std::vector<double>;


void print(matrix const & a, matrix::size_type n)
{
    for (matrix::size_type i{0} ; i < n; ++i) {
	for (matrix::size_type j{0}; j < n; ++j) {
	    std::cout << a[i+n*j] << ' ';
	}
	std::cout << '\n';	
    }
    std::cout << '\n';
}

bool iszero(matrix::value_type x)
{
    return std::abs(x) < 1e-10;
}


void normalize(matrix & a, matrix & b,
	       matrix::size_type n, matrix::size_type k)
{
    for (matrix::size_type i{k}; i < n; ++i) {
	matrix::value_type maxv = -1.0;

	for (matrix::size_type j{0}; j < n; ++j) {
	    auto const v = std::abs(b[i+n*j]);
	    if (v > maxv) {
		maxv = v;
	    }
	}

	for (matrix::size_type j{k}; j < n; ++j) {
	    auto const v = std::abs(a[i+n*j]);
	    if (v > maxv) {
		maxv = v;
	    }
	}
	
	for (matrix::size_type j{0}; j < n; ++j) {
	    b[i+n*j] /= maxv;
	}

	for (matrix::size_type j{k}; j < n; ++j) {
	    a[i+n*j] /= maxv;
	}
    }
}


void swap(matrix & a, matrix & b, matrix::size_type n, 
	  matrix::size_type k, matrix::size_type p)
{
    for (matrix::size_type j{k}; j < n; ++j) {
	std::swap(a[p + n*j], a[k + n *j]);
    }

    for (matrix::size_type j{0}; j < n; ++j) {
	std::swap(b[p + n*j], b[k + n *j]);
    }
}


void pivot(matrix & a, matrix & b,
	   matrix::size_type n, matrix::size_type k)
{
    normalize(a,b,n,k);

    auto maxp = k;
    for (auto p = k+1; p < n; ++p) {
	if (std::abs(a[p+n*k]) > std::abs(a[maxp+n*k])) {
	    maxp = p;
	}
    }

    if (maxp != k) {
	swap(a,b,n,k,maxp);
    }
}



void triang(matrix & a, matrix & b, matrix::size_type n)
{    
    for (matrix::size_type k{0}; k < n-1; ++k) {
	pivot(a,b,n,k);
	
	for (matrix::size_type i{k+1}; i < n; ++i) {
	    auto const ell = -a[i+n*k]/a[k+n*k];
	    for (matrix::size_type j{k}; j < n; ++j) {
		a[i+n*j] += ell * a[k+n*j];
	    }
	    for (matrix::size_type j{0}; j < n; ++j) {
		b[i+n*j] += ell * b[k+n*j];
	    }
	}
    }
}


void backsub(matrix const & a, matrix const & b, matrix::size_type n,
	     matrix & x)
{
    for (matrix::size_type k{0}; k < n; ++k) {
	for (auto i = n; i-- > 0; ) {
	    x[i+n*k] = b[i+n*k];
	    for (matrix::size_type j{i+1}; j < n; ++j) {
		x[i+n*k] -= a[i+n*j] * x[j+n*k];
	    }

	    if (iszero(a[i+n*i])) {
		std::cerr << u8"Ο πίνακας δεν έχει αντίστροφο\n";
		std::exit(-1);
	    }
	    x[i+n*k] /= a[i+n*i];	
	}
    }
}




void gauss(matrix & a, matrix & b, matrix::size_type n, matrix & x)
{
    triang(a, b, n);
    backsub(a, b, n, x);
}




void identity(matrix & b)
{    for (auto & z : b) {
	z = 0.0;	
    }
    
    for (matrix::size_type i{0}; i < n; ++i) {
	b[i+n*i] = 1.0;
    }
}    


/*
  a : Ο πίνακας (δεδομένο)
  b : Ο αντίστροφός του (το "αποτέλεσμα" της υπορουτίνας).
*/
void invert(matrix & a, matrix & b, matrix::size_type n)
{
    identity(b);
    matrix x(n*n);
    gauss(a,b, n, x);
    b=x;
}


int main()
{
    matrix::size_type constexpr n{4};
    matrix a(n*n);
    matrix ainv(n*n);

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

    invert(a, ainv, n);
    print(ainv,n);
}
