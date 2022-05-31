#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

using matrix=std::vector<double>;
using vec=std::vector<double>;

bool iszero(matrix::value_type x)
{
    return std::abs(x) < 1e-10;
}


// maxval in row i, from k column to end
matrix::value_type rowmaxval(const matrix & a, const vec & b, 
			     vec::size_type k, vec::size_type i)
{
    const auto n = b.size();

    auto maxv = std::abs(b[i]);

    for (auto j = k; j < n; ++j) {
	const auto v = std::abs(a[i+n*j]);
	if (v > maxv) {
	    maxv = v;
	}
    }

    return maxv;
}

void normalize(matrix & a, vec & b, vec::size_type k)
{
    const auto n = b.size();

    for (auto i = k; i < n; ++i) {
	const auto maxv = rowmaxval(a,b, k, i);

	if (iszero(maxv)) {
	    continue;
	}
	
	for (auto j = k; j < n; ++j) {
	    a[i+n*j] /= maxv;
	}
	b[i] /= maxv;
    }
}


void swap(matrix & a, vec & b, vec::size_type k, vec::size_type i)
{
    const auto n = b.size();
    
    for (auto j = k; j < n; ++j) {
	std::swap(a[i + n*j], a[k + n *j]);
    }
    std::swap(b[i], b[k]);
}


void pivot(matrix & a, vec & b, vec::size_type k)
{
    const auto n = b.size();

    normalize(a,b,k);

    auto maxp = k;
    for (auto i = k+1; i < n; ++i) {
	if (std::abs(a[i+n*k]) > std::abs(a[maxp+n*k])) {
	    maxp = i;
	}
    }

    if (maxp != k) {
	swap(a,b,k,maxp);
    }
}

void triang(matrix & a, vec & b)
{
    const auto n = b.size();
    
    for (vec::size_type k{0}; k < n-1; ++k) {
	pivot(a,b,k);
	
	for (auto i = k+1; i < n; ++i) {
	    const auto ell = -a[i+n*k]/a[k+n*k];
	    a[i+n*k]=0;
	    for (auto j = k+1; j < n; ++j) {
		a[i+n*j] += ell * a[k+n*j];
	    }
	    b[i] += ell * b[k];
	}
    }
}

// ax=b
matrix::value_type solve_linear(matrix::value_type a,
				matrix::value_type b)
{
    if (iszero(a)) {
	if (iszero(b)) {
	    std::cerr << u8"Το σύστημα έχει άπειρες λύσεις\n";
	} else {
	    std::cerr << u8"Το σύστημα δεν έχει λύση\n";
	}
    }

    return b/a;
}

void backsub(const matrix & a, const vec & b, vec & x)
{
    const auto n = b.size();

    for (auto i = n; i-- > 0; ) {
	x[i] = b[i];
	for (auto j = i+1; j < n; ++j) {
	    x[i] -= a[i+n*j] * x[j];
	}
	x[i] = solve_linear(a[i+n*i], x[i]);
    }
}


void gauss(matrix & a, vec & b, vec & x)
{
    triang(a, b);
    backsub(a,b,x);
}


int main()
{
    std::ifstream in{"gauss2.txt"};
    vec::size_type n;
    in >> n;

    matrix a(n*n);
    vec b(n);
    vec x(n);

    for (vec::size_type i{0}; i<n; ++i) {
	for (vec::size_type j{0}; j<n; ++j) {
	    in >> a[i+j*n];
	}
    }

    for (auto & z : b) {
	in >> z;
    }

    gauss(a,b,x);

    for (const auto & z : x) {
	std::cout << z << '\n';
    }
    
}
