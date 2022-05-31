#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>

using vec = std::vector<double>;


void normalize(vec & a, vec & b, std::size_t k)
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
	
	for (auto  p = k; p < n; ++p) {
	    a[j+n*p] /= maxv;
	}
	b[j] /= maxv;
    }
}


void swap(vec & a, vec & b, std::size_t k, std::size_t p)
{
    auto const n = b.size();
    
    for (auto j = k; j < n; ++j) {
	std::swap(a[p + n*j], a[k + n *j]);
    }
    std::swap(b[p], b[k]);
}


void pivot(vec & a, vec & b, std::size_t k)
{
    auto const n = b.size();

    normalize(a,b,k);

    auto maxp = k;
    for (auto p = k+1; p < n; ++p) {
	if (std::abs(a[p+n*k]) > std::abs(a[maxp+n*k])) {
	    maxp = p;
	}
    }

    if (maxp != k) {
	swap(a,b,k,maxp);
    }
}

void triang(vec & a, vec & b)
{
    auto const n = b.size();
    
    for (std::size_t k{0}; k < n-1; ++k) {
	pivot(a,b,k);
	
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


void gauss(vec & a, vec & b, vec &x)
{
    triang(a, b);
    backsub(a,b,x);
}



void nonlinear(vec const & x, vec & f)
{
    f[0] = x[0]             + x[1]             + x[2] - 3.0;
    f[1] = x[0]*x[0] * x[1] + x[1]*x[1] * x[2] + x[2]*x[2] * x[0] - 4.0;
    f[2] = x[0]*x[0]        + x[1]*x[1]        + x[2]*x[2]  - 5.0;
}

void linear(vec const & x, vec & a)
{
    auto const n = x.size();
    
    a[0+n*0] = 1.0;
    a[0+n*1] = 1.0;
    a[0+n*2] = 1.0;

    a[1+n*0] = 2.0 * x[0] * x[1] + x[2]*x[2];
    a[1+n*1] = 2.0 * x[1] * x[2] + x[0]*x[0];
    a[1+n*2] = 2.0 * x[2] * x[0] + x[1]*x[1];

    a[2+n*0] = 2.0 * x[0];
    a[2+n*1] = 2.0 * x[1];
    a[2+n*2] = 2.0 * x[2];
}


bool check(vec const & x, double tol)
{
    for (auto const & z : x) {
	if (std::abs(z) > tol) {
	    return false;
	}
    }
    return true;
}


int main()
{
    int constexpr n{3};
    double constexpr tol{1e-7};
    
    vec x{1.0,1.5, 2.1};
    
    while (true) {
	vec f(n);
	nonlinear(x, f);
	
	if (check(f, tol)) {
	    break;
	}

	vec a(n*n);
	linear(x,a);

	vec y(n);
	gauss(a,f,y);

	for (int i{0}; i < n; ++i) {
	    x[i] -= y[i];
	}
    }

    for (auto const & z : x) {
	std::cout << z << '\n';
    }
}
