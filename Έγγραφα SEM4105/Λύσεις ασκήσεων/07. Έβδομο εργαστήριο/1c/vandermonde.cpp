#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>


using vec=std::vector<double>;


void normalize(vec & a, vec & b, vec::size_type k)
{
    auto const n = b.size();

    for (auto i = k; i < n; ++i) {
	auto maxv = std::abs(b[i]);
	for (auto j = k; j < n; ++j) {
	    auto const v = std::abs(a[i+n*j]);
	    if (v > maxv) {
		maxv = v;
	    }
	}
	
	for (auto j = k; j < n; ++j) {
	    a[i+n*j] /= maxv;
	}
	b[i] /= maxv;
    }
}


void swap(vec & a, vec & b, vec::size_type k, vec::size_type p)
{
    auto const n = b.size();
    
    for (auto j = k; j < n; ++j) {
	std::swap(a[p + n*j], a[k + n *j]);
    }
    std::swap(b[p], b[k]);
}


void pivot(vec & a, vec & b, vec::size_type k)
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


void backsub(vec const & a, vec const & b, vec & x)
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



double pol(vec const & a, double xbar)
{
    double r{0.0};
    for (vec::size_type i{0}; i < a.size(); ++i) {
	r += a[i] * std::pow(xbar, i);
    }
    return r;
}

int main()
{
    std::ifstream in{"points.dat"};
    vec::size_type n;
    in >> n;

    vec x(n);
    vec y(n);
    for (vec::size_type i{0}; i<n; ++i) {
	in >> x[i]  >> y[i];
    }
    
    vec a(n*n);
    vec coef(n);

    for (vec::size_type i{0}; i<n; ++i) {
	for (vec::size_type j{0}; j<n; ++j) {
	    a[i+j*n] = std::pow(x[i], j);
	}
    }

    triang(a,y);
    backsub(a,y,coef);

    int constexpr m{100};

    auto const p = std::minmax_element(x.cbegin(), x.cend());

    double const & min = *p.first;
    double const & max = *p.second;
    double const step = (max-min)/(m-1);

    
    for (int i{0}; i < m; ++i) {
	double xout = min + step * i;
	double yout = pol(coef,xout);

	std::cout << xout << ' ' << yout << '\n';
    }
    
}
