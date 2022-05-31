#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>



using matrix=std::vector<double>;
using vec=std::vector<double>;

/* Gauss elimination with pivoting: start */

void normalize(matrix & a, vec & b, vec::size_type k)
{
    auto const n = b.size();

    for (auto i = k; i < n; ++i) {
	auto maxv = std::abs(b[i]);
	for (auto p = k; p < n; ++p) {
	    auto const v = std::abs(a[i+n*p]);
	    if (v > maxv) {
		maxv = v;
	    }
	}
	
	for (auto p = k; p < n; ++p) {
	    a[i+n*p] /= maxv;
	}
	b[i] /= maxv;
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

void triang(matrix & a, vec & b)
{
    auto const n = b.size();
    
    for (vec::size_type k{0}; k < n-1; ++k) {
	pivot(a,b,k);
	
	for (vec::size_type i{k+1}; i < n; ++i) {
	    auto ell = -a[i+n*k]/a[k+n*k];
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




void read(std::string const fname, vec & x, vec & y)
{
    std::ifstream in{fname};

    int n;
    in >> n;

    x.reserve(n);
    y.reserve(n);
    
    for (int i{0}; i < n; ++i) {
	double a,b;
	in >> a >> b;

	x.push_back(a);
	y.push_back(b);
    }
    
}


// f^{(order)}_i (x)
double f(int i, int order, double x)
{
    double v;
    
    switch (order) {
    case 0:    // => function value
	v = std::pow(x,i);
	break;
	
    case 1:    // => first derivative
	if (i == 0) {
	    v = 0.0;
	} else {
	    v= i*std::pow(x,i-1);
	}
	break;

    case 2: // => second derivative
	if (i == 0 || i == 1) {
	    v = 0.0;
	} else {
	    v= i*(i-1)*std::pow(x,i-2);
	}
	break;
    default:
	std::cerr << "not implemented\n";
	std::exit(-1);

    }

    return v;
}


void sys(vec const & x, double xbar, int order, matrix & a, vec & b)
{
    auto n = x.size();
    
    for (int i{0}; i < n; ++i) {
	for (int j{0}; j < n; ++j) {
	    a[i+j*n] = f(i, 0, x[j]);
	}
	b[i] = f(i,order,xbar);
    }

}

double correct(int order, double x)
{
    switch (order) {
    case 0:
	return std::pow(std::sin(x),2);
    case 1:
	return 2*std::sin(x)*std::cos(x);
    case 2:
	return 2*std::cos(2*x);
    default:
	std::cerr << "not implemented\n";
	std::exit(-1);
    }
}

int main()
{
    vec x,y;
    read("points2.dat", x, y);

    auto n = x.size();
    
    matrix a(n*n);
    vec b(n);

    vec sol(n);
    double constexpr xbar{2.5};
    
    for (int order{1}; order <= 2; ++order) {
	sys(x, xbar, order, a, b);
	
	gauss(a,b,sol);
	
	double s{0.0};
	for (int i{0}; i < n; ++i) {
	    s+=sol[i] * y[i];
	}
	std::cout.precision(15);
	
	std::cout << order << " derivative at " << xbar
		  << " is " << s << '\n';
	std::cout << "difference from correct value "
		  << s - correct(order,xbar) << '\n';
    }
    
}
