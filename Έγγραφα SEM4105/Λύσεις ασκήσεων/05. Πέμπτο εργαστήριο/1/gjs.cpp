#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


using matrix = std::vector<double>;
using vec    = std::vector<double>;


bool check(vec const & y, double tol)
{
    return std::all_of(y.cbegin(), y.cend(),
		       [&tol] (auto const & z) {return std::abs(z) < tol;});
}


void gauss_jacobi(matrix const & a, vec const & b, vec & x)
{
    auto const n = b.size();
    
    for (auto & z : x) {
	z = 0.0;
    }
    
    while (true){
	vec y{b};
	
	for (vec::size_type i{0}; i < n; ++i) {
	    for (vec::size_type j{0}; j < n; ++j) {
		y[i] -= a[i+n*j] * x[j];
	    }
	    y[i] /= a[i+i*n];
	}
	
	for (vec::size_type i{0}; i < n; ++i) {
	    x[i] += y[i];
	}
	
	if (check(y,1e-7)) {
	    break;
	}
    }

}



void gauss_seidel(matrix const & a, vec const & b, vec & x)
{
    auto const n = b.size();
    
    for (auto & z : x) {
	z = 0.0;
    }

    while (true) {
	vec y{b};
	
	for (vec::size_type i{0}; i < n; ++i) {
	    for (vec::size_type j{0}; j < n; ++j) {
		y[i] -= a[i+n*j] * x[j];
	    }
	    y[i] /= a[i+i*n];
	    
	    x[i] += y[i];
	}
	
	if (check(y,1e-7)) {
	    break;
	}
    }
    
    
}


int main()
{
    std::size_t constexpr n{4};

    matrix a(n*n);
    
    a[0+n*0] = 12.1;
    a[0+n*1] = 3.9;
    a[0+n*2] = 0.3;
    a[0+n*3] = -4.1;

    a[1+n*0] = 4.3;
    a[1+n*1] = -11.3;
    a[1+n*2] = 0.8;
    a[1+n*3] = 1.5;

    a[2+n*0] = 1.0;
    a[2+n*1] = -2.8;
    a[2+n*2] = 14.3;
    a[2+n*3] = -8.1;
    
    a[3+n*0] = 2.4;
    a[3+n*1] = 6.1;
    a[3+n*2] = -1.1;
    a[3+n*3] = 12.5;


    vec b{1.2, 2.3, 3.4, 4.5};

    vec x(n);

    gauss_jacobi(a,b,x);
    
    std::cout << "solution with jacobi:\n";
    for (auto const & z : x) {
	std::cout << z << ' ';
    }
    std::cout << '\n';

    gauss_seidel(a,b,x);
    
    std::cout << "solution with seidel:\n";
    for (auto const & z : x) {
	std::cout << z << ' ';
    }
    std::cout << '\n';

}
