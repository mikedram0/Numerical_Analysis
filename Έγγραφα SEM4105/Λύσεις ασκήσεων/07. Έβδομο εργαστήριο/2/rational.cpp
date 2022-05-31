#include <vector>
#include <cmath>
#include <iostream>

using matrix=std::vector<double>;
using vec=std::vector<double>;

/*
  R(x_i) = y_i => P(x_i) = y_i Q(x_i) =>
  a_0 + a_1 x_i + a_2 x_i^2 = y_i (1 + b_1 x_i + b_2 x_i^2 + b_3 x_i^3 =>
  1 a_0 + x_i a_1 + x_i^2 a_2 - y_i x_i b_1 - y_i x_i^2 b_2 - y_i x_i^3 b_3 = y_i

  | 1  x_1   x_1^2   -y_1*x_1   -y_1*x_1^2  -y_1*x_1^3 |  | a_0 |    | y_1 |
  | 1  x_2   x_2^2   -y_2*x_2   -y_2*x_2^2  -y_2*x_2^3 |  | a_1 |    | y_2 |
  | 1  x_3   x_3^2   -y_3*x_3   -y_3*x_3^2  -y_3*x_3^3 |  | a_2 |    | y_3 |
  | 1  x_4   x_4^2   -y_4*x_4   -y_4*x_4^2  -y_4*x_4^3 |* | b_1 | =  | y_4 |
  | 1  x_5   x_5^2   -y_5*x_5   -y_5*x_5^2  -y_5*x_5^3 |  | b_2 |    | y_5 |
  | 1  x_6   x_6^2   -y_6*x_6   -y_6*x_6^2  -y_6*x_6^3 |  | b_3 |    | y_6 |

*/

double R(const vec & c, double x)
{
    double p = c[0] + c[1] * x + c[2] *x*x;
    double q = 1.0 + c[3] * x + c[4] * x*x + c[5] * x*x*x;

    return p/q;
}


bool iszero(matrix::value_type x)
{
    return std::abs(x) < 1e-10;
}


// maxval in row i, from k column to end
matrix::value_type rowmaxval(matrix const & a, vec const & b, 
			     vec::size_type k, vec::size_type i)
{
    auto const n = b.size();

    auto maxv = std::abs(b[i]);

    for (auto j = k; j < n; ++j) {
	auto const v = std::abs(a[i+n*j]);
	if (v > maxv) {
	    maxv = v;
	}
    }

    return maxv;
}

void normalize(matrix & a, vec & b, vec::size_type k)
{
    auto const n = b.size();

    for (auto i = k; i < n; ++i) {
	auto const maxv = rowmaxval(a,b, k, i);

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
    auto const n = b.size();
    
    for (auto j = k; j < n; ++j) {
	std::swap(a[i + n*j], a[k + n *j]);
    }
    std::swap(b[i], b[k]);
}


void pivot(matrix & a, vec & b, vec::size_type k)
{
    auto const n = b.size();

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
    auto const n = b.size();
    
    for (vec::size_type k{0}; k < n-1; ++k) {
	pivot(a,b,k);
	
	for (auto i = k+1; i < n; ++i) {
	    auto const ell = -a[i+n*k]/a[k+n*k];
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

void backsub(matrix const & a, vec const & b, vec & x)
{
    auto const n = b.size();

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
    constexpr std::size_t n{6};

    matrix a(n*n);
    vec x{0.9, 1.1, 1.5, 2.0, 2.9, 3.5};
    vec y{5.607, 4.576, 3.726, 3.354, 3.14, 3.087};
    vec b{y};
    vec c(n);
	
    for (std::size_t i{0}; i < n; ++i) {
	a[i+0*n] = 1.0;
	a[i+1*n] = x[i];
	a[i+2*n] = x[i]*x[i];
	a[i+3*n] = -y[i]*x[i];
	a[i+4*n] = -y[i]*x[i]*x[i];
	a[i+5*n] = -y[i]*x[i]*x[i]*x[i];
    }

	
    // solve A*C=B to get C, the coefficients a,b
    gauss(a,b,c);

    std::cout << "a_0  a_1 a_2 = "
	      << c[0] << ' ' << c[1] << ' ' << c[2] <<'\n';
    std::cout << "b_1  b_2 b_3 = "
	      << c[3] << ' ' << c[4] << ' ' << c[5] <<'\n';

    for (vec::size_type i{0}; i < n; ++i) {
	std::cout << R(c,x[i]) << ' ' <<  y[i] <<'\n';
    }
}


