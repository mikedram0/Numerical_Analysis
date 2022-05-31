#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

#include <vector>
#include <cassert>
#include <algorithm>



using matrix = std::vector<double>;
using vec=std::vector<double>;


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
	    auto const ell = -a[i+n*k]/a[k+n*k];
	    a[i+n*k]=0;
	    for (vec::size_type j{k+1}; j < n; ++j) {
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



 /*

Για i=0,...,n-1 :   
    dx_i = x_{i+1}-x_i

Για i=0,...,n-1 :

    a_i dx_i^2 + b_i dx_i + c_i = (f_{i+1} - f_i) / dx_i


Για i=0,...,n-2 :

    3 a_i dx_i^2 + 2b_i dx_i + c_i - c_{i+1} = 0


Για i=0,...,n-2 :

    3a_i dx_i + b_i - b_{i+1} = 0



b_0 = 0

3 a_{n-1} dx_{n-1} + b_{n-1} = 0

 */
void fill(vec const & x, vec const & y, matrix & A, vec & b)
{
    std::size_t const n{x.size()-1};
    std::size_t const m{3*n};

    assert(A.size() == m*m);
    assert(b.size() == m);
    
    vec dx(n,0.0);

    for (std::size_t i{0}; i < n; ++i) {
	dx[i] = x[i+1]-x[i];
    }
    
    for (std::size_t i{0}; i < n; ++i) {
	std::size_t const I{i};
	std::size_t const J{3*i};

	A[I+m*J]     = dx[i]*dx[i];
	A[I+m*(J+1)] = dx[i];
	A[I+m*(J+2)] = 1.0;

	b[I] = (y[i+1]-y[i])/dx[i];
    }

    for (std::size_t i{0}; i < n-1; ++i) {
	std::size_t const I{n+i};
	std::size_t const J{3*i};

	A[I+m*J]     = 3.0*dx[i]*dx[i];
	A[I+m*(J+1)] = 2.0*dx[i];
	A[I+m*(J+2)] = 1.0;
	A[I+m*(J+5)] = -1.0;
    }
    

    for (std::size_t i{0}; i < n-1; ++i) {
	std::size_t const I{2*n-1+i};
	std::size_t const J{3*i};
	
	A[I+m*J]     = 3.0*dx[i];
	A[I+m*(J+1)] = 1.0;
	A[I+m*(J+4)] = -1.0;
    }

    A[(m-2)+m*1] = 1.0;    
    
    A[(m-1)+m*(m-1-2)] = 3.0*dx[n-1];
    A[(m-1)+m*(m-1-1)] = 1.0;
}


void spline_eval(vec const & x, vec const & y, vec & sol)
{
    std::size_t const n{x.size()-1};
    std::size_t const m{3*n};

    matrix A(m*m,0.0);
    vec b(m,0.0);
    
    fill(x,y,A,b);

    gauss(A,b,sol);

}


std::size_t locate(vec const & x, double v)
{
    std::size_t const n{x.size()-1};
    assert((v >= x[0]) && (v <= x[n]));

    std::size_t i{0};
    while (v > x[i]) {
	++i;
    }
    if (i!=0) {
	--i;
    }

    return i;
}

double spline(vec const & x, vec const & y, vec const & sol, double v)
{
    auto i = locate(x,v);

    auto const & ca = sol[3*i];
    auto const & cb = sol[3*i+1];
    auto const & cc = sol[3*i+2];

    auto const vx = v-x[i];
    return ((ca*vx+cb)*vx+cc)*vx+y[i];
}


int main()
{
    std::size_t n;
    std::ifstream in{"points.dat"};

    in >> n;
    if (!in) {
	return -1;
    }
    
    vec x(n), y(n);
    for (std::size_t i{0}; i < n; ++i) {
	in >> x[i] >> y[i];
    }
    --n;

    vec sol(3*n,0.0);
	
    spline_eval(x,y,sol);

    
    size_t constexpr m{15};
    
    std::cout.precision(16);

    double z = x[0];
    double const h = (x[n]-x[0]) / (m-1);

    for (size_t i=0; i < m;++i) {
	if (i==m-1) z = x[n];
	
	double w = spline(x,y,sol, z);

	std::cout << z << ' ' << w << '\n';
	z+=h;
    }
}
