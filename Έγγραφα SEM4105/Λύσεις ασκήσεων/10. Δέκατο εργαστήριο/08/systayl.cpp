#include <iostream>
#include <cmath>
#include <vector>

using vec = std::vector<double>;


void difeq(double x, double y, double z, vec & dy, vec & dz)
{
    using std::sin;
    using std::cos;

    dy[0] = y + z*z - x*x*x;
    dz[0] = y*y*y + z + cos(x);

    dy[1] =-3.0 * x*x + dy[0] + 2.0 * z * dz[0];
			 
    dz[1] = 3.0 * y*y * dy[0] + dz[0] - sin(x);
 
    dy[2] = -6.0 * x + dy[1] + 2.0 * (dz[0]*dz[0] + z * dz[1]);
    dz[2] = 3.0 * y * (2.0 * dy[0]*dy[0] + y * dy[1]) + dz[1] - cos(x);

    dy[3] = -6.0 + dy[2] + 2.0 * (3.0 * dz[0] * dz[1] + z * dz[2]);
    dz[3] = dz[2] + sin(x) + 6.0 * dy[0] * dy[0] * dy[0]
	+ 3.0 * y * (6.0 * dy[0] * dy[1] + y * dy[2]);
}


void taylorstep(double x0, double y0, double z0, double x1,
		double & y1, double & z1)
{
    vec dy(5), dz(5);
    
    difeq(x0, y0, z0, dy, dz);

    y1 = y0;
    z1 = z0;
    double h = x1-x0;
    double term = h;
    
    for (std::size_t i = 0; i < dy.size(); ++i) {
	y1 += dy[i] * term;
	z1 += dz[i] * term;
	term *= h / (i+2);
    }
}

void taylor(double a, double b, double h, double ya, double za)
{
    double xold = a;
    double yold = ya;
    double zold = za;
    double ynew, znew;
    vec dy(4), dz(4);

    std::cout.precision(15);
    
    while (xold < b) {
	double xnew = xold + h;
	if (xnew > b) {
	    xnew = b;
	}
	
	taylorstep(xold, yold, zold, xnew, ynew, znew);

	xold = xnew;
	yold = ynew;
	zold = znew;

	std::cout << xnew << ' '<< ynew << ' ' << znew << '\n';
    }
    
}


int main()
{
    double constexpr a{0.0};
    double constexpr b{1.0};
    double constexpr h{0.01};
    
    double constexpr yinit{0.3};
    double constexpr zinit{0.1};

    taylor(a,b,h,yinit,zinit);
}
