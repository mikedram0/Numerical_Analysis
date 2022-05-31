#include <cmath>
#include <vector>
#include <iostream>


void difeq(double x, double y, std::vector<double> & dy)
{    
    using std::sin;
    using std::cos;
    
    dy[0] = cos(x)-sin(y) + x*x;
    
    dy[1] = 2.0 * x - sin(x) - cos(y) * dy[0];
    
    dy[2] = 2.0 - cos(x) - cos(y) * dy[1] + sin(y) * dy[0]*dy[0];
    
    dy[3] = sin(x) + 3.0 * sin(y) * dy[0] * dy[1] +
	cos(y) * (std::pow(dy[0],3) - dy[2]);	

    dy[4] = cos(x) + cos(y) * (6.0 * dy[0] * dy[0] * dy[1] - dy[3]) 
	+ (3.0 * dy[1] * dy[1] + dy[0] * (4.0 * dy[2] - std::pow(dy[0],3)))
	* sin(y);
}

void taylorstep(double x0, double y0, double x1, double & y1)
{
    std::vector<double> dy(5);
    
    difeq(x0, y0, dy);
    
    y1 = y0;
    double h = x1-x0;
    double term = h;
    
    for (std::size_t i = 0; i < dy.size(); ++i) {
	y1 += dy[i] * term;
	term *= h / (i+2);
    }

}


void taylor(double a, double b, double h, double ya, double & yb)
{
    double xold = a;
    double yold = ya;
    double ynew;

    while (xold < b) {
	double xnew = xold + h;
	if (xnew > b) {
	    xnew = b;
	}
	
	taylorstep(xold, yold, xnew, ynew);

	xold = xnew;
	yold = ynew;
    }

    yb = ynew;
}


int main()
{
    double xa = -1.0;    
    double xb = 1.0;
    double h = 0.01;
    double ya = 3.0;
    double yb;
    
    taylor(xa, xb, h, ya, yb);

    std::cout.precision(15);    
    std::cout << xb << ' ' << yb << '\n'; 
}
