#include <vector>
#include <cmath>
#include <iostream>
#include <functional>


double constexpr pi{3.14159265358979323846};

using vec = std::vector<double>;


void rk4(double x, vec const & y, double h,
	 std::function<double (int, double , vec const &)> f,
	 double & xnew, vec & ynew)
{
    auto const n = y.size();
    
    vec k1(n);
    vec k2(n);
    vec k3(n);
    vec k4(n);
    
    vec p(n);

    for (int i{0}; i < n; ++i) {
	k1[i] = h * f(i,x,y);
    }

    for (int i{0}; i < n; ++i) {
	p[i] = y[i] + k1[i]/2.0;
    }

    for (int i{0}; i < n; ++i) {	
 	k2[i] = h * f(i,x+h/2.0,p);

    }

    for (int i{0}; i < n; ++i) {
	p[i] = y[i] + k2[i]/2.0;
    }

    for (int i{0}; i < n; ++i) {	
 	k3[i] = h * f(i,x+h/2.0,p);
    }

    for (int i{0}; i < n; ++i) {
	p[i] = y[i] + k3[i];
    }
    
    for (int i{0}; i < n; ++i) {	
 	k4[i] = h * f(i,x+h,p);
    }

    for (int i{0}; i < n; ++i) {	
	ynew[i] = y[i] + (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;
    }
    
    xnew = x + h;
}



/*
  The system is 
  D theta / D t  = z
  D z / D t      = -sin(theta)

  theta -> y[0], z-> y[1].  D theta / D t -> f[0], D z / D t -> f[1]. t -> x
*/

double f(int i, double x, vec const & y)
{
    switch (i) {
    case 0:
	return y[1];
    case 1:
	return -std::sin(y[0]);
    default:
	std::cerr << "error in f\n";
	return -1.0;
    }
}


/*
  solution of theta'' = - theta:   theta = A cos(t) + B sin(t)
  theta(0) = 45 deg. , theta'(0) = 0 
*/
void approx(double t, vec & y)
{
    double constexpr a{45.0 * pi/180.0}; 

    y[0] =  a * std::cos(t);
    y[1] = -a * std::sin(t);
}


int main()
{
    double constexpr tinit{0.0};
    double constexpr tfin{10.0};
//    double constexpr tfin{200*pi};
    int constexpr niter{100};
    double constexpr h{(tfin-tinit)/niter};


    vec const yinit{45.0 / 180.0 * pi, 0.0};

    double told{tinit};
    vec yold{yinit};
    vec ynew(yinit.size());
    vec yapprox(yinit.size());

    
    std::cout.precision(15);

    for (int k{0}; k < niter; ++k) {
	double tnew;
	rk4(told,yold,h,f,tnew,ynew);
	told = tnew;
	for (int i{0}; i < yold.size(); ++i) {
	    yold[i] = ynew[i];
	}

	approx(tnew, yapprox);

	std::cout << tnew << ' ';
	for (auto & z : ynew) {
	    std::cout  << z << ' ';
	}
	for (auto & z : yapprox) {
	    std::cout  << z << ' ';
	}
	std::cout << '\n';
	
    }
}
