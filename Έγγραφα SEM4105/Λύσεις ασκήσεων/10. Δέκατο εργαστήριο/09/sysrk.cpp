#include <vector>
#include <fstream>
#include <functional>

using vec = std::vector<double>;


double force(double x)
{
    return x  - 0.01*x*x*x;
}



/*
  m x'' = F(x) =>    x' = v     =>    y0 = x, y1 = v
                     v' = F/m
*/
void f(double t, vec const & y, vec & dy)
{
    double constexpr mass{2.0};
    
    dy[0] = y[1];
    dy[1] = force(y[0]) / mass;
}


void heun(double t0, vec const & y0, double t1, vec & y1, 
	  std::function<void (double, vec const &, vec &)> f)    
{
    auto h = t1-t0;

    auto const n = y0.size();

    vec k1(n), k2(n);

    f(t0, y0, k1);
    for (auto & z : k1) {
	z *= h;
    }

    vec p(y0);

    for (vec::size_type i{0}; i < n; ++i) {
	p[i] += k1[i];
    }

    f(t1, p, k2);
    for (auto & z : k2) {
	z *= h;
    }

    for (vec::size_type i{0}; i < n; ++i) {
	y1[i] = y0[i] + (k1[i] + k2[i])/2.0;
    }
}

void ralston(double t0, vec const & y0, double t1, vec & y1,
	     std::function<void (double, vec const &, vec &)> f)
{
    auto h = t1-t0;

    auto const n = y0.size();

    vec k1(n), k2(n);

    f(t0, y0, k1);
    for (auto & z : k1) {
	z *= h;
    }

    vec p(y0);

    for (vec::size_type i{0}; i < n; ++i) {
	p[i] += 2.0/3.0 * k1[i];
    }

    f(t0+2.0/3.0*h, p, k2);
    for (auto & z : k2) {
	z *= h;
    }

    for (vec::size_type i{0}; i < n; ++i) {
	y1[i] = y0[i] + (k1[i] + 3.0 * k2[i])/4.0;
    }
}



int main()
{
    double t0{0.0};
    double x0{2.5e-2};
    double v0{0.0};
    
    std::size_t constexpr nsteps{100000};
    double constexpr h{0.001};
    
    vec y0{x0,v0};
    vec y1(y0.size());

    std::ofstream out{"sysrk.txt"};
    for (std::size_t i{0}; i <= nsteps; ++i) {
	out << t0 << ' ' << y0[0]  << ' ' << y0[1] << '\n';
	
	//heun(t0, y0, t0+h, y1, f);
	ralston(t0, y0, t0+h, y1, f);
	
	t0 += h;
	y0 = y1;
    }    
}
