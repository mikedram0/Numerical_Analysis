#include <iostream>
#include <cmath>
#include <functional>


double f(double x)
{
    return x*x -std::pow(1.0-x,5);
}



int bisection(double a, double b, //     Initial limits
	      std::function<double (double)> f, // function
	      double toler, // tolerance
	      double & x, // final approximation
	      int & iter, // Number of iterations needed
	      int & eval // Number of function evaluations
    )
{
    double fa{f(a)};
    double fb{f(b)};

    eval = 2;
    
    if (fa * fb > 0.0) {
	return -1;
    }

    for (iter = 1; ; ++iter) {
	x = (a+b)/2.0;
	double fnval{f(x)};
	++eval;
	
	//  Check if root is found
	if (std::abs(fnval) < toler) {
	    break;
	}

	if (fa * fnval > 0.0) {
	    a = x;
	    fa = fnval;
	} else {
	    b = x;
	    fb = fnval;
	}	
    }

    return 0;
}




//     Regula falsi (false position) method

int regfal(double a, double b, //     Initial limits
       std::function<double (double)> f, // function
       double toler, // tolerance
       double & x, // final approximation
       int & iter, // Number of iterations needed
       int & eval // Number of function evaluations
    )
{
    double fa{f(a)}, fb{f(b)};
    eval = 2;
    
    if (fa * fb >= 0.0) {
	std::cerr << "limits are wrong\n";
	return -1;
    }
    
    for (iter = 1; ; ++iter) {
  	x = a - fa * (b-a) / (fb-fa);
	double fx{f(x)};
	++eval;
	
	//     check if root is found
	if (std::abs(fx) < toler) {
	    break;
	}

	//     new limits
	if (fa * fx > 0.0) {
	    a = x;
	    fa = fx;
	} else {
	    b = x;
	    fb = fx;
	}
    }

    return 0;
}

//  Regula falsi (false position) method; illinois algorithm.
int illinois(double a, double b, //     Initial limits
	     std::function<double (double)> f, // function
	     double toler, // tolerance
	     double & x, // final approximation
	     int & iter, // Number of iterations needed
	     int & eval // Number of function evaluations
    )
{
    double fa{f(a)}, fb{f(b)};
    eval = 2;

    if (fa * fb >= 0.0) {
	std::cerr << "limits are wrong\n";
	return -1;
    }
    

    int chnga{0};
    int chngb{0};
    
    for (iter = 1; ; ++iter) {
	if (chnga < 2 && chngb < 2) {
	    x = a - fa * (b-a) / (fb-fa);
	}

	if (chnga == 2) {
	    chnga = 0;
	    x = a - 2*fa * (b-a) / (fb-2*fa);
	}
	
	if (chngb == 2) {
	    chngb = 0;
	    x = a - fa * (b-a) / (2*fb-fa);
	}

	
	double fx{f(x)};
	++eval;
	
	//     check if root is found
	if (std::abs(fx) < toler) {
	    break;
	}

	//     new limits
	if (fa * fx > 0.0) {
	    a = x;
	    fa = fx;

	    chngb = 0;
	    ++chnga;
	} else {
	    b = x;
	    fb = fx;

	    chnga = 0;
	    ++chngb;
	}
    }

    return 0;
}







int main()
{
    double x;
    double constexpr toler = 1e-6;
    int iter;
    int eval;

    int err = bisection(0.0, 1.0, f, toler, x, iter, eval);
    if (err != 0) {
	std::cerr << "error in bisection\n";
	return -1;
    }

    std::cout << "x is " << x << " by bisection method,"
	      << " in " << iter << " iterations "
	      << "and "<< eval << " function evaluations\n";


    err = regfal(0.0, 1.0, f, toler, x, iter, eval);
    if (err != 0) {
	std::cerr << "error in regfal\n";
	return -1;
    }
    std::cout << "x is " << x << " by false position method,"
	      << " in " << iter << " iterations "
	      << "and "<< eval << " function evaluations\n";



    err = illinois(0.0, 1.0, f, toler, x, iter, eval);
    if (err != 0) {
	std::cerr << "error in illinois\n";
	return -1;
    }
    std::cout << "x is " << x << " by illinois method,"
	      << " in " << iter << " iterations "
	      << "and "<< eval << " function evaluations\n";

}

