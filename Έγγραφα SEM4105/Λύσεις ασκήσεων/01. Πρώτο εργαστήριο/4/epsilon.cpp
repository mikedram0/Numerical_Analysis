#include <iostream>
#include <limits>


int main()
{
    float epsre{1.0f};
    
    while (1.0f + epsre != 1.0f) {
	epsre *= 0.5f;
    }
  
    //    here epsre is half the epsilon.
    std::cout << "float epsilon is " <<  2.0f * epsre << '\n';
    std::cout << "float epsilon from C++ is "
	      << std::numeric_limits<float>::epsilon() << '\n';
  
  double epsdo = 1.0;
  
  while (1.0 + epsdo != 1.0) {
    epsdo *= 0.5;
  }
  
  //     here epsdo is half the epsilon.
  std::cout << "double epsilon is " <<  2.0 * epsdo << '\n';
  
  std::cout << "double epsilon from C++ is "
	    << std::numeric_limits<double>::epsilon() << '\n';

}
