#include <iostream>
#include <cmath>

// expseries 
int main()
{
    std::cout << u8"Δώσε πραγματικό αριθμό: ";
    double x;
    std::cin >> x;
    
    double sum{0.0};
    double a{1.0}; // oros gia n=0
    for (int n{1}; sum +a != sum; ++n) {  
	sum += a;
	
	//    next term is ...
	a *= x / n;
    }

    std::cout.precision(12);
    
    std::cout << "To athroisma sto " << x << " einai " << sum << '\n';
    std::cout << "H swsti timi: "<< std::exp(x) << '\n';
}
