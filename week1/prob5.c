#include <stdio.h>

double exp( double x){
	double result = 1;
	double term = 1;
	int i = 1;
	while(result + term > result){
		term = term * x / i;
		i++ ;
		result += term;
}

	return result;
}



int main(){

	printf("%.15g", exp(1.23));

	return 0;
}
