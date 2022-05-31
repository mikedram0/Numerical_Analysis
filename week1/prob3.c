#include <stdio.h>

int main(){
	int i =1;
	double sum = 0;
	for(i = 1 ; i <= 19 ; ++i)
		sum += i/10.0;
	printf("%g", sum);
}
