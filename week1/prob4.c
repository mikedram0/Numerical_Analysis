#include <stdio.h>
#include <float.h>

int main(){

	double e = 1;
	while(1+e != 1)
		e /= 2;
	e = e*2; // We want the last value BEFORE 1 + e = 1
	printf("Algorithm result: %.15g \n",e);
	printf("Machine epsilon: %g \n", DBL_EPSILON);
	


	return 0;
}
