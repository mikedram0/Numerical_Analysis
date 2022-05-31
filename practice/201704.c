#include <stdio.h>
#include "../mylib.h"
#include "complex.h"

typedef double complex COMPL;

int main(){

	int N = 3;
	COMPL A[] = {1,3,4,2,5,-1,-2,1,3};
	COMPL roots[3];
	eigenval(N,A,roots,1.4+2.5I, 2+0.5I);
	//COMPLprint(N,1,roots);

	return 0;
}