#include <stdio.h>
#include "../mylib.h"
#include <complex.h>
#define N 3

int main(){
    typedef double complex COMPL;

    double A[N*N] = {7,1,-2,1,4,-4,-2,-4,2} ;
    double B[N] = {9.0/2.0, -2.0/3.0, 3.0};

    gauss(N,A,1,B,1e-14);
    print(1,N,B);
    
    COMPL A2[N*N] = {7,1,-2,1,4,-4,-2,-4,2};
    COMPL roots[3];
    eigenval(N,A2,roots,1,2);
    putchar('\n');
    COMPLprint(1,N,roots);

    return 0;
}
