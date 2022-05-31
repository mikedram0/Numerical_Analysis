#include <stdio.h>
#include <complex.h>
#include "../mylib.h"

typedef double complex COMPL;

COMPL f( COMPL x){
    return cpow(x,3) -5*cpow(x,2) + 13*x - 22;
}

int main(){

    printc(temn(1,2,f,1e-6));

    return 0;
}
