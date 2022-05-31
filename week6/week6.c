#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <complex.h>
#include "../mylib.h"

typedef double complex CMPL ;

CMPL f( CMPL z){
    double x = creal(z);
    double y = cimag(z);
    return  4*pow(x,2) - pow(y,3) + 28 + I*( 3*pow(x,3) + 4*pow(y,2) - 145 );

}
void printc( CMPL  x){
    printf("%g%+gi\n", creal(x),cimag(x));
}

double f1( double* a, int i){
    double x = a[0];
    double y = a[1];
    double z = a[2];
    switch(i){
        case(0):
            return (x + y + z - 3);
            break;
        case(1):
            return 1;
            break;
        case(2):
            return 1;
            break;
        case(3):
            return 1;
            break;
    }
}
double f2( double* a , int i){
    double x = a[0];
    double y = a[1];
    double z = a[2];
    switch(i){
        case(0):
            return (pow(x,2)*y + pow(y,2)*z + pow(z,2)*x - 4);
            break;
        case(1):
            return (2*x*y + z*z);
            break;
        case(2):
            return (x*x + 2*y*z);
            break;
        case(3):
            return ( y*y + 2*z*x);
            break;
    }
}
double f3( double* a, int i){
    double x = a[0];
    double y = a[1];
    double z = a[2];
    switch(i){
        case(0):    
            return (x*x + y*y + z*z - 5);
            break;
        case(1):
            return (2*x);
            break;
        case(2):
            return (2*y);
            break;
        case(3):
            return (2*z);
    }
}

double (*F[])( double*, int) = {f1,f2,f3};
double absmax( int M, int N, double A[]){
    int i = 0;
    double max = 0;
    double temp;
    for(i = 0; i<= M*N - 1; ++i){
        if( (temp = fabs(A[i]))  > max)
            max = temp;
    }
    return max;
}


void matdiff( int M, int N, double A[], double B[]){
    int i = 0;
    for(i == 0; i<= N*M - 1 ; ++i){
        A[i] -= B[i];
    }

}
int nonlinear(double *xi,  double (*F[])(double*, int), int N, double epsilon){
    double A[N*N];
    double b[N];
    int i,j;
    int k = 0;

    for( i = 0; i< N; ++i){
        b[i] = F[i](xi,0);
    }
    while( absmax(1,N,b) > epsilon){
        for(i = 0; i < N; ++i){
            for(j = 0; j < N; ++j){
                A[i+N*j] = F[i](xi,j+1);
            }
        }
        gauss(N,A,1,b);
        matdiff(1,N,xi,b);
        for( i = 0; i< N; ++i){
            b[i] = F[i](xi,0);
        }
        k++;
    }
    printf("%d\n",k);
}

int main(){
    
    printc(temn(-1+1I,2+3I,f,1e-12 ));
    double xi[] = {1,2,3};
    nonlinear(xi, F,3,1e-10);
    print(1,3,xi);

    return 0;
}
