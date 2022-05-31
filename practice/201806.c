#include <stdio.h>
#include "../mylib.h"

double f1( double* a, int i){
    double x1 = a[0];
    double x2 = a[1];
    double x3 = a[2];
	double x4 = a[3];
    switch(i){
        case(0):    
            return (x1  + x2 - 2);
            break;
        case(1):
            return (1);
            break;
        case(2):
            return (1);
            break;
        case(3):
            return (0);
        case(4):
            return (0);
    }
}
double f2( double* a, int i){
    double x1 = a[0];
    double x2 = a[1];
    double x3 = a[2];
	double x4 = a[3];
    switch(i){
        case(0):    
            return (x1*x3  + x2*x4);
            break;
        case(1):
            return (x3);
            break;
        case(2):
            return (x4);
            break;
        case(3):
            return (x1);
        case(4):
            return (x2);
    }
}
double f3( double* a, int i){
    double x1 = a[0];
    double x2 = a[1];
    double x3 = a[2];
	double x4 = a[3];
    switch(i){
        case(0):    
            return (x1*x3*x3  + x2*x4*x4 - 2.0/3.0);
            break;
        case(1):
            return (x3*x3);
            break;
        case(2):
            return (x4*x4);
            break;
        case(3):
            return (2*x1*x3);
        case(4):
            return (2*x2*x4);
    }
}
double f4( double* a, int i){
    double x1 = a[0];
    double x2 = a[1];
    double x3 = a[2];
	double x4 = a[3];
    switch(i){
        case(0):    
            return (x1*x3*x3*x3  + x2*x4*x4*x4);
            break;
        case(1):
            return (x3*x3*x3);
            break;
        case(2):
            return (x4*x4*x4);
            break;
        case(3):
            return (3*x1*x3*x3);
        case(4):
            return (3*x2*x4*x4);
    }
}


int main(){
	double (*F[])( double*, int) = {f1,f2,f3,f4};
	double xi[] = {1,2,3,4};
	nonlinear(xi, F,3,1e-10);
	print(4,1,xi);
	return 0;
}