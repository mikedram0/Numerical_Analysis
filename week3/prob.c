#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "../mylib.h"
//#define printc(n) printf("%g%+gi\n", creal(n),cimag(n))

typedef double complex COMPL;

void printc( COMPL x){
    printf("%g%+gi\n", creal(x),cimag(x));
}

COMPL f0(COMPL x){
	return cpow(x,2) +1;
}

COMPL f1( COMPL x){
    return csin(x) - cpow(x,2);
}

double f2(double x){
	return (pow(x,2)+5)/6;
}
double f3(double x){
	return sqrt(6*x -5);
}
double f4(double x){
	return acos(pow(x,1/3.0));
}
double f5(double x){
	return sin(x) - pow(x,2);
}
double df5(double x){
	return cos(x) - 2*x;
}
double f6(double x){
	return 3*x*exp(x)-1;
}
double df6(double x){
	return 3*(exp(x) + x*exp(x));
}
double f7(double x){
	return 4*cos(x) - exp(-x);
}
double df7(double x){
	return -4*sin(x) + exp(-x);
}
double gf7(double x){
	return acos(exp(-x)/4);
}
COMPL muller( COMPL x0, COMPL x1, COMPL x2, COMPL (*func)( COMPL),double epsilon){
    
    int counter = 0;
    COMPL x,w0,w1,a,b,c,f0,f1,f2,D,d;
 
    while(cabs(f2 = func(x2)) > epsilon ){

        f1 = func(x1);
        f0 = func(x0);
        w0 = (f2-f0)/(x2-x0);
        w1 = (f2-f1)/(x2-x1);

        a = (w1-w0)/(x1-x0);
        b = w0 + a*(x2-x0);
        c = f2;
   
        D = csqrt(cpow(b,2)-4*a*c);

        if (cabs(b+D) > cabs(b-D))
           d = b-D;
        else
           d = b+D;

        x = x2 - 2*c/d;  

        x0 = x1;
        x1 = x2;
        x2 = x;
        
        ++counter;
    }
    
    printf("%d\n",counter);
    return x2;

}

double stablepoint(double a, double  (*func)(double), double epsilon){

	double g;
	int count = 0;
	while( ( fabs(a - (g = func(a)) ) ) > epsilon  ){
		a = g;
		count++;
}
	printf("%d\n", count);
	return a;

}

double NR(double a, double (*func)( double), double (*dfunc)(double), double  epsilon){
	int count = 0;
	double f;
	while( (fabs(f = func(a))) > epsilon ){
		a = a - f/(dfunc(a));
		count++;
}
	printf("%d\n",count);
	return a;	
}
int main(){

    printc(muller(1,2,3,f1,1e-15));
    printc(muller(-1-2I,0.4+4I,1+0I,f0,1e-15));
	printf("%g\n", stablepoint(2,f2,1e-15));
	printf("%g\n", stablepoint(10,f3,1e-15));
	printf("%g\n", stablepoint(0.6,f4,1e-15));
	printf("%g\n", NR(2,f5,df5,1e-15));
	printf("%g\n", NR(2,f6,df6,1e-15));
	putchar('\n');
	printf("%g\n", temn(1,2,f7,1e-8));
	printf("%g\n", bisec(1,2,f7,1e-8));
	printf("%g\n", NR(1,f7,df7,1e-8));
	printf("%g\n", stablepoint(1,gf7,1e-8));
    return 0;

}
