#include "mylib.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

typedef double complex COMPL;

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

double falsepoint( double a , double b , double (*func)(double)  , double epsilon){
    double funca = func(a);
    double funcb = func(b);
    if( funca*funcb >= 0 ){
       fprintf(stderr,"Invalid initial parameters for bisec: f(a)f(b) >= 0\n") ;
       exit(-1);
    }
    int count = 0;
    double x = (a+b)/2.0;
    double f;
    while(fabs(f = func(x)) > epsilon){
        if( funca*f > 0){
            a = x;
	    funca = f;
       } else  {
            b = x;
	    funcb = f;
	}
	x =(b*funca - a*funcb)/(funca-funcb); 
        count++;
    }
    printf("%d\n",count);
    return x;
}
double illinois( double a , double b , double (*func)(double)  , double epsilon){

    double funca = func(a);
    double funcb = func(b);
    if( funca*funcb >= 0 ){
       fprintf(stderr,"Invalid initial parameters for bisec: f(a)f(b) >= 0\n") ;
       exit(-1);
    }
    int count = 0;
    int acount = 0;
    int bcount = 0;
    double x = (a+b)/2.0;
    double f;
    while(fabs(f = func(x)) > epsilon){
        if( funca*f > 0){
            a = x;
	    funca = f;
	    acount++;
            bcount = 0;
        } else { 
            b = x;
	    funcb = f;
            bcount++;
            acount = 0;
        }
        if( acount >= 2){
            x = (2*b*funca - a*funcb)/(2*funca - funcb);
            acount = 0;
        } else if (bcount >= 2){
            x = (b*funca - 2*a*funcb)/(funca - 2*funcb);
            bcount = 0;
        }else  
            x =(b*funca - a*funcb)/(funca-funcb); 
        count++;
    }
    printf("%d\n",count);
    return x;
}

double bisec( double a , double b , double (*func)(double)  , double epsilon){

    if( func(a)*func(b) >= 0 ){
       fprintf(stderr,"Invalid initial parameters for bisec: f(a)f(b) >= 0\n") ;
       exit(-1);
    }

    double x = (a+b)/2.0;
    int count = 0;
    double f;
    while(fabs(f = func(x)) > epsilon){
        if( func(a)*f > 0)
            a = x;
        else 
            b = x;

        x = (a+b)/2.0;
        count++;
    }
    printf("%d\n",count);
    return x;
}
COMPL temn( COMPL a, COMPL b, COMPL (*func)(COMPL), double epsilon){
	int count = 0;
	COMPL f;
	COMPL funca = func(a);
	COMPL funcb = func(b);
	COMPL c = (b*funca - a*funcb)/(funca - funcb);
	printc(c);
	while( cabs(f = func(c)) > epsilon){
		a = b;
		b = c;		
		funca = func(a);
		funcb = f;
		c = (b*funca - a*funcb)/(funca - funcb);
		printc(c);
		count++;
	}
	printf("%d\n",count);
	return c;
}

COMPL mod_temn( COMPL a, COMPL b, COMPL (*func)(COMPL,COMPL*,int), double epsilon, COMPL* roots, int n){
	int count = 0;
	COMPL f;
	COMPL funca = func(a,roots,n);
	COMPL funcb = func(b,roots,n);
	COMPL c = (b*funca - a*funcb)/(funca - funcb);
	while( cabs(f = func(c,roots,n)) > epsilon){
		a = b;
		b = c;		
		funca = func(a,roots,n);
		funcb = f;
		c = (b*funca - a*funcb)/(funca - funcb);
		count++;
	}
	printf("%d\n",count);
	return c;
}