#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double testf(double x){
    return pow(x,3) - 5;
}

double f1( double x){
    return (pow(x,3) + 4*pow(x,2) - 10);
}

double f2( double x ){
    return (sqrt(x) - cos(x));
}
double f3 (double x){
    return -2 + 6.2*x - 4.0*pow(x,2)+0.7*pow(x,3);
}   
double f4 (double x, double* roots, int n){
	int i = 0;
	double result = pow(x,10) - 0.95;
	for(i = 0 ; i < n ; ++i)
		result /= (x-roots[i]);

	return result;
}
double f5 ( double x){
	return pow(x,2)-pow(1-x,5);
}
double f6 (double x){
	return 3*log(x) + 5;
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
double temn( double a, double b, double (*func)(double), double epsilon){
	int count = 0;
	double f;
	double funca = func(a);
	double funcb = func(b);
	double c = (b*funca - a*funcb)/(funca - funcb);
	while( fabs(f = func(c)) > epsilon){
		a = b;
		b = c;		
		funca = func(a);
		funcb = f;
		c = (b*funca - a*funcb)/(funca - funcb);
		count++;
	}
	printf("%d\n",count);
	return c;
}
double mod_temn( double a, double b, double (*func)(double,double*,int), double epsilon, double* roots, int n){
	int count = 0;
	double f;
	double funca = func(a,roots,n);
	double funcb = func(b,roots,n);
	double c = (b*funca - a*funcb)/(funca - funcb);
	while( fabs(f = func(c,roots,n)) > epsilon){
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
int main(){

//   printf("%.15g\n",bisec(1,2,f1,1e-15));
 //  printf("%.15g\n",bisec(0,1,f2,1e-15));

//   printf("%.15g\n",falsepoint(0.4,0.6,f3,1e-15));
//   printf("%g\n",falsepoint(0,1.4,f4,1e-6) );
//   printf("%g\n",bisec(0,1.4,f4,1e-6) );
//   putchar('\n');
 //  printf("%g\n",bisec(0,1,f5,1e-9) );
//   printf("%g\n",falsepoint(0,1,f5,1e-9) );
//   printf("%g\n",illinois(0,1,f5,1e-9) );
//   printf("%g\n",temn(0.1,0.2,f6,1e-6));
    printf("%g\n",temn(10,2,testf,1e-20)); 

   int n = 2;
   int i = 0;
   double roots[n]; 
   for(i = 0; i < n ; ++i){
 	roots[i] = mod_temn(0,1,f4,1e-6,roots,i);        
	printf("%g",roots[i]);
   }
 
    return 0;
}
