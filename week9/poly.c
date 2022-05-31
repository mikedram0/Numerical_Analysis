#include <stdio.h>
#include "../mylib.h"
#include <assert.h>
#include <math.h>

double MOD_temn( double a, double b, double (*func)(double,int,double*,int), double epsilon, double* roots, int n, int m){
	double f;
	double funca = func(a,n,roots,m);
	double funcb = func(b,n,roots,m);
	double c = (b*funca - a*funcb)/(funca - funcb);
	while( fabs(c - b) > epsilon){
		a = b;
		b = c;		
		funca = funcb;
		funcb = func(b,n,roots,m);
		c = (b*funca - a*funcb)/(funca - funcb);
	}
	return c;
}

double hermite( double x, int n, double* roots, int m){
    assert( n > 0 );
    //Using re-normalized formulas to avoid overflow. See Numeral recipies p186 for details.
    double p0 = 0;
    double p1 = 1/(pow(M_PI,0.25));
    double p2;
    int i;
    for(i = 0; i < n; ++i){
        p2 = x*sqrt(2.0/(i+1.0))*p1 - sqrt( (double) i/(i+1.0))*p0    ; 
        p0 = p1;
        p1 = p2;
        printf("%g\n", p2 );
    }
    for(i = 0 ; i < m ; ++i)
		p2 /= (x-roots[i]);


    return p2;
}    
double legendre( double x, int n, double* roots, int m){
    assert( n > 0 );
    double p0 = 1;
    double p1 = x;
    double p2 = x;
    int i;
    for(i = 1; i < n; ++i){
        p2 = (1.0/(i + 1)) * ( (2*i+1)*x*p1 - i*p0  );
        p0 = p1;
        p1 = p2;
    }

    for(i = 0 ; i < m ; ++i)
		p2 /= (x-roots[i]);


    return p2;
}    

void GH( double roots[], double weights[], int n , double epsilon ){
    double invpi4 = pow(M_PI,-0.25);
    double xn = 0;
    double x,p0,p1,p2,diff;
    int i,j;
    int k = 0;
    int ITERMAX = 100;
    int m = (n+1)/2;
    for(i = 0; i < m; ++i){
        switch(i){
            case 0:
                x = sqrt( (2.0*n+1)) -1.85575*pow((2.0*n+1.0),-0.16667 );
                break;
            case 1:
                x -= 1.14*pow( (double) n, 0.426)/x;
                break;
            case 2:
                x = 1.86*x-0.86*roots[0];
                break;
            case 3:
                x = 1.91*x-0.91*roots[1];
                break;
            default:
                x = 2.0*x-roots[i-2];
                break;
        }
        //printf("%g\n",x);
        k = 0;
        while( ++k < ITERMAX ){
            p0 = 0;
            p1 = invpi4;
            for(j = 0; j < n; ++j){
                p2 = x*sqrt(2.0/(j+1.0))*p1 - sqrt( (double) j/(j+1.0))*p0    ; 
                p0 = p1;
                p1 = p2;
            }
            diff = sqrt(2.0*n)*p0;
            xn = x - p2/diff;
            if( fabs(xn-x) < epsilon ){break;}
            x = xn;
        }    

        printf("%d\n", k);
        x = xn;
        roots[i] = x;
        roots[n-1-i] = -x;
        weights[i] = 2.0/(diff*diff);
        weights[n-1-i] = weights[i];
    }
    
}



int main(){
    
  // int i = 0;
   int n = 3;
   //double roots[n];
   //for(i = 0; i < n; ++i){
   //     roots[i] = MOD_temn(1.323,2.712,hermite,1e-6,roots,n,i);
   //} 
   // print(n,1,roots);
   // return 0;
    double roots[n];
    double weights[n];
    GH( roots, weights, n, 1e-6);
    print(n,1,roots);
    print(n,1,weights);
}

