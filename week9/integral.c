#include <stdio.h>
#include "../mylib.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

double trapezoid( double *f, double a, double b, int n){

    double result = 0;
    double h = (b-a)/n;
    int i =0;
    for(i = 1; i <= n-1 ; ++i)
        result += f[i];

    result = (h/2)*(f[0]+f[n]+2*result);
    return result;
}

double romberg( double *f, double a, double b, int n){
    
    assert(n % 4 == 0);
    double result1 = 0;
    double result2 = 0;
    double result3 = 0;
    double h = 4*(b-a)/n;
    int i =0;
    for(i = 4; i <= n-4 ; i+=4)
        result1 += f[i];
    for(i = 2; i <= n-2 ; i+=2)
        result2 += f[i];
    for(i = 1; i <= n-1 ; ++i)
        result3 += f[i];

    result1 = (h/2)*(f[0]+f[n]+2*result1);
    result2 = (h/4)*(f[0]+f[n]+2*result2);
    result3 = (h/8)*(f[0]+f[n]+2*result3);
    return (1.0/45)*(result1 - 20*result2 + 64*result3 );
}

double simpson( double *f, double a, double b, int n){

    assert(n  % 2 == 0);
    double h = (b-a)/n;
    double result1 = 0;
    double result2 = 0;
    int i =0;
    for(i = 1; i  <= n - 1 ; i+=2 )
        result1 += f[i];
    for(i = 2; i <= n -2 ; i+= 2)
        result2 += f[i];

    result1 = (h/3)*(f[0]+f[n]+ 4*result1 + 2*result2);
    return result1;
}

double quad( double *f, double a, double b, int n  ){
    if( n % 2 != 0 && n > 3 ){
        double h = (b-a)/n;
        double result1 = 0;
        double result2 = 0;
        int i =0;
        for(i = 4; i  <= n - 1 ; i+=2 )
            result1 += f[i];
        for(i = 5; i <= n -2 ; i+= 2)
            result2 += f[i];

        result1 = (h/3)*(f[3]+f[n]+ 4*result1 + 2*result2) + (3*h/8)*(f[0]+f[3]+3*(f[1]+f[2]));
    
        return result1;
    }

    if( n % 2 == 0 ){
        return simpson(f,a,b,n);
    }
}

void hermite( double roots[], double weights[], int n , double epsilon ){
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

        x = xn;
        roots[i] = x;
        roots[n-1-i] = -x;
        weights[i] = 2.0/(diff*diff);
        weights[n-1-i] = weights[i];
    }
    
}

double GH ( double (*f)( double ), int n ){
    double x[n];
    double w[n];
    double y[n];
    hermite(x,w,n,1e-10);
    int i;
    double integral = 0;
    for( i = 0; i < n; ++i ){
        integral += w[i]*f(x[i]);
    }
    return integral;

}



double f(double x){
    return sin(x);
}
double f2(double x){
    return x*x;
}
int main(){

    int i,j;
    int N;
    double trap = 0;
    double sim,rom;
    double sol = 2;
    double diff;
    for(i = 1; i <= 9 ; ++i){
        N = pow(2,i)*4;
        double *x = malloc( (N+1) * sizeof(double) );
        double *y = malloc( (N+1) * sizeof(double) );
        linspace(0,M_PI,N+1,x);
        printf("h = %g\n", x[1] - x[0]);
        for(j = 0; j <= N; ++j  ){
            y[j] = sin(x[j]);
        }
        double h = (M_PI)/N;
        trap = trapezoid(y,0,M_PI,N);
        sim = simpson(y,0,M_PI,N);
        rom = romberg(y,0,M_PI,N);
        printf("For N = %d, trapezoid : %.16g, exact : %g, diff: %.16g\n", N,trap,sol, fabs(trap-sol) );
        printf("For N = %d, simpson : %.16g, exact : %g, diff: %.16g\n", N,sim,sol, fabs(sim - sol));
        printf("For N = %d, romberg : %.16g, exact : %g, diff: %.16g\n", N,rom,sol, fabs(rom - sol));
        putchar('\n');
        free(x);
        free(y);
    }

    FILE *fptr = fopen("points.txt", "r");
    int N2;
    fscanf(fptr,"%d", &N2);
    double *y2 = malloc(N*sizeof(double));
    double *x2 = malloc(N*sizeof(double));
    for(i = 0; i < N2 ; ++i){
        fscanf(fptr, "%lf", &x2[i]);
        fscanf(fptr, "%lf", &y2[i]);
    } 
    double h = x2[1] - x2[0];
    printf("%g\n", h);
    double integral = quad(y2,x2[0],x2[N2-1],N2-1 );
    double diff2 = fabs(integral - 0.5*(1+ pow(M_E,3)*(sin(3)-cos(3) )   ));
    printf("%g, %g\n", integral,diff2);

    free(x2);
    free(y2);

    double integ = GH(f2,4);
    printf("%g, error: %g\n", integ, fabs(integ-sqrt(M_PI)/2 ) );

    return 0;
}
