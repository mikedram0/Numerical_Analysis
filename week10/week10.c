#include <stdio.h>
#include "../mylib.h"
#include <math.h>

void feuler( double (*f)(double, double), double x0, double y0, double x1, double h, double*out  ){

    int n = (x1 - x0)/h;
    int i = 0;
    out[0] = y0;
    double y = y0;
    double x = x0;
    for(i = 0; i < n ; ++i){
        y += h*f(x,y); 
        x += h;
        out[i+1] = y;
    }


void beuler( double (*f)(double, double), double x0, double y0, double x1, double h, double*out  ){
    
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
	}
	return c;


}
double f( double x, double y){
    return cos(x) - x*sin(x) + 0*y;
}

int main(){
    int i = 0;
    double h = 0.01;
    double x0 = 0;
    double x1 = 3;

    int n = (x1-x0)/h;
    printf("%d\n",n);
    double out[n + 1];
    
    feuler(f, x0,2,x1,h,out);
    //print(1,n,out);

    FILE *fptr = fopen("data.txt", "w");
    for( i = 0; i <= n ; ++i){
        fprintf(fptr,"%g %g\n", x0+i*h,out[i]  );
        }
        
    fclose(fptr);


    return 0;
}
