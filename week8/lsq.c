#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mylib.h"

double mean(double* x, int n){
    int i = 0;
    double sum = 0;
    for(i = 0; i < n ; ++i){
        sum += x[i];
    }
    sum /= n;
    return sum;
}

void lsq(double* x, double* y, int n, double *a, double *b){
    int i = 0;
    double xy = 0;
    double x2 = 0;
    double y2 = 0;
    double xm = mean(x,n);
    double ym = mean(y,n);


    for(i = 0; i < n; ++i){
        xy += x[i]*y[i];
        x2 += x[i]*x[i];
        y2 += y[i]*y[i];
    }
    xy /= n;
    x2 /= n;
    y2 /= n;

    *a = (xy - xm*ym)/(x2 - xm*xm);
    *b = (ym-(*a)*xm);

    double r2 = (*a)*(*a)*(x2 - xm*xm)/(y2 - ym*ym);
    printf("r^2 = %g\n", r2);
   
}

void fit( double *x, double *y, int n, double *a, double *b, char type){
    int i;
    double yy[n];
    double xx[n];
    switch(type){
        case 'l':
            lsq(x,y,n,a,b);
            break;

        case 'e':
            for(i = 0; i < n; ++i){
                xx[i] = exp(x[i]);
            }
            lsq(xx,y,n,a,b);
            break;

        case 'p':
            for(i = 0; i<n ; ++i){
                yy[i] = log(y[i]);
                xx[i] = log(x[i]);
            }
            lsq(xx,yy,n,a,b);
            *b = exp(*b);
            break;

        case 'g':
            for(i = 0; i<n; ++i){
                xx[i] = log(x[i]);
            }
            lsq(xx,y,n,a,b);
            break;
    
    }

}
void polynomial( double *x, double *y, double *res,int M, int N ){
    
    double *A = malloc((M+1)*(M+1)*sizeof(double));
    //double *B = malloc((M+1)*sizeof(double));

    int i,j,k;
    double sum1 = 0;
    double sum2 = 0;

    for(i = 0; i <= M ; ++i ){
        sum1 = 0;
        for(k = 0; k <= N-1 ; ++k){
            sum1 += pow(x[k],i)*y[k];
        }
        res[i] = sum1;
        for(j = 0; j <= M; ++j){
            sum2 = 0;
            for( k = 0; k <= N-1 ; ++k){
                sum2 += pow(x[k],i+j);
            }
            A[i + (M+1)*j] = sum2; 
        }
    }
    
    gauss((M+1),A,1,res,1e-14);
    

}
void linspace(double a, double b, int N, double *x){

    int i = 0;
    double h = (b-a)/(N-1);
    for( i = 0 ; i < N ; ++i ){
    
        x[i] = a + i*h;
    }
}
void range(double a, double b, double h, double *x){

    int i = 0;
    int N = (int) (b-a)/h + 1;
    for( i = 0 ; i < N ; ++i ){
    
        x[i] = a + i*h;
    }
}

int main(){
    
    FILE *fpt = fopen("points.txt","r");
    if(fpt == NULL){
       printf("Could not open file, exiting");
       return -1; 
    }
    int N;
    fscanf(fpt, "%d", &N);
    double *x = malloc(N*sizeof(double));
    double *y = malloc(N*sizeof(double));
    int i;
    for( i = 0 ; i < N ; ++i){
        fscanf(fpt,"%lf", &x[i]); 
        fscanf(fpt,"%lf", &y[i]); 
        }

    fclose(fpt);

    double a;
    double b;
    fit(x,y,N,&a,&b,'l');
    printf("%g*x+%g\n",a,b);
    fit(x,y,N,&a,&b,'p');
    printf("%g*x^%g\n",b,a);
    fit(x,y,N,&a,&b,'e');
    printf("%g*e^x+%g\n",a,b);
    fit(x,y,N,&a,&b,'g');
    printf("%g*ln(x)+%g\n",a,b);

    double T[21];
    double P[] = {0.0013,0.0162,0.0297,0.0318,0.0484,0.0965,0.1357,0.2947,0.4563,0.5398,0.8884,1.0031,1.4193,1.9052,2.4026,2.5031,3.9072,4.3156,5.5060,6.9044,7.6370};

    range(300,2300,100,T);
    fit(T,P,21,&a,&b,'p');
    printf("%g*x^%g\n",b,a);

    double x2[] = {0,0.25,0.5,0.75,1};
    double y2[] = {1,1.284,1.6487,2.117,2.7183};

    double res2[3];

    polynomial(x2,y2,res2,2,5);
    print(1,3,res2);
    return 0;
}
